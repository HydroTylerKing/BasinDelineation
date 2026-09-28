#!/usr/bin/env python3
"""
Basin Delineator - click a map, get the watershed and its river network.

First run downloads the Copernicus GLO-30 DEM for Lesotho (plus a buffer so basins that
cross into South Africa are complete), conditions it once, and caches the result.
Every launch after that opens in a couple of seconds.

    python basin_app.py                    # Lesotho, 30 m Copernicus DEM (default)
    python basin_app.py --res 90           # lighter 90 m version (low-RAM machines)
    python basin_app.py --dem my_dem.tif   # use your own DEM (any CRS)
    python basin_app.py --bbox 27 -30 28.5 -29   # custom area (W S E N, degrees)
"""
import argparse
import json
import math
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import hydro  # noqa: E402

# Lesotho plus buffer: Caledon (Mohokare) headwaters on the Free State side, Tele river, escarpment.
LESOTHO_BBOX = (26.9, -30.85, 29.6, -28.15)
R_EARTH = 6371008.8
DEG = math.pi / 180.0


def log(msg):
    print(msg, flush=True)


# --------------------------------------------------------------------------------------
# DEM acquisition
# --------------------------------------------------------------------------------------
def _tile_label(lat, lon):
    ns = ("S" if lat < 0 else "N") + f"{abs(lat):02d}_00"
    ew = ("W" if lon < 0 else "E") + f"{abs(lon):03d}_00"
    return ns, ew


def _ssl_context():
    import ssl
    try:  # python.org builds on macOS often lack system certificates; certifi comes with rasterio
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def download_copernicus(bbox, res, tile_dir):
    """Download Copernicus DEM COG tiles (AWS open data, no login) covering bbox."""
    w, s, e, n = bbox
    code = 10 if res == 30 else 30
    bucket = f"copernicus-dem-{res}m"
    tile_dir.mkdir(parents=True, exist_ok=True)
    files = []
    lats = range(math.floor(s), math.ceil(n))
    lons = range(math.floor(w), math.ceil(e))
    total = len(lats) * len(lons)
    i = 0
    for lat in lats:
        for lon in lons:
            i += 1
            ns, ew = _tile_label(lat, lon)
            name = f"Copernicus_DSM_COG_{code}_{ns}_{ew}_DEM"
            out = tile_dir / f"{name}.tif"
            if out.exists() and out.stat().st_size > 0:
                files.append(out)
                continue
            url = f"https://{bucket}.s3.amazonaws.com/{name}/{name}.tif"
            log(f"  downloading tile {i}/{total}: {name}")
            try:
                tmp = out.with_suffix(".part")
                with urllib.request.urlopen(url, timeout=120, context=_ssl_context()) as r, open(tmp, "wb") as f:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        f.write(chunk)
                tmp.replace(out)
                files.append(out)
            except urllib.error.HTTPError as ex:
                if ex.code in (403, 404):
                    log("    (no tile - ocean or outside coverage)")
                else:
                    raise
    if not files:
        raise SystemExit("No DEM tiles could be downloaded for that area.")
    return files


# --------------------------------------------------------------------------------------
# Terrain model
# --------------------------------------------------------------------------------------
class Terrain:
    """Conditioned DEM held in memory (memory-mapped from the cache)."""

    def __init__(self, cache_dir):
        meta = json.loads((cache_dir / "meta.json").read_text())
        self.meta = meta
        self.fdir = np.load(cache_dir / "fdir.npy", mmap_mode="r")
        self.acc = np.load(cache_dir / "acc.npy", mmap_mode="r")
        self.dem = np.load(cache_dir / "dem.npy", mmap_mode="r")
        # numba wants real arrays; loading into RAM once keeps clicks fast
        self.fdir = np.ascontiguousarray(self.fdir)
        self.acc = np.ascontiguousarray(self.acc)
        self.nr, self.nc = self.fdir.shape
        self.a, self.b, self.c0, self.d, self.e, self.f0 = meta["transform"]
        self.geographic = meta["geographic"]
        self.crs = meta["crs"]
        self.dy_m = meta["dy_m"]
        self.cell_km2_row = np.array(meta["cell_km2_row"], dtype=np.float64) if len(meta["cell_km2_row"]) > 1 \
            else np.full(self.nr, meta["cell_km2_row"][0])
        self.buf = None
        self.lock = threading.Lock()
        self.basins = {}
        self.next_id = 1
        self._to_ll = self._from_ll = None
        if not self.geographic:
            from pyproj import Transformer
            self._to_ll = Transformer.from_crs(self.crs, "EPSG:4326", always_xy=True)
            self._from_ll = Transformer.from_crs("EPSG:4326", self.crs, always_xy=True)

    # coordinate helpers --------------------------------------------------------------
    def rc_to_lonlat(self, rows, cols):
        x = self.c0 + self.a * (np.asarray(cols) + 0.5)
        y = self.f0 + self.e * (np.asarray(rows) + 0.5)
        if self._to_ll is not None:
            x, y = self._to_ll.transform(x, y)
        return np.asarray(x), np.asarray(y)

    def xy_to_lonlat(self, x, y):
        if self._to_ll is not None:
            x, y = self._to_ll.transform(x, y)
        return np.asarray(x), np.asarray(y)

    def lonlat_to_rc(self, lon, lat):
        x, y = lon, lat
        if self._from_ll is not None:
            x, y = self._from_ll.transform(lon, lat)
        col = int(math.floor((x - self.c0) / self.a))
        row = int(math.floor((y - self.f0) / self.e))
        return row, col

    def bounds_ll(self):
        rows = np.array([0, 0, self.nr, self.nr]) - 0.5
        cols = np.array([0, self.nc, 0, self.nc]) - 0.5
        lon, lat = self.rc_to_lonlat(rows, cols)
        return [float(lon.min()), float(lat.min()), float(lon.max()), float(lat.max())]

    # queries ----------------------------------------------------------------------------
    def delineate(self, lon, lat, snap_m, thr):
        r, c = self.lonlat_to_rc(lon, lat)
        if not (0 <= r < self.nr and 0 <= c < self.nc):
            raise ValueError("That point is outside the DEM area.")
        if self.fdir[r, c] < 0 and self.acc[r, c] <= 0:
            raise ValueError("No elevation data at that point.")
        clicked = (r, c)
        if snap_m > 0:
            rad = max(0, int(round(snap_m / self.dy_m)))
            r, c = hydro.snap_to_channel(self.acc, r, c, rad, np.float32(max(thr, 1e-9)))
        with self.lock:
            if self.buf is None:
                self.buf = np.empty(self.nr * self.nc, dtype=np.int32)
            m = hydro.upstream_cells(self.fdir, r, c, self.buf)
            r0, r1, c0, c1, mask, emin, emax, esum, ecnt = hydro.basin_window(self.buf, m, self.nc, self.dem)

        bid = self.next_id
        self.next_id += 1
        self.basins[bid] = dict(r0=int(r0), c0=int(c0), mask=mask, outlet=(int(r), int(c)))
        if len(self.basins) > 20:
            self.basins.pop(min(self.basins))

        warnings = []
        if r0 == 0 or c0 == 0 or r1 == self.nr - 1 or c1 == self.nc - 1:
            warnings.append("The basin reaches the edge of the DEM, so it is probably truncated. "
                            "Rebuild with a larger --bbox to include the full upstream area.")

        has_e = ecnt > 0
        area = float(self.acc[r, c])
        olon, olat = self.rc_to_lonlat(r, c)
        clon, clat = self.rc_to_lonlat(*clicked)
        moved = _haversine(float(clon), float(clat), float(olon), float(olat))

        poly = self._polygon(mask, r0, c0)
        streams, sstats = self.streams(bid, thr)
        stats = dict(
            area_km2=area,
            cells=int(m),
            outlet_elev_m=float(self.dem[r, c]),
            min_elev_m=float(emin) if has_e else None,
            mean_elev_m=float(esum / ecnt) if has_e else None,
            max_elev_m=float(emax) if has_e else None,
            relief_m=float(emax - self.dem[r, c]) if has_e else None,
            perimeter_km=poly["properties"]["perimeter_km"],
            snap_distance_m=moved,
            **sstats,
        )
        poly["properties"].update({k: v for k, v in stats.items() if v is not None})
        return dict(id=bid, outlet=[float(olon), float(olat)], clicked=[lon, lat],
                    basin=poly, streams=streams, stats=stats, warnings=warnings)

    def streams(self, bid, thr):
        b = self.basins.get(bid)
        if b is None:
            raise ValueError("Basin expired - click again.")
        return self._trace(b["r0"], b["c0"], b["mask"], thr, with_order=True)

    def view_streams(self, w, s, e, n, thr, max_cells=40_000_000):
        corners_lon = [w, e, w, e]
        corners_lat = [s, s, n, n]
        rc = [self.lonlat_to_rc(lo, la) for lo, la in zip(corners_lon, corners_lat)]
        rs = [p[0] for p in rc]
        cs = [p[1] for p in rc]
        r0, r1 = max(0, min(rs)), min(self.nr - 1, max(rs))
        c0, c1 = max(0, min(cs)), min(self.nc - 1, max(cs))
        if r1 < r0 or c1 < c0:
            return dict(type="FeatureCollection", features=[]), {}
        if (r1 - r0 + 1) * (c1 - c0 + 1) > max_cells:
            raise ValueError("Zoom in to show the network for the whole view.")
        mask = np.ones((r1 - r0 + 1, c1 - c0 + 1), dtype=np.uint8)
        return self._trace(r0, c0, mask, thr, with_order=False)

    def _trace(self, r0, c0, mask, thr, with_order):
        thr = max(float(thr), 1e-9)
        cells, offs = hydro.trace_streams(self.fdir, self.acc, int(r0), int(c0), mask, np.float32(thr))
        n_pts = len(cells)
        if n_pts > 3_000_000:
            raise ValueError(f"Too many channel cells ({n_pts:,}) - raise the threshold.")
        order = hydro.strahler(cells, offs) if with_order else None
        rows = cells // self.nc
        cols = cells - rows * self.nc
        lon, lat = self.rc_to_lonlat(rows, cols)
        lon = np.round(lon, 6)
        lat = np.round(lat, 6)
        accv = self.acc[rows, cols]
        # segment lengths, vectorised: step lengths summed per segment
        step = np.zeros(len(cells))
        if len(cells) > 1:
            step[1:] = _haversine(lon[:-1], lat[:-1], lon[1:], lat[1:]) / 1000.0
            step[offs[:-1]] = 0.0  # first vertex of each segment has no incoming step
        starts = offs[:-1]
        seglen = np.add.reduceat(step, starts) if len(starts) else np.zeros(0)
        # drop vertices in the middle of straight D8 runs (same direction in and out)
        dcell = self.fdir[rows, cols]
        keep = np.ones(len(cells), dtype=bool)
        if len(cells) > 2:
            keep[1:] = dcell[:-1] != dcell[1:]
            keep[offs[:-1]] = True
            keep[offs[1:] - 1] = True
        koffs = np.concatenate([[0], np.cumsum(keep)])[offs]
        xy = np.column_stack([lon[keep], lat[keep]]).tolist()
        feats = []
        total_len = 0.0
        for s in range(len(offs) - 1):
            a, b = offs[s], offs[s + 1]
            if b - a < 2:
                continue
            L = float(seglen[s])
            total_len += L
            props = dict(acc_km2=round(float(accv[b - 2]), 3), length_km=round(L, 3))
            if order is not None:
                props["order"] = int(order[s])
            feats.append(dict(type="Feature", properties=props,
                              geometry=dict(type="LineString", coordinates=xy[koffs[s]:koffs[s + 1]])))
        stats = dict(stream_length_km=round(total_len, 3), stream_segments=len(feats),
                     threshold_km2=thr)
        if order is not None and len(order):
            stats["max_strahler_order"] = int(order.max())
        return dict(type="FeatureCollection", features=feats), stats

    def _polygon(self, mask, r0, c0):
        from rasterio import features
        from affine import Affine
        tr = Affine(self.a, self.b, self.c0 + self.a * c0, self.d, self.e, self.f0 + self.e * r0)
        best = None
        best_area = -1
        for geom, val in features.shapes(mask, mask=mask.astype(bool), transform=tr, connectivity=8):
            ring = np.array(geom["coordinates"][0])
            x, y = ring[:, 0], ring[:, 1]
            ar = 0.5 * abs(np.dot(x[:-1], y[1:]) - np.dot(x[1:], y[:-1]))
            if ar > best_area:
                best_area, best = ar, ring
        lon, lat = self.xy_to_lonlat(best[:, 0], best[:, 1])
        per = float(_haversine(lon[:-1], lat[:-1], lon[1:], lat[1:]).sum()) / 1000.0
        coords = np.column_stack([np.round(lon, 6), np.round(lat, 6)]).tolist()
        return dict(type="Feature", properties=dict(perimeter_km=round(per, 3)),
                    geometry=dict(type="Polygon", coordinates=[coords]))


def _haversine(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = (np.asarray(v, dtype=np.float64) * DEG for v in (lon1, lat1, lon2, lat2))
    h = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * R_EARTH * np.arcsin(np.sqrt(h))


# --------------------------------------------------------------------------------------
# Build cache
# --------------------------------------------------------------------------------------
def build(src_files, bbox, cache_dir):
    import rasterio
    from rasterio.merge import merge

    t0 = time.time()
    log("Reading DEM ...")
    if len(src_files) == 1 and bbox is None:
        with rasterio.open(src_files[0]) as ds:
            dem = ds.read(1, masked=True)
            transform, crs = ds.transform, ds.crs
        valid = ~np.ma.getmaskarray(dem)
        dem = np.ma.filled(dem.astype(np.float32), np.nan)
    else:
        srcs = [rasterio.open(p) for p in src_files]
        arr, transform = merge(srcs, bounds=bbox, nodata=-32767.0, dtype="float32")
        crs = srcs[0].crs
        for s_ in srcs:
            s_.close()
        dem = arr[0]
        valid = dem > -32000
        dem[~valid] = np.nan
    valid &= np.isfinite(dem)
    nr, nc = dem.shape
    geographic = bool(crs.is_geographic)
    a, e = transform.a, transform.e
    if geographic:
        lat = transform.f + e * (np.arange(nr) + 0.5)
        dy = abs(e) * DEG * R_EARTH
        dx_row = a * DEG * R_EARTH * np.cos(lat * DEG)
    else:
        dy = abs(e)
        dx_row = np.full(nr, abs(a))
    cell_km2_row = dx_row * dy / 1e6
    log(f"  {nr:,} x {nc:,} cells ({nr * nc / 1e6:.1f} M), cell ~{math.sqrt(cell_km2_row.mean()) * 1000:.0f} m")
    log("Conditioning DEM (one-time; first run also compiles the engine) ...")
    fdir, acc = hydro.condition(dem, valid, dx_row, dy, cell_km2_row, log=log)
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.save(cache_dir / "fdir.npy", fdir)
    np.save(cache_dir / "acc.npy", acc)
    np.save(cache_dir / "dem.npy", dem)
    meta = dict(transform=list(transform)[:6], crs=crs.to_wkt(), geographic=geographic,
                dy_m=float(dy), cell_km2_row=cell_km2_row.tolist() if geographic else [float(cell_km2_row[0])],
                shape=[nr, nc], built=time.strftime("%Y-%m-%d %H:%M"),
                source=[Path(p).name for p in src_files][:50])
    (cache_dir / "meta.json").write_text(json.dumps(meta))
    log(f"Done in {time.time() - t0:.0f} s. Cached in {cache_dir}")


# --------------------------------------------------------------------------------------
# Web server
# --------------------------------------------------------------------------------------
def make_handler(terrain, info):
    index_html = (HERE / "index.html").read_bytes()

    class H(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _send(self, code, body, ctype="application/json"):
            if not isinstance(body, (bytes, bytearray)):
                body = json.dumps(body, separators=(",", ":")).encode()
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            u = urllib.parse.urlparse(self.path)
            q = {k: v[0] for k, v in urllib.parse.parse_qs(u.query).items()}
            try:
                if u.path in ("/", "/index.html"):
                    return self._send(200, index_html, "text/html; charset=utf-8")
                if u.path.startswith("/static/"):
                    f = (HERE / u.path.lstrip("/")).resolve()
                    if f.is_file() and (HERE / "static") in f.parents:
                        ctype = {".js": "text/javascript", ".css": "text/css", ".png": "image/png"}.get(f.suffix,
                                                                                                       "application/octet-stream")
                        return self._send(200, f.read_bytes(), ctype)
                    return self._send(404, {"error": "not found"})
                if u.path == "/api/info":
                    return self._send(200, info)
                if u.path == "/api/delineate":
                    t = time.time()
                    res = terrain.delineate(float(q["lon"]), float(q["lat"]),
                                            float(q.get("snap", 0)), float(q.get("thr", 1)))
                    res["ms"] = round((time.time() - t) * 1000)
                    return self._send(200, res)
                if u.path == "/api/streams":
                    t = time.time()
                    fc, st = terrain.streams(int(q["id"]), float(q["thr"]))
                    return self._send(200, dict(streams=fc, stats=st, ms=round((time.time() - t) * 1000)))
                if u.path == "/api/view_streams":
                    fc, st = terrain.view_streams(float(q["w"]), float(q["s"]), float(q["e"]),
                                                  float(q["n"]), float(q["thr"]))
                    return self._send(200, dict(streams=fc, stats=st))
                return self._send(404, {"error": "not found"})
            except ValueError as ex:
                return self._send(400, {"error": str(ex)})
            except Exception as ex:  # keep the app alive on any error
                import traceback
                traceback.print_exc()
                return self._send(500, {"error": f"{type(ex).__name__}: {ex}"})

    return H


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dem", help="Your own DEM GeoTIFF (any CRS) instead of Copernicus")
    ap.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"),
                    help="Area in degrees (default: Lesotho + buffer)")
    ap.add_argument("--res", type=int, choices=[30, 90], default=30, help="Copernicus resolution in m (default 30)")
    ap.add_argument("--data", default=str(HERE / "data"), help="Folder for downloaded tiles and cache")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--rebuild", action="store_true", help="Recondition the DEM even if cached")
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args()

    data = Path(args.data)
    if args.dem:
        src = Path(args.dem)
        bbox = tuple(args.bbox) if args.bbox else None
        key = src.stem + ("" if bbox is None else "_" + "_".join(f"{v:g}" for v in bbox))
    else:
        bbox = tuple(args.bbox) if args.bbox else LESOTHO_BBOX
        key = f"cop{args.res}_" + "_".join(f"{v:g}" for v in bbox)
    cache = data / "cache" / key.replace("-", "m")

    if args.rebuild or not (cache / "meta.json").exists():
        if args.dem:
            files = [src]
        else:
            log(f"Getting Copernicus GLO-{args.res} DEM for {bbox} ...")
            files = download_copernicus(bbox, args.res, data / f"tiles_{args.res}m")
        build(files, bbox, cache)

    log("Loading terrain ...")
    t = time.time()
    terrain = Terrain(cache)
    # warm up the compiled functions so the first click is instant
    hydro.snap_to_channel(terrain.acc, 1, 1, 1, np.float32(1.0))
    hydro.trace_streams(terrain.fdir, terrain.acc, 0, 0, np.ones((2, 2), np.uint8), np.float32(1.0))
    hydro.upstream_cells(terrain.fdir, 0, 0, np.empty(terrain.nr * terrain.nc, np.int32))
    hydro.basin_window(np.zeros(1, np.int32), 1, terrain.nc, terrain.dem)
    info = dict(bounds=terrain.bounds_ll(), shape=[terrain.nr, terrain.nc],
                cell_km2=float(np.median(terrain.cell_km2_row)), cell_m=float(terrain.dy_m),
                source=("Your DEM: " + Path(args.dem).name) if args.dem else f"Copernicus GLO-{args.res}",
                built=terrain.meta.get("built"))
    log(f"  ready ({time.time() - t:.1f} s)")

    srv = ThreadingHTTPServer(("127.0.0.1", args.port), make_handler(terrain, info))
    url = f"http://127.0.0.1:{args.port}/"
    log(f"\nBasin Delineator running at {url}\nClose this window (or press Ctrl+C) to quit.\n")
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
