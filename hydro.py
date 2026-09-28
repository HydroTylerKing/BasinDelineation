"""
Hydrology engine for the basin delineation app.

All heavy lifting is done once, when the DEM is prepared:
  1. Priority-Flood depression filling (Barnes et al. 2014, with a FIFO "pit queue"
     so that filled depressions and flats drain toward their true spill point).
  2. D8 flow directions using true steepest descent on the filled surface, with
     distances corrected for latitude on geographic DEMs. Cells with no lower
     neighbour (flats / filled pits) drain along the flood path to the spill point,
     which guarantees every cell reaches the DEM edge and no loops can form.
  3. Flow accumulation in km² (cell areas corrected for latitude), computed in a
     single pass using the flood's processing order as a topological order.

After that, every click is cheap:
  * watershed = walk up the D8 tree from the outlet (no visited array needed),
  * streams   = cells in the watershed with accumulation >= threshold, traced into
    polylines between channel heads and junctions, with Strahler order.
"""
import heapq
import numpy as np
from numba import njit

# D8 neighbour offsets. Code k points from a cell to neighbour (DR[k], DC[k]).
# 0=E 1=SE 2=S 3=SW 4=W 5=NW 6=N 7=NE ; the opposite direction is (k + 4) % 8.
DR = np.array([0, 1, 1, 1, 0, -1, -1, -1], dtype=np.int64)
DC = np.array([1, 1, 0, -1, -1, -1, 0, 1], dtype=np.int64)


# --------------------------------------------------------------------------------------
# Pre-processing
# --------------------------------------------------------------------------------------
@njit(cache=True)
def _flood(dem, valid):
    """Priority-Flood with pit queue. Returns (filled, parent_dir, order, n_ordered)."""
    nr, nc = dem.shape
    n = nr * nc
    filled = dem.copy()
    closed = np.zeros((nr, nc), dtype=np.uint8)
    parent = np.full((nr, nc), -1, dtype=np.int8)
    order = np.empty(n, dtype=np.int32)
    k = 0

    heap = [(np.float64(0.0), np.int64(0), np.int64(0))]
    heap.pop()
    seq = 0

    # growable ring buffer for the pit queue
    cap = 1 << 20
    pit = np.empty(cap, dtype=np.int64)
    ph = 0  # head
    pn = 0  # number of items

    for r in range(nr):
        for c in range(nc):
            if not valid[r, c]:
                closed[r, c] = 1
    # seed: valid cells on the raster edge or touching nodata
    for r in range(nr):
        for c in range(nc):
            if not valid[r, c]:
                continue
            edge = r == 0 or c == 0 or r == nr - 1 or c == nc - 1
            if not edge:
                for j in range(8):
                    if not valid[r + DR[j], c + DC[j]]:
                        edge = True
                        break
            if edge:
                closed[r, c] = 1
                heapq.heappush(heap, (np.float64(filled[r, c]), np.int64(seq), np.int64(r * nc + c)))
                seq += 1

    while pn > 0 or len(heap) > 0:
        if pn > 0:
            idx = pit[ph]
            ph = (ph + 1) % cap
            pn -= 1
        else:
            item = heapq.heappop(heap)
            idx = item[2]
        order[k] = idx
        k += 1
        r = idx // nc
        c = idx - r * nc
        e = filled[r, c]
        for j in range(8):
            rr = r + DR[j]
            cc = c + DC[j]
            if rr < 0 or cc < 0 or rr >= nr or cc >= nc:
                continue
            if closed[rr, cc]:
                continue
            closed[rr, cc] = 1
            parent[rr, cc] = (j + 4) % 8  # neighbour drains back toward the cell that reached it
            if filled[rr, cc] <= e:
                filled[rr, cc] = e
                if pn == cap:  # grow ring buffer, unrolling it
                    newp = np.empty(cap * 2, dtype=np.int64)
                    for t in range(pn):
                        newp[t] = pit[(ph + t) % cap]
                    pit = newp
                    ph = 0
                    cap *= 2
                pit[(ph + pn) % cap] = rr * nc + cc
                pn += 1
            else:
                heapq.heappush(heap, (np.float64(filled[rr, cc]), np.int64(seq), np.int64(rr * nc + cc)))
                seq += 1
    return filled, parent, order, k


@njit(cache=True)
def _d8(filled, parent, valid, dx_row, dy):
    """Steepest-descent D8 on the filled surface; flats fall back to the flood path."""
    nr, nc = filled.shape
    fdir = np.full((nr, nc), -1, dtype=np.int8)
    for r in range(nr):
        dx = dx_row[r]
        dd = np.sqrt(dx * dx + dy * dy)
        for c in range(nc):
            if not valid[r, c]:
                continue
            e = filled[r, c]
            best = -1
            bs = 0.0
            for j in range(8):
                rr = r + DR[j]
                cc = c + DC[j]
                if rr < 0 or cc < 0 or rr >= nr or cc >= nc:
                    continue
                if not valid[rr, cc]:
                    continue
                drop = e - filled[rr, cc]
                if drop <= 0:
                    continue
                if DR[j] == 0:
                    s = drop / dx
                elif DC[j] == 0:
                    s = drop / dy
                else:
                    s = drop / dd
                if s > bs:
                    bs = s
                    best = j
            if best >= 0:
                fdir[r, c] = best
            else:
                fdir[r, c] = parent[r, c]  # -1 for edge outlets
    return fdir


@njit(cache=True)
def _accumulate(fdir, order, k, cell_km2_row):
    nr, nc = fdir.shape
    acc = np.zeros((nr, nc), dtype=np.float32)
    for t in range(k):
        idx = order[t]
        r = idx // nc
        acc[r, idx - r * nc] = cell_km2_row[r]
    for t in range(k - 1, -1, -1):
        idx = order[t]
        r = idx // nc
        c = idx - r * nc
        j = fdir[r, c]
        if j >= 0:
            acc[r + DR[j], c + DC[j]] += acc[r, c]
    return acc


def condition(dem, valid, dx_row, dy, cell_km2_row, log=print):
    """Run fill -> D8 -> accumulation. Returns (fdir int8, acc float32 km²)."""
    import time
    t = time.time()
    dem = np.ascontiguousarray(dem, dtype=np.float32)
    valid = np.ascontiguousarray(valid, dtype=np.bool_)
    filled, parent, order, k = _flood(dem, valid)
    log(f"  depressions filled   ({time.time() - t:.1f} s)")
    t = time.time()
    fdir = _d8(filled, parent, valid, dx_row.astype(np.float64), float(dy))
    del filled, parent
    log(f"  flow directions      ({time.time() - t:.1f} s)")
    t = time.time()
    acc = _accumulate(fdir, order, k, cell_km2_row.astype(np.float32))
    log(f"  flow accumulation    ({time.time() - t:.1f} s)")
    return fdir, acc


# --------------------------------------------------------------------------------------
# Per-click queries
# --------------------------------------------------------------------------------------
@njit(cache=True)
def upstream_cells(fdir, r0, c0, buf):
    """All cells draining to (r0, c0). Fills buf[0:m] with flat indices, returns m."""
    nr, nc = fdir.shape
    buf[0] = r0 * nc + c0
    head = 0
    tail = 1
    while head < tail:
        idx = buf[head]
        head += 1
        r = idx // nc
        c = idx - r * nc
        for j in range(8):
            rr = r + DR[j]
            cc = c + DC[j]
            if rr < 0 or cc < 0 or rr >= nr or cc >= nc:
                continue
            if fdir[rr, cc] == (j + 4) % 8:
                buf[tail] = rr * nc + cc
                tail += 1
    return tail


@njit(cache=True)
def basin_window(buf, m, nc, dem):
    """Bounding window, mask and elevation stats of the m cells in buf."""
    r0 = 1 << 30
    c0 = 1 << 30
    r1 = -1
    c1 = -1
    for t in range(m):
        idx = buf[t]
        r = idx // nc
        c = idx - r * nc
        r0 = min(r0, r)
        r1 = max(r1, r)
        c0 = min(c0, c)
        c1 = max(c1, c)
    mask = np.zeros((r1 - r0 + 1, c1 - c0 + 1), dtype=np.uint8)
    emin = np.inf
    emax = -np.inf
    esum = 0.0
    ecnt = 0
    for t in range(m):
        idx = buf[t]
        r = idx // nc
        c = idx - r * nc
        mask[r - r0, c - c0] = 1
        e = dem[r, c]
        if np.isfinite(e):
            emin = min(emin, e)
            emax = max(emax, e)
            esum += e
            ecnt += 1
    return r0, r1, c0, c1, mask, emin, emax, esum, ecnt


@njit(cache=True)
def trace_streams(fdir, acc, r0, c0, mask, thr):
    """
    Trace channel cells (acc >= thr) inside a window starting at (r0, c0) whose shape is
    mask.shape, restricted to mask==1. Returns (cells, offsets): segment s is
    cells[offsets[s]:offsets[s+1]] as flat indices of the full raster.
    Segments run from a channel head or junction down to the next junction / domain exit.
    """
    h, w = mask.shape
    nr, nc = fdir.shape
    ndon = np.zeros((h, w), dtype=np.int8)
    isstr = np.zeros((h, w), dtype=np.uint8)
    nstr = 0
    for i in range(h):
        for j in range(w):
            if mask[i, j] and acc[r0 + i, c0 + j] >= thr:
                isstr[i, j] = 1
                nstr += 1
    for i in range(h):
        for j in range(w):
            if not isstr[i, j]:
                continue
            d = fdir[r0 + i, c0 + j]
            if d < 0:
                continue
            ii = i + DR[d]
            jj = j + DC[d]
            if ii >= 0 and jj >= 0 and ii < h and jj < w and isstr[ii, jj]:
                ndon[ii, jj] += 1
    nstart = 0
    for i in range(h):
        for j in range(w):
            if isstr[i, j] and ndon[i, j] != 1:
                nstart += 1
    cells = np.empty(nstr + nstart + 1, dtype=np.int64)
    offsets = np.empty(nstart + 1, dtype=np.int64)
    p = 0
    s = 0
    for i in range(h):
        for j in range(w):
            if not (isstr[i, j] and ndon[i, j] != 1):
                continue
            offsets[s] = p
            s += 1
            ci = i
            cj = j
            cells[p] = (r0 + ci) * nc + (c0 + cj)
            p += 1
            while True:
                d = fdir[r0 + ci, c0 + cj]
                if d < 0:
                    break
                ni = ci + DR[d]
                nj = cj + DC[d]
                if ni < 0 or nj < 0 or ni >= h or nj >= w or not isstr[ni, nj]:
                    break
                cells[p] = (r0 + ni) * nc + (c0 + nj)
                p += 1
                if ndon[ni, nj] != 1:
                    break  # reached a junction; it starts its own segment
                ci = ni
                cj = nj
    offsets[s] = p
    return cells[:p], offsets[: s + 1]


@njit(cache=True)
def snap_to_channel(acc, r, c, radius, thr):
    """Nearest cell within radius whose accumulation >= thr (i.e. on the drawn network).
    Ties in distance go to the larger accumulation. Falls back to the clicked cell."""
    nr, nc = acc.shape
    best_d = radius * radius + 1
    best_a = -1.0
    br = r
    bc = c
    for i in range(max(0, r - radius), min(nr, r + radius + 1)):
        for j in range(max(0, c - radius), min(nc, c + radius + 1)):
            d = (i - r) * (i - r) + (j - c) * (j - c)
            a = acc[i, j]
            if d > radius * radius or a < thr:
                continue
            if d < best_d or (d == best_d and a > best_a):
                best_d = d
                best_a = a
                br = i
                bc = j
    return br, bc


def strahler(cells, offsets):
    """Strahler order per segment (segments come out of trace_streams)."""
    ns = len(offsets) - 1
    if ns == 0:
        return np.zeros(0, dtype=np.int32)
    first = cells[offsets[:-1]]
    last = cells[offsets[1:] - 1]
    start_of = {int(f): s for s, f in enumerate(first)}
    down = np.full(ns, -1, dtype=np.int64)
    for s in range(ns):
        d = start_of.get(int(last[s]), -1)
        if d != s:
            down[s] = d
    # topological sweep from channel heads downstream
    indeg = np.zeros(ns, dtype=np.int64)
    for s in range(ns):
        if down[s] >= 0:
            indeg[down[s]] += 1
    order_ = np.ones(ns, dtype=np.int32)
    maxin = np.zeros(ns, dtype=np.int32)
    cnt = np.zeros(ns, dtype=np.int32)
    stack = [s for s in range(ns) if indeg[s] == 0]
    while stack:
        s = stack.pop()
        if cnt[s] >= 2:
            order_[s] = maxin[s] + 1
        elif maxin[s] > 0:
            order_[s] = maxin[s]
        d = down[s]
        if d >= 0:
            o = order_[s]
            if o > maxin[d]:
                maxin[d] = o
                cnt[d] = 1
            elif o == maxin[d]:
                cnt[d] += 1
            indeg[d] -= 1
            if indeg[d] == 0:
                stack.append(d)
    return order_
