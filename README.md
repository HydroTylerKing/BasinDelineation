# Basin Delineator

Click a point on the map and get its watershed boundary and river network in well under a second. A slider sets the minimum contributing area that defines a stream, and the network redraws as you move it.

Built for Lesotho, but works anywhere with `--bbox` or your own DEM.

## Start it

**Windows:** double-click `Start Basin Delineator (Windows).bat`
**macOS / Linux:** double-click `start.command` (or run `./start.command` in a terminal)

If Python 3.10+ is installed, the launcher uses it. If not, the Windows launcher downloads a private copy of Python into the app's `tools` folder (nothing is installed system-wide). On macOS/Linux, install Python 3.10+ from python.org first.

**First run** downloads the Copernicus GLO-30 DEM for Lesotho plus a buffer (16 tiles, roughly 400–600 MB, free, no login) and conditions it once. On a typical laptop that takes 1–3 minutes after the download and needs about 2 GB of free RAM. Every later start takes a couple of seconds, because the result is cached in `data/cache`.

The app then opens in your browser at http://127.0.0.1:8765. It runs only on your computer; closing the terminal window stops it.

## Using it

- **Click the map** to place an outlet. The basin, its streams and its statistics appear at once.
- **Drag the pin** to move the outlet and re-delineate.
- **Minimum contributing area** (slider or type a value in km²) sets where channels begin. The network updates live.
- **Snap to nearest stream within** moves your click to the nearest cell of the drawn network, so an outlet placed a few pixels off a river still lands on it. Set it to "off" to use the exact cell you clicked (e.g. for a surveyed gauge with the matching threshold).
- **Enter coordinates** for a gauge or dam site.
- **Also show network across the whole view** draws the channel network everywhere on screen, which helps you pick outlets.
- **Keep previous basins** leaves earlier basins outlined for comparison.
- **Download GeoJSON** saves the basin polygon, outlet point (with all statistics) and streams (with Strahler order, upstream area and length). Opens directly in QGIS or ArcGIS.

Statistics: area, perimeter, outlet / mean / max elevation, relief, total stream length, drainage density, number of segments, and maximum Strahler order.

## Options

```
python basin_app.py --res 90                 # 90 m DEM: ~9x less memory, for older machines
python basin_app.py --dem my_lidar.tif       # your own DEM, any projection
python basin_app.py --bbox 27 -30.2 28.6 -29 # different area (W S E N in degrees)
python basin_app.py --port 9000              # different port
python basin_app.py --rebuild                # recondition even if cached
```

On Windows add these after the .bat name in a terminal, or run them with `.venv\Scripts\python basin_app.py ...`.

## How it works

Conditioning (once):

1. **Depression filling** with Priority-Flood (Barnes et al., 2014). A FIFO pit queue makes filled pits and flats drain toward their real spill point.
2. **D8 flow directions** by true steepest descent on the filled surface. East-west cell spacing is corrected for latitude (a Copernicus cell is about 31 m N-S but only 27 m E-W in Lesotho). Flat cells follow the flood path to their outlet, so every cell drains off the DEM with no loops.
3. **Flow accumulation in km²**, using the true area of each cell rather than a cell count.

Every click then walks up the flow-direction tree from the outlet, and traces channel cells (accumulation ≥ threshold) into polylines between channel heads and junctions.

Checks: the sum of accumulation at all DEM outlets equals the DEM's total area, every downstream cell has more accumulation than its donor (so no loops), and basin polygon areas match the accumulated area.

## Accuracy notes

- **Copernicus is a surface model.** In the Lesotho highlands this matters little, but in towns and plantations buildings and trees can raise the surface. In the lowlands along the Mohokare, check small basins against imagery.
- **Road embankments and dam walls** block flow in any DEM, and filling floods the area behind them. If a basin looks wrong just upstream of a road crossing, place the outlet downstream of it or use a DEM with the culvert burned in.
- **Reservoirs and very flat areas** (Katse, Mohale, floodplains) get parallel, straight flow paths across the flat surface. Basin boundaries are not affected, but stream lines drawn across a lake are schematic.
- A yellow warning appears if a basin reaches the edge of the DEM, which means it is truncated. The default area already includes the South African parts of the Caledon and Tele catchments that drain into or along Lesotho. Rebuild with a larger `--bbox` for outlets further downstream (e.g. the Orange below Lesotho).

## Files

- `basin_app.py` - launcher, DEM download, local web server
- `hydro.py` - the hydrology engine (compiled with Numba)
- `index.html`, `static/` - the map interface (Leaflet, bundled)
- `data/` - downloaded tiles and the conditioned cache (created on first run; delete it to start fresh)

Basemaps (Esri, OpenTopoMap, OpenStreetMap) need an internet connection. Delineation itself works offline once the DEM is cached.

DEM: Copernicus GLO-30 / GLO-90, © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018, provided under COPERNICUS by the European Union and ESA.
