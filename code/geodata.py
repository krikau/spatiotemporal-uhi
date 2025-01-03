from osgeo import gdal, osr
import numpy as np


def get_bbox(path_raster):
    """Get bounding box of raster."""
    # code based on: https://gis.stackexchange.com/questions/57710/determining-coordinates-of-corners-of-raster-layer-using-pyqgis/57711#57711

    bag = gdal.Open(path_raster)
    bag_gtrn = bag.GetGeoTransform()
    bag_proj = bag.GetProjectionRef()
    bag_srs = osr.SpatialReference(bag_proj)
    geo_srs = bag_srs.CloneGeogCS()  # New srs obj to go from x,y -> φ,λ
    transform = osr.CoordinateTransformation(bag_srs, geo_srs)

    bag_bbox_cells = (
        (0.0, 0.0),
        (0, bag.RasterYSize),
        (bag.RasterXSize, bag.RasterYSize),
        (bag.RasterXSize, 0),
    )

    geo_pts = []
    for x, y in bag_bbox_cells:
        x2 = bag_gtrn[0] + bag_gtrn[1] * x + bag_gtrn[2] * y
        y2 = bag_gtrn[3] + bag_gtrn[4] * x + bag_gtrn[5] * y
        geo_pt = transform.TransformPoint(x2, y2)[:2]
        geo_pts.append(geo_pt)

        # Print each step of transformation
        print(
            f"Pixel Coord: ({x}, {y}) -> Proj Coords: ({x2}, {y2}) -> (φ,λ) coords: {geo_pt}"
        )

    print(geo_pts)

    # Get bounding box
    north = geo_pts[0][1]
    south = geo_pts[1][1]
    east = geo_pts[3][0]
    west = geo_pts[0][0]

    lat_range = (south, north)
    lon_range = (west, east)

    return lat_range, lon_range


def export_GTiff(img, path_output, geoTransform, geoCS="WGS84"):
    """Export array to GeoTiff."""
    # the image information does not get reprojected only the metadata changes

    # parameter settings
    driver = gdal.GetDriverByName("GTiff")
    try:
        nYSize, nXSize, b = img.shape
        nBands = 1
    except ValueError:
        nYSize, nXSize = img.shape
        nBands = 1
    eType = gdal.GDT_Float64

    # Returns:NULL on failure, or a new GDALDataset.
    dataset_output = driver.Create(path_output, nXSize, nYSize, nBands, eType)

    # set gepgraphic transformation
    dataset_output.SetGeoTransform(geoTransform)

    # set spatial reference
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS(geoCS)
    dataset_output.SetProjection(srs.ExportToWkt())

    # write array to raster
    export = dataset_output.GetRasterBand(1).WriteArray(img)

    # free memory of driver
    dataset_output = None

    return export


def combine_tiles(tiles, h, w):
    """Combine tiles to recreate the full image."""
    z, x, y = np.array(tiles).shape

    # Create an empty image of the same shape as the original image
    full_image = np.zeros((h, w))

    idx = 0
    for i in range(0, h, x):
        for j in range(0, w, y):
            tile = tiles[idx]
            full_image[i : i + x, j : j + y] = tile
            idx += 1

    return full_image


def get_devisor(x):
    """Get devisor of x."""
    a = []

    i = 1
    while i <= x:
        if x % i == 0:
            d = x // i
            a.append(d)
            i = i + 1

        else:
            i = i + 1

    a.reverse()

    return a


def get_tile_size(img):
    """Get tile size."""
    h, w = img.shape[:2]

    a = get_devisor(h)
    b = get_devisor(w)

    x = a[int(len(a) / 2 + 1)]
    y = b[int(len(b) / 2 + 1)]

    return x, y


def create_tiles(img):
    """Create tiles."""
    h, w = img.shape[:2]
    x, y = get_tile_size(img)

    tiles = []
    for i in range(0, h, x):
        for j in range(0, w, y):
            tiles.append(img[i : i + x, j : j + y])

    return tiles
