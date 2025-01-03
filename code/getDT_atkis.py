""" Calculation of temperature gradients (DT) for a given pixel (u, v).

- variable buffer zone (circular mask) is used to define the size of the data subset
- edge pixels are currently include in the calculation
- rural areas are defined by atkis data and not night lights (could be changed)
- might need to manually change the raster band number used for the calculation

# Todo:
# threshold elevation currently not used ( manual setting in code)
# Check CRS of all files beeing used! They need to be the same.
"""

import os
import math
import numpy as np
from osgeo import gdal
import scipy.io as spio
import multiprocessing as mp
from multiprocessing import Pool, Array
import yaml
import matplotlib.pyplot as plt
import time
import tifffile
import ctypes as c
from tqdm import tqdm
import rasterio
from datetime import datetime

# from code import visual, geodata
import geodata, visual

gdal.DontUseExceptions()
pool = Pool()

import warnings

warnings.filterwarnings("ignore")


def create_mask_circular(h, w, center=None, radius=None):
    """Create a circular mask for a numpy array of shape (h, w).

    based on:  https://stackoverflow.com/questions/44865023/how-can-i-create-a-circular-mask-for-a-numpy-array


    Args:
        h (int): Height of the array.
        w (int): Width of the array.
        center (tuple, optional): Center coordinates of the circle. Default is None (center of the image).
        radius (int, optional): Radius of the circle. Default is None (minimum distance between center and image walls).

    Returns:
        numpy.ndarray: Circular mask.
    """

    if center is None:  # use the middle of the image
        center = (int(w / 2), int(h / 2))

    if radius is None:  # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w - center[0], h - center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)

    mask = dist_from_center <= radius

    return mask


def getDT(u, v):
    """Temperature Gradient (DT) calculation for a given pixel (u, v).
    Calculate DT for a given pixel (u, v) using data subsets and masks.

    Args:
        u (int): Row index of the pixel.
        v (int): Column index of the pixel.
        LST (numpy.ndarray): Land surface temperature data.
        aspect (numpy.ndarray): Aspect data.
        elevation (numpy.ndarray): Elevation data.
        urban_areas (numpy.ndarray): Urban areas data.
        threshold_urban_areas (int): Urban areas threshold.
        threshold_pixel (int): Pixel count threshold.
        circle (numpy.ndarray): Circular mask.

    Returns:
        float: DT value for the pixel (u, v).
    """

    # create DT array from shared memory object
    DT = np.frombuffer(DT_mp.get_obj()).reshape((h, w))

    # check if the pixel+buffer zone is within the image
    if u - buffer < 0:
        u1 = 0
    else:
        u1 = u - buffer

    if u + buffer + 1 > h:
        u2 = h
    else:
        u2 = u + buffer + 1

    if v - buffer < 0:
        v1 = 0
    else:
        v1 = v - buffer

    if v + buffer + 1 > w:
        v2 = w
    else:
        v2 = v + buffer + 1

    # create data subset
    l = LST[u1:u2, v1:v2]
    nl = urban_areas[u1:u2, v1:v2]
    e = elevation[u1:u2, v1:v2]
    asp = aspect[u1:u2, v1:v2]
    r = roi[u1:u2, v1:v2]

    # create circular mask for edge pixels
    if v - buffer < 0 or u - buffer < 0 or v + buffer + 1 > w or u + buffer + 1 > h:
        circle_partial = circle.copy()[0 : u2 - u1, 0 : v2 - v1]

    # create masks for the data selection
    masked_nl = np.ma.masked_where(nl == 1, nl)  #
    mask_nl = masked_nl.mask  # rural areas

    # chose points of similar elevation
    # create mask from elev where values are between -100 and 100
    masked_elev = np.ma.masked_where(
        (e < elevation[u, v] - 100) | (e > elevation[u, v] + 100), e
    )
    mask_elev = ~masked_elev.mask

    # select region of interest
    masked_roi = np.ma.masked_where(r != 1, r)  # ==1
    mask_roi = masked_roi.mask

    # chose points of similar aspect
    # create mask from Asp where values are between -90 and 90
    if aspect[u, v] < 91:
        asp[asp > 269] = asp[asp > 269] - 360
    elif aspect[u, v] > 269:
        asp[asp < 91] = asp[asp < 91] + 360

    masked_asp = np.ma.masked_where(
        (asp < aspect[u, v] - 90) | (asp > aspect[u, v] + 90), asp
    )
    mask_asp = ~masked_asp.mask

    if v - buffer < 0 or u - buffer < 0 or v + buffer + 1 > w or u + buffer + 1 > h:
        # for edge pixels
        try:
            mask_rural = mask_elev & mask_asp & mask_nl & circle_partial  # & mask_roi
        except:
            print(
                mask_elev.shape,
                mask_asp.shape,
                mask_nl.shape,
                circle_partial.shape,
                # mask_roi.shape,
            )
    else:
        # for all other pixels
        mask_rural = mask_elev & mask_asp & mask_nl & circle  # & mask_roi

    L = np.ma.array(l, mask=~mask_rural)

    # only calculate DT if there are enough pixels left
    if L.count() > threshold_pixels:
        # calculate DT from local land surface temperature and rural surroundings
        DT[u, v] = LST[u, v] - np.nanmedian(L.compressed())
    else:
        DT[u, v] = np.nan

    # only used for testing
    if __name__ != "__main__":
        if u == 500 and v == 1427:
            print(
                "row %s, col %s done with %s and %s\n"
                % (u, v, LST[u, v], np.nanmedian(L.compressed()))
            )
            a = LST[u, v]
            b = np.nanmedian(L.compressed())
            return [a, b]
        if u == 500 and v == 1426:
            print(
                "row %s, col %s done with %s and %s\n"
                % (u, v, LST[u, v], np.nanmedian(L.compressed()))
            )
            a = LST[u, v]
            b = np.nanmedian(L.compressed())
            return [a, b]
        if u == 500 and v == 1425:
            print(
                "row %s, col %s done with %s and %s\n"
                % (u, v, LST[u, v], np.nanmedian(L.compressed()))
            )
            a = LST[u, v]
            b = np.nanmedian(L.compressed())  # np.ma.median(l)
            return [a, b]
        if u == 500 and v == 1424:
            print(
                "row %s, col %s done with %s and %s\n"
                % (u, v, LST[u, v], np.nanmedian(L.compressed()))
            )
            a = LST[u, v]
            b = np.nanmedian(L.compressed())  # np.ma.median(l)
            return [a, b]
        if u == 500 and v == 1423:
            print(
                "row %s, col %s done with %s and %s\n"
                % (u, v, LST[u, v], np.nanmedian(L.compressed()))
            )
            a = LST[u, v]
            b = np.nanmedian(L.compressed())  # np.ma.median(l)
            return [a, b]


def load_mock():
    """Load mock data for testing."""

    # data import
    file_in = "/data/Frankfurt/FrankfurtLST_summermean_2015-20.tif"
    ds = gdal.Open(file_in)
    LST = np.array(ds.GetRasterBand(1).ReadAsArray())

    h, w = LST.shape[:2]

    # Test data
    # same as in matlab

    aspect = np.zeros([h, w]) - 1

    mat = spio.loadmat("/data/elev.mat", squeeze_me=True)
    elevation = mat["elev"]

    # elevation= np.zeros([h,w])+50

    mat = spio.loadmat("/data/NL.mat", squeeze_me=True)
    urban_areas = mat["NL"]

    # Create the circular mask
    buffer = 500  # --> 50 km radius da 500 * 100 = 50.000 m
    circle = create_mask_circular(h=(buffer * 2) + 1, w=(buffer * 2) + 1, radius=buffer)

    threshold_urban_areas = 15
    threshold_pixels = 30

    DT_mp = Array(c.c_double, h * w)  # for shared memory

    return (
        DT_mp,
        h,
        w,
        LST,
        elevation,
        urban_areas,
        aspect,
        circle,
        threshold_urban_areas,
        threshold_pixels,
        buffer,
    )


def load_config(file_path):
    """Load configuration file."""
    print(os.getcwd())

    with open(file_path, "r") as config_file:
        config = yaml.safe_load(config_file)

    path_LST = config["path_LST"]
    ds = gdal.Open(path_LST)
    print(ds)

    # normaly raster band 1 (mean) MODIS 3(90th)
    LST = np.array(ds.GetRasterBand(1).ReadAsArray())
    LST[LST == 0] = np.nan
    LST[LST == 9999] = np.nan
    LST[LST == -99999] = np.nan

    threshold_urban_areas = config["threshold_urban_areas"]
    threshold_pixels = config["threshold_pixels"]
    threshold_elevation = config["threshold_elevation"]

    roi = tifffile.imread(config["path_roi"])
    roi[roi != 1] = np.nan

    buffer = config["buffer"]
    circle = create_mask_circular(h=(buffer * 2) + 1, w=(buffer * 2) + 1, radius=buffer)

    elevation = tifffile.imread(config["path_elevation"])
    elevation[elevation == 0] = np.nan

    aspect = tifffile.imread(config["path_aspect"])
    aspect[aspect == -9999] = np.nan

    urban_areas = tifffile.imread(config["path_urban_areas"]).astype(np.float32)
    urban_areas[urban_areas != 1] = np.nan  # for the definition of urban areas

    h, w = LST.shape[:2]
    DT_mp = Array(c.c_double, h * w)  # for shared memory

    return (
        path_LST,
        ds,
        DT_mp,
        LST,
        elevation,
        urban_areas,
        aspect,
        circle,
        threshold_urban_areas,
        threshold_pixels,
        buffer,
        roi,
    )


if __name__ == "__main__":
    print("Processing started...")

    config_file_path = "./params.yaml"
    (
        path_LST,
        ds,
        DT_mp,
        LST,
        elevation,
        urban_areas,
        aspect,
        circle,
        threshold_urban_areas,
        threshold_pixels,
        buffer,
        roi,
    ) = load_config(config_file_path)

    h, w = LST.shape[:2]

    # without edge pixels
    # U = range(buffer, h - buffer)
    # V = range(buffer, w - buffer)

    # with edge pixels
    U = range(0, h)
    V = range(0, w)

    start_time = time.perf_counter()

    for u in tqdm(U):
        with Pool() as p:
            # calculation of temperature gradient for each pixel (multiprocessing)
            p.starmap(getDT, [(u, v) for v in V])

            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            if elapsed_time == 3600:
                start_time = time.perf_counter()
                src = rasterio.open(path_LST)

                filename = path_LST.split("/")[-1].split(".")[0]
                path_in = "/".join(path_LST.split("/")[:-1]) + "/"
                try:
                    # get the current time
                    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

                    with rasterio.open(
                        path_in
                        + "DT_"
                        + filename
                        + "_radius_"
                        + str(buffer)
                        + "_"
                        + now
                        + ".tif",
                        "w",
                        height=h,
                        width=w,
                        count=1,
                        dtype=src.dtypes[0],
                        crs=src.crs,
                        transform=src.transform,
                    ) as dst:
                        dst.write(DT, 1)
                except Exception as e:
                    print(e)
                    print("Error while writing DT to file.")

                    # write tiff file not using rasterio
                    tifffile.imwrite(
                        path_in
                        + "DT_"
                        + filename
                        + "_radius_"
                        + str(buffer)
                        + "_tif"
                        + now
                        + ".tif",
                        DT,
                    )

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time

    print("Processing complete.")
    print(f"Elapsed time: {elapsed_time:.4f} seconds")

    # recreate DT array from shared memory object
    DT = np.frombuffer(DT_mp.get_obj()).reshape((h, w))

    # data export and visulatisation

    visual.temperature_gradients(DT)

    # bb for Frankfurt
    # north = 50.67755843831869
    # south = 49.56454580129461
    # west = 7.771325522917981
    # east = 9.503277390700418

    # visual.temperature_gradients_on_map(
    #   DT, lon_range=[west, east], lat_range=[south, north]
    # )

    # export DT
    # proj = ds.GetProjection()
    # srs = osr.SpatialReference(wkt=proj)

    # if srs.IsProjected:
    #     CS = srs.GetAttrValue("projcs")
    # CS = srs.GetAttrValue("geogcs")
    # geodata.export_GTiff(
    #     DT, "./data/tmp/DT_200_buffer_100km.tif", ds.GetGeoTransform(), CS
    # )
    src = rasterio.open(path_LST)

    filename = path_LST.split("/")[-1].split(".")[0]

    path_in = "/".join(path_LST.split("/")[:-1]) + "/"
    try:
        # get the current time
        now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        with rasterio.open(
            path_in + "DT_" + filename + "_radius_" + str(buffer) + "_" + now + ".tif",
            "w",
            height=h,
            width=w,
            count=1,
            dtype=src.dtypes[0],
            crs=src.crs,
            transform=src.transform,
        ) as dst:
            dst.write(DT, 1)
    except Exception as e:
        print(e)
        print("Error while writing DT to file.")

        # write tiff file not using rasterio
        tifffile.imwrite(
            path_in
            + "DT_"
            + filename
            + "_radius_"
            + str(buffer)
            + "_tif"
            + now
            + ".tif",
            DT,
        )
else:
    (
        DT_mp,
        h,
        w,
        LST,
        elevation,
        urban_areas,
        aspect,
        circle,
        threshold_urban_areas,
        threshold_pixels,
        buffer,
    ) = load_mock()
