# import the necessary packages

# https://geopandas.org/en/stable/docs/user_guide/set_operations.html

import geopandas as gpd
import os

import rasterio
import rasterio.plot
import rasterio.mask

# tutorial: https://automating-gis-processes.github.io/CSC/notebooks/L5/zonal-statistics.html

from rasterstats import zonal_stats

import numpy as np
import pandas as pd

import warnings

warnings.filterwarnings("ignore")


def calc_stats_urban_areas(parameter, raster, urban, roi):
    """Calculate zonal statistics for the urban areas (per district) within the roi.

    "Hotpsots"

    Args:
        parameter (str): Name of the parameter.
        raster (rasterio.io.DatasetReader): Raster file.
        urban (geopandas.geodataframe.GeoDataFrame): Urban areas.
        roi (geopandas.geodataframe.GeoDataFrame): Region of interest.

    Returns:
        geopandas.geodataframe.GeoDataFrame: Urban areas with zonal statistics.

    """

    # load the input data
    shape_raster = roi.copy()

    affine = raster.transform

    # mask out the non urban areas
    out_image, transformed = rasterio.mask.mask(
        raster, urban.geometry
    )  # or raster_geometry_mask?
    out = out_image[0]
    out[out == raster.nodata] = np.nan  # default no data values
    out[out == 9999] = np.nan  # because PET has 9999 as no data value
    out[out == -99999] = np.nan  # because Ta has -99999 as no data value
    out[out == 0] = np.nan

    # calculate zonal statistics for the urban areas within the roi
    zs = zonal_stats(
        roi, out, affine=affine, stats=["min", "max", "mean", "median", "count"]
    )

    # Extract mean values from the results
    values_min = [feature["min"] for feature in zs]
    values_max = [feature["max"] for feature in zs]
    values_mean = [feature["mean"] for feature in zs]
    values_median = [feature["median"] for feature in zs]
    values_count = [feature["count"] for feature in zs]

    # Add a new column 'MeanValue' to the GeoDataFrame
    shape_raster["min_" + parameter] = values_min
    shape_raster["max_" + parameter] = values_max
    shape_raster["mean_" + parameter] = values_mean
    shape_raster["median_" + parameter] = values_median
    shape_raster["count_" + parameter] = values_count

    # Export the new GeoDataFrame
    sorted_gdf = shape_raster.sort_values(
        by="mean_" + parameter, ascending=False
    )  # Sort the DataFrame based on a specific column
    sorted_gdf = sorted_gdf.reset_index(drop=True)  # reset the index

    # write_dataframe(
    #    sorted_gdf, path_out + ".gpkg"
    # )  # write the GeoDataFrame to a new geopackage

    return sorted_gdf


def calc_stats_urban_areas_type(parameter, raster, urban, path_out, col="objektart"):
    """Calculate zonal statistics for the urban areas (per district) within the roi based on the col value.


    "Hotpsots"

    Args:
        parameter (str): Name of the parameter.
        raster (rasterio.io.DatasetReader): Raster file.
        urban (geopandas.geodataframe.GeoDataFrame): Urban areas.
        path_out (str): Output path.
        col (str): Column name.

    Returns:
        geopandas.geodataframe.GeoDataFrame: Urban areas with zonal statistics.

    """

    affine = raster.transform

    atkis_data = []
    for atkis in urban[col].unique():
        test = urban[urban[col] == atkis].dissolve()
        # write the shape file to a new geopackage

        # Create a Mask
        out_image, transformed = rasterio.mask.mask(
            raster, test.geometry
        )  # or raster_geometry_mask?
        out = out_image[0]
        out[out == raster.nodata] = np.nan  # default no data values
        out[out == 0] = np.nan  # default no data values
        out[out == 9999] = np.nan  # because PET has 9999 as no data value
        out[out == -99999] = np.nan  # because Ta has -99999 as no data value

        # calculate zonal statistics
        zs = zonal_stats(
            test, out, affine=affine, stats=["min", "max", "mean", "median", "count"]
        )

        # Extract mean values from the results
        values_min = [feature["min"] for feature in zs]
        values_max = [feature["max"] for feature in zs]
        values_mean = [feature["mean"] for feature in zs]
        values_median = [feature["median"] for feature in zs]
        values_count = [feature["count"] for feature in zs]

        # Add a new column 'MeanValue' to the GeoDataFrame
        test["min_" + parameter + "_" + atkis] = values_min
        test["max_" + parameter + "_" + atkis] = values_max
        test["mean_" + parameter + "_" + atkis] = values_mean
        test["median_" + parameter + "_" + atkis] = values_median
        test["count_" + parameter + "_" + atkis] = values_count

        # correct fid
        test["fid"] = test.index
        atkis_data.append(test)

    i = 0
    for data in atkis_data:
        data.to_file(
            path_out,
            layer=urban[col].unique()[i],
            driver="GPKG",
        )
        i = i + 1

    dict_atkis = dict()
    for name in urban["objektart"].unique():
        dict_atkis[name] = gpd.read_file(path_out, layer=name)

    atkis_mean = list()
    for name in urban["objektart"].unique():
        atkis_mean.append(
            [name, dict_atkis[name]["mean_" + parameter + "_" + name].values[0]]
        )
    atkis_mean = pd.DataFrame(atkis_mean, columns=["objektart", "mean_DT"]).sort_values(
        "mean_DT"
    )

    return atkis_mean


def calc_stats_lcz(path_raster, LCZ, path_out):
    r = rasterio.open(path_raster)

    # create new data frame
    lcz_stats = pd.DataFrame(
        columns=[
            "parameter",
            "LCZ",
            "Mean",
            "Median",
            "Min",
            "Max",
            "Std",
            "Var",
            "Quantile_10",
            "Quantile_90",
            "Quantile_95",
            "N",
        ]
    )

    for i, lcz in enumerate(LCZ):
        # print("LCZ: " + str(i))
        try:
            out_image, transformed = rasterio.mask.mask(
                r, lcz.geometry
            )  # or raster_geometry_mask?
            out = out_image[0]
            out[out == r.nodata] = np.nan  # default no data values
            out[out == 9999] = np.nan
            out[out == -99999] = np.nan
            out[out == 0] = np.nan

            # write statistics to file

            row_values = [
                i + 1,
                np.nanmean(out),
                np.nanmedian(out),
                np.nanmin(out),
                np.nanmax(out),
                np.nanstd(out),
                np.nanvar(out),
                np.nanquantile(out, 0.1),
                np.nanquantile(out, 0.9),
                np.nanquantile(out, 0.95),
                np.count_nonzero(~np.isnan(out)),
            ]

            # Add to data frame
            lcz_stats = pd.concat(
                [
                    lcz_stats,
                    pd.DataFrame(
                        {
                            "parameter": [r.name.split("/")[-1].split(".")[0]],
                            "LCZ": [i + 1],
                            "Mean": [np.nanmean(out)],
                            "Median": [np.nanmedian(out)],
                            "Min": [np.nanmin(out)],
                            "Max": [np.nanmax(out)],
                            "Std": [np.nanstd(out)],
                            "Var": [np.nanvar(out)],
                            "Quantile_10": [np.nanquantile(out, 0.1)],
                            "Quantile_90": [np.nanquantile(out, 0.9)],
                            "Quantile_95": [np.nanquantile(out, 0.95)],
                            "N": [np.count_nonzero(~np.isnan(out))],
                        }
                    ),
                ],
                ignore_index=True,
            )
        except:
            print("Error with LCZ: " + str(i + 1))

    if not os.path.exists(path_out + "/statistics_lcz_all.csv"):
        lcz_stats.to_csv(path_out + "/statistics_lcz_all.csv", index=False)
    else:
        lcz_stats.to_csv(
            path_out + "/statistics_lcz_all.csv", mode="a", header=False, index=False
        )

    r.close()
