import ee  # Import the Earth Engine library.
import geemap  # Visualisierung
import yaml  # Configuration file

# Trigger the authentication flow.
ee.Authenticate()

# Initialize the library.
ee.Initialize()


def load_config(file_path):
    """Load configuration file."""
    with open(file_path, "r") as config_file:
        config = yaml.safe_load(config_file)

    # Access configuration parameters
    # google drive folder
    folder_drive = config["folder_drive"]

    # Dates
    date_start = config["date_start"]
    date_end = config["date_end"]
    date_median = config["date_median"]

    # buffer zone around the region of interest
    buffer_value = config["buffer"]

    # region of interest
    north = config["north"]
    south = config["south"]
    west = config["west"]
    east = config["east"]
    bb = ee.Geometry.BBox(west, south, east, north)

    # EPSG code
    EPSG_code = str(config["EPSG_code"])

    # administrative district
    distric_name = str(config["distric_name"])
    distric = (
        ee.FeatureCollection("FAO/GAUL/2015/level2")
        .filter(ee.Filter.eq("ADM1_NAME", distric_name))
        .union()
    )

    print("Configuration file loaded.")
    print("Google Drive folder: ", folder_drive)
    print("Date start: ", date_start)
    print("Date end: ", date_end)
    print("Date median: ", date_median)
    print("Buffer value: ", buffer_value)
    print("Region of interest: ", bb)
    print("Administrative district: ", distric_name)

    return (
        folder_drive,
        date_start,
        date_end,
        date_median,
        buffer_value,
        bb,
        distric,
        EPSG_code,
    )


def get_MODIS(
    daytime="Night", date_start="2013-1-1", date_end="2023-1-1", months=[6, 8]
):
    MOD = (
        ee.ImageCollection(
            "MODIS/061/MYD11A1"
        )  # terra: "MODIS/061/MOD11A1" # aqua: MODIS/061/MYD11A1
        .select(
            [
                "LST_" + daytime + "_1km",
                "QC_" + daytime,
                daytime + "_view_time",
                daytime + "_view_angle",
            ]
        )
        .filterDate(date_start, date_end)
        .filter(ee.Filter.calendarRange(months[0], months[1], "month"))
    )

    return MOD


def filter_MODIS(image, daytime="Night"):
    """Filter MODIS data to remove pixels with low quality."""

    def filter(image):
        bitMask7 = 1 << 7  # binary with 1s as the 7th bit zeros elsewhere
        bitMask6 = 1 << 6  # binary with 1s as the 6th bit zeros elsewhere
        qa = image.select("QC_" + daytime)
        mymask = (
            qa.bitwiseAnd(bitMask6).And(qa.bitwiseAnd(bitMask7)).eq(0)
        )  # thus is 0 when sixth and seventh bit are 1
        image = image.updateMask(mymask)
        temp = (
            image.select("LST_" + daytime + "_1km")
            .multiply(0.02)
            .subtract(273.15)
            .rename("LST")
        )  # Transfering to Celsius
        time = (
            image.select("" + daytime + "_view_time")
            .multiply(0.1)
            .rename("observation_time")
        )
        zenith = (
            image.select("" + daytime + "_view_angle")
            .subtract(65)
            .rename("zenith_angle")
        )

        K1error = qa.bitwiseAnd(bitMask6).Or(qa.bitwiseAnd(bitMask7)).eq(0).multiply(1)
        K2error = qa.bitwiseAnd(bitMask6).neq(0).multiply(2)
        K3error = qa.bitwiseAnd(bitMask7).neq(0).multiply(3)
        uncertainty = K1error.add(K2error).add(K3error).rename("uncertainty").double()
        return temp.addBands([time, zenith, uncertainty]).set(
            "system:time_start", image.get("system:time_start")
        )

    MOD_clean = image.map(filter)

    return MOD_clean


def reduce_MODIS(image, roi):
    def getallthedata(imagecol):
        LST_mean = imagecol.select("LST").mean().rename("LST_mean")
        LST_std = imagecol.select("LST").reduce(ee.Reducer.stdDev()).rename("LST_std")
        LST_median = imagecol.select("LST").median().rename("LST_median")
        LST_p90 = (
            imagecol.select("LST").reduce(ee.Reducer.percentile([90])).rename("LST_p90")
        )
        Zenith = imagecol.select("zenith_angle").mean().rename("angle_mean")
        Zenith_STD = (
            imagecol.select("zenith_angle")
            .reduce(ee.Reducer.stdDev())
            .rename("angle_std")
        )
        obs = imagecol.select("LST").count().rename("observation_num").double()
        time = imagecol.select("observation_time").mean().rename("time_mean")
        time_STD = (
            imagecol.select("observation_time")
            .reduce(ee.Reducer.stdDev())
            .multiply(60)
            .rename("time_std")
        )  # in minutes
        uncertainty = imagecol.select("uncertainty").mean().rename("uncertainty_mean")
        uncertainty_max = imagecol.select("uncertainty").max().rename("uncertainty_max")
        return LST_mean.addBands(
            [
                LST_std,
                LST_median,
                LST_p90,
                Zenith,
                Zenith_STD,
                obs,
                time,
                time_STD,
                uncertainty,
                uncertainty_max,
            ]
        ).clip(roi)

    MOD_reduced = getallthedata(image)

    return MOD_reduced


def get_L8(roi, date_start, date_end, months=[6, 8]):
    L8 = (
        ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
        .filterBounds(roi)
        .select(["ST_B10", "QA_PIXEL", "ST_QA"])
        .filterDate(date_start, date_end)
        .filter(ee.Filter.calendarRange(months[0], months[1], "month"))
        .filter(ee.Filter.eq("IMAGE_QUALITY_TIRS", 9))
    )

    return L8


def filter_L8(image):
    def filter(image):
        qa = image.select("QA_PIXEL")
        cloudBitMask = 1 << 3
        cirrusBitMask = 1 << 2
        dilatedcloudBitMask = 1 << 1

        mymask = (
            qa.bitwiseAnd(cloudBitMask)
            .eq(0)
            .And(qa.bitwiseAnd(cirrusBitMask).eq(0))
            .And(qa.bitwiseAnd(dilatedcloudBitMask).eq(0))
        )
        mymask2 = image.select("ST_QA").lt(300)
        image = image.updateMask(mymask)  # .updateMask(mymask2)
        temp = (
            image.select("ST_B10")
            .multiply(0.00341802)
            .add(149)
            .subtract(273.15)
            .rename("LST")
        )

        time_string = ee.String(image.get("SCENE_CENTER_TIME"))
        hour = ee.Number.parse(time_string.slice(0, 2)).add(2)
        minute = ee.Number.parse(time_string.slice(3, 5)).divide(60)
        time = ee.Image.constant(hour.add(minute)).rename("observation_time").double()
        zenith = (
            ee.Image.constant(ee.Number(image.get("SUN_ELEVATION")).add(ee.Number(-90)))
            .rename("zenith_angle")
            .double()
        )
        uncertainty = image.select("ST_QA").multiply(0.01).rename("uncertainty")
        return temp.addBands([time, zenith, uncertainty]).set(
            "system:time_start", image.get("system:time_start")
        )

    # applying the functions to all Data
    L8_clean = image.map(filter)

    return L8_clean


def reduce_L8(image, roi):
    # Reduce Data L8
    def getallthedata(imagecol):
        LST_mean = imagecol.select("LST").mean().rename("LST_mean")
        LST_std = imagecol.select("LST").reduce(ee.Reducer.stdDev()).rename("LST_std")
        LST_median = imagecol.select("LST").median().rename("LST_median")
        LST_p90 = (
            imagecol.select("LST").reduce(ee.Reducer.percentile([90])).rename("LST_p90")
        )
        Zenith = imagecol.select("zenith_angle").mean().rename("angle_mean").double()
        Zenith_STD = (
            imagecol.select("zenith_angle")
            .reduce(ee.Reducer.stdDev())
            .rename("angle_std")
            .double()
        )
        obs = imagecol.select("LST").count().rename("observation_num").double()
        time = imagecol.select("observation_time").mean().rename("time_mean").double()
        time_STD = (
            imagecol.select("observation_time")
            .reduce(ee.Reducer.stdDev())
            .multiply(60)
            .rename("time_std")
            .double()
        )  # in minutes
        uncertainty = imagecol.select("uncertainty").mean().rename("uncertainty_mean")
        uncertainty_max = imagecol.select("uncertainty").max().rename("uncertainty_max")

        return LST_mean.addBands(
            [
                LST_std,
                LST_median,
                LST_p90,
                Zenith,
                Zenith_STD,
                obs,
                time,
                time_STD,
                uncertainty,
                uncertainty_max,
            ]
        ).clip(roi)

    L8_reduced = getallthedata(image)

    return L8_reduced


def export_ee(folder_drive, img, name, px_scale, roi, EPSG_code="4326"):
    task = ee.batch.Export.image.toDrive(
        image=img,
        description=name,
        region=roi,
        scale=px_scale,
        maxPixels=10000000000000,
        folder=folder_drive,
        crs="EPSG:" + EPSG_code,
    )
    task.start()


def process_MODIS(date_start, date_end, date_median, month_summer, roi, folder_drive):
    # Load MODIS data from Earth Engine
    MOD_night = get_MODIS("Night", date_start, date_end, month_summer)
    MOD_day = get_MODIS("Day", date_start, date_end, month_summer)

    print("MOD images:", MOD_night.size().getInfo())

    orig_scale = MOD_day.first().projection().nominalScale().getInfo()
    print("Projection, crs, and crs_transform:", MOD_day.first().projection())
    print("Resolution MOD:", orig_scale)

    # Filter MODIS data
    MOD_night_clean = filter_MODIS(MOD_night, "Night")
    MOD_day_clean = filter_MODIS(MOD_day, "Day")

    print("MOD images:", MOD_night_clean.size().getInfo())

    # Reduce MODIS data
    MOD_night_reduced_start_end = reduce_MODIS(MOD_night_clean, roi)
    MOD_day_reduced_start_end = reduce_MODIS(MOD_day_clean, roi)

    MOD_night_reduced_start_med = reduce_MODIS(
        MOD_night_clean.filterDate(date_start, date_median), roi
    )
    MOD_day_reduced_start_med = reduce_MODIS(
        MOD_day_clean.filterDate(date_start, date_median), roi
    )

    MOD_night_reduced_med_end = reduce_MODIS(
        MOD_night_clean.filterDate(date_median, date_end), roi
    )
    MOD_day_reduced_med_end = reduce_MODIS(
        MOD_day_clean.filterDate(date_median, date_end), roi
    )

    print("final Image:", MOD_day_reduced_start_end.bandTypes().getInfo())

    # Export MODIS data
    export_ee(
        folder_drive,
        MOD_night_reduced_start_end.reproject("EPSG:" + EPSG_code, None, None),
        f"MODIS_night_{date_start}_{date_end}",
        1000,
        roi,
    )
    export_ee(
        folder_drive,
        MOD_day_reduced_start_end.reproject("EPSG:" + EPSG_code, None, None),
        f"MODIS_day_{date_start}_{date_end}",
        1000,
        roi,
    )
    export_ee(
        folder_drive,
        MOD_night_reduced_start_med.reproject("EPSG:" + EPSG_code, None, None),
        f"MODIS_night_{date_start}_med",
        1000,
        roi,
    )
    export_ee(
        folder_drive,
        MOD_day_reduced_start_med.reproject("EPSG:" + EPSG_code, None, None),
        f"MODIS_day_{date_start}_med",
        1000,
        roi,
    )
    export_ee(
        folder_drive,
        MOD_night_reduced_med_end.reproject("EPSG:" + EPSG_code, None, None),
        f"MODIS_night_med_{date_end}",
        1000,
        roi,
    )
    export_ee(
        folder_drive,
        MOD_day_reduced_med_end.reproject("EPSG:" + EPSG_code, None, None),
        f"MODIS_day_med_{date_end}",
        1000,
        roi,
    )

    print("Export MODIS Done.")


def process_L8(bb, date_start, date_end, date_median, month_summer, roi, folder_drive):
    # Load Landsat 8 data from Earth Engine
    L8 = get_L8(bb, date_start, date_end, month_summer)

    orig_scale = L8.first().projection().nominalScale().getInfo()
    print("Resolution L8:", orig_scale)

    print("L8 images:", L8.size().getInfo())

    # Filter Landsat 8 data
    L8_clean = filter_L8(L8)

    # Reduce Landsat 8 data
    L8_reduced_start_end = reduce_L8(L8_clean, roi)
    L8_reduced_start_med = reduce_L8(L8_clean.filterDate(date_start, date_median), roi)
    L8_reduced_med_end = reduce_L8(L8_clean.filterDate(date_median, date_end), roi)

    print("final Image:", L8_reduced_start_end.bandTypes().getInfo())

    # Export Landsat 8 data
    export_ee(
        folder_drive,
        L8_reduced_start_end.reproject("EPSG:" + EPSG_code, None, None),
        f"L8_{date_start}_{date_end}",
        30,
        roi,
    )
    export_ee(
        folder_drive,
        L8_reduced_start_med.reproject("EPSG:" + EPSG_code, None, None),
        f"L8_{date_start}_med",
        30,
        roi,
    )
    export_ee(
        folder_drive,
        L8_reduced_med_end.reproject("EPSG:" + EPSG_code, None, None),
        f"L8_med_{date_end}",
        30,
        roi,
    )

    print("Export L8 Done.")


if __name__ == "__main__":
    # Load configuration file and set paramters
    (
        folder_drive,
        date_start,
        date_end,
        date_median,
        buffer_value,
        bb,
        distric,
        EPSG_code,
    ) = load_config("code/ee_data_config.yml")
    month_summer = [6, 8]

    # Load ROI from Earth Engine
    roi = distric  # or bbox
    roi = roi.geometry().buffer(
        buffer_value
    )  # add buffer zone around the region of interest

    # Process MODIS and Landsat 8 data
    process_MODIS(date_start, date_end, date_median, month_summer, roi, folder_drive)
    process_L8(bb, date_start, date_end, date_median, month_summer, roi, folder_drive)

    # Process the Digital Elevation Model
    # load water (to mask surface water)
    water = ee.ImageCollection("GLCF/GLS_WATER").min().neq(2)
    # elevation
    dem = ee.Image("USGS/GMTED2010").updateMask(water).clip(bb)
    # export elevation data
    export_ee(
        folder_drive,
        dem.reproject("EPSG:" + EPSG_code, None, 1000),
        "dem_GMTED2010_1000",
        1000,
        bb,
    )

    # Process the Night Lights
    # load night lights
    lights = (
        ee.ImageCollection("NOAA/DMSP-OLS/NIGHTTIME_LIGHTS")
        .select("stable_lights")
        .filterDate(date_start, date_end)
        .mean()
        .updateMask(water)
        .clip(bb)
    )
    # export night lights data
    task = ee.batch.Export.image.toDrive(
        image=lights.reproject("EPSG:" + EPSG_code, None, 1000),
        description="night_lights" + str(date_start + "_" + str(date_end) + "_1000"),
        region=roi,
        scale=100,
        maxPixels=10000000000000,
        folder=folder_drive,
    )
    task.start()

    print("Done.")
