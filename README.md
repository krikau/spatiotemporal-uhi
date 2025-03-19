# About
- Calculation of temperature anomalies using land surface temperatures (LST), air temperature (Ta), and thermal comfort parameters such as Humidex.  
- Spatiotemporal analysis of heat extremes in Hesse, Germany, considering local climate zones (LCZ) and population density.  

## Input Data

All datasets, except for the PET data, are openly and freely available. Therefore, they are not included in this repository. To ensure the traceability of results, a CSV file containing the combined datasets required for reproducing the analysis is provided. For the PET dataset, dummy values are included in place of the original data.  

### Heat Metrics

- **Physiological Equivalent Temperature (PET)**: Distribution of the datset is restricted, however access can be requested at the "Hessisches Ministerium für Wirtschaft, Energie, Verkehr und Wohnen" (https://landesplanung.hessen.de/klima/landesweite-klimaanalyse/karten-simulationsergebnisse)
- **Land Surface Temperatures (LST)**:  LST from MODIS [Aqua](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MYD11A1) and [Terra](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD11A1) was downloaded using Google Earth Engine.   
- **Air Temperature (Ta) and Relative Humidity (RH)**: The [HOSTRADA dataset](https://opendata.dwd.de/climate_environment/CDC/grids_germany/hourly/hostrada/) contains meteorological information for Germany.
- **Humidex**: Was calculated based on the HOSTRADA dataset.

### Additional Spatial Data
- **Local Climate Zones (LCZ)**: The [global map of LCZ](https://doi.org/10.5281/zenodo.6364593) was used.
- **Urban / Rural Distinction**: Based on the [CORINE Land Cover 2018](https://doi.org/10.2909/960998c1-1870-4e82-8051-6485205ebbac).
- **Population Density**: Based on the [census data](https://www.zensus2022.de/EN/What-is-the-census/grid_cells_results_2011.html) from 2011.
- **Elevation Data**: [Copernicus Global Digital Elevation Model]( https://doi.org/10.5069/G9028PQB) from the European Space Agency (2024) downloaded via Google Earth Engine.
- **Region of Interest (ROI)**: Additional data such as the region of interest are available from the Authorative Topographic-Cartographic Information System ([ATKIS](https://www.adv-online.de/Products/Geotopography/ATKIS/))

## Workflow

- The data paths in the Jupyter notebooks ([01_processing_creation.ipynb](./code/01_processing_creation.ipynb), [02_processing_combination.ipynb](./code/02_processing_combination.ipynb)) will need to be adjusted to match your specific configuration. The notebook [03_analysis.ipynb](./code/03_analysis.ipynb) should run out of the box, although the PET values were replaced by random values.

### 0. Download the data. 
- Code for the data download from Google Earth Engine is provided in the notebook [ee_data.ipynb](./code/ee/ee_data.ipynb) change configuration parameters in [ee_data_config.yaml](./code/ee/ee_data_config.yaml) if needed. Everything else was manually downloaded. Although the ROI is Hesse, Germany, a 100 km buffer zone is required for calculating the temperature anomalies (ΔT). Therefore temperature values from areas outside of the ROI need to be included as well.
- All datasets are resampled to a spatial resolution of 1 km (using neares neighbour interpolation); ensure uniform image dimensions (width and height) for further analysis. This was done using [QGIS](https://www.qgis.org/).


### 1. [01_processing_creation.ipynb](./code/01_processing_creation.ipynb)
- Processing of the raw data from MODIS and HOSTRADA.
- Caculation of derived datasets such as thermal comfort paramters and extrem heat events.
- Caculation of temperature anomalies.

### 2. [02_processing_combination.ipynb](./code/02_processing_combination.ipynb)
- Combining the datasets in different datframes for quick analysis.

### 3. [03_analysis.ipynb](./code/03_analysis.ipynb)
- Creation of figures. 

### 4. [04_extreme.ipynb](./code/04_extreme.ipynb)
- heat extreme events such as heat waves and (very) hot days in relation to LST and Ta