#################################################################################################
#		
# Title:	An integrated dataset of observed national inventories of fires
# 				
# Authors: 	Andrina N. Gincheva, ..., Marco Turco (marco.turco@um.es)
#
#################################################################################################

#################################################################################################
# A. General instructions 
#################################################################################################

This project is designed to be executed with QGIS and with R codes. 
Execute script files in the order they are listed.

Data sources:

- Fire in Australia: 
	New South Wales: https://datasets.seed.nsw.gov.au/dataset/fire-history-wildfires-and-prescribed-burns-1e8b6
	Queensland: https://qldspatial.information.qld.gov.au/catalogue/custom/
	Northern Territory: https://firenorth.org.au/nafi3/downloads/firescars/
	South Australia: https://data.sa.gov.au/data/dataset/fire-history
	Tasmania: http://listdata.thelist.tas.gov.au/opendata/
	Victoria: https://discover.data.vic.gov.au/dataset/
	Western Australia: https://catalogue.data.wa.gov.au/dataset/dbca-fire-history

- Australian regions: https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files#downloads-for-gda2020-digital-boundary-files

- EFFIS data: contact https://effis.jrc.ec.europa.eu/apps/data.request.form/

- European regions (NUTS3): http://ec.europa.eu/eurostat/web/nuts/history

- NFDB dataset for Canada: https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point_large_fires.zip

- NBAC dataset for Canada: https://cwfis.cfs.nrcan.gc.ca/downloads/nbac/nbac_1986_to_2020_20210810.zip

- CONAF dataset for Chile: https://www.conaf.cl/incendios-forestales/incendios-forestales-en-chile/estadisticas-historicas/ (download the document "Estadísticas - Daño (Superficie Afectada) por Incendios Forestales según Mes 1985 - 2022")

- Chile regions shapefile: https://www.bcn.cl/obtienearchivo?id=repositorio/10221/10398/2/Regiones.zip

- FPA-FOD dataset for the U.S.: https://www.fs.usda.gov/rds/archive/catalog/RDS-2013-0009.6

- MTBS dataset for the U.S.: https://mtbs.gov/direct-download

- 1ºx1º grid: https://github.com/SantanderMetGroup/ATLAS/raw/main/reference-grids/land_sea_mask_1degree.nc4



Notes regarding reproducibility:

Script files starting with "1_" in their name are for data preprocessing. 
Most of these script files will NOT run because we do not include the 
raw data because files are simply too large to conveniently share. 
We suggest you run scripts starting with "2_" which directly reproduces the 
results in the paper. 

If you have any questions or wish to express any comment to the authors, please 
contact Dr. Marco Turco at the email indicated above.


#################################################################################################
# B. Description of script files
#################################################################################################

Scripts for dataset development.

- 1_1_x_prepare_data_australia-REGION.R
Read the Australian data for the specific region called "REGION", and aggregate BA values at monthly resolution and onto the 1ºx1º grid. 

- 1_1_8_merge_regional_data_australia.R
Merge the regional data for Australia and export the data.

- 1_2_1_prepare_data_canada_NBAC.R
Read the NBAC dataset, and aggregate BA values at monthly resolution and onto the 1ºx1º grid. Export the data.

- 1_2_2_prepare_data_canada_NFDB.R
Read the NFDB dataset, and aggregate BA values at monthly resolution and onto the 1ºx1º grid. Export the data.

- 1_3_prepare_data_chile.R
Read the CONAF dataset, and aggregate BA values at monthly resolution and onto the 1ºx1º grid. Export the data.

- 1_4_prepare_data_europe.R
Read the EFFIS dataset, and aggregate BA values at monthly resolution and onto the 1ºx1º grid. Export the data.

- 1_5_1_prepare_data_us_fpa_fod.R
Read the FPA-FOD dataset, and aggregate BA values at monthly resolution and onto the 1ºx1º grid. Export the data.

- 1_5_2_prepare_data_us_mtbs.R
Read the MTBS dataset, and aggregate BA values at monthly resolution and onto the 1ºx1º grid. Export the data.

- 1_6_prepare_data_firecci51.R
Read the FIRECCI51 dataset, and filter and aggregate BA values at monthly resolution and onto the 1ºx1º grid.


Scripts for dataset validation.

- 2_1_compare_us_canada.R
Script to analyze data for figures 1 and 3.

- 2_X_compare_REGION.R
Scripts to analyze data for figures 2 and 4.

- 3_X_compare_FWI_REGION.R
Scripts to analyze data for figure 5.


Scripts for plotting the final figures.

- figure_X.R
Script to produce figure X.