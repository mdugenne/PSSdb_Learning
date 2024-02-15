README file for taxa-specific PSSdb products (First created on December 2023):

This directory contains the taxa-specific products (1a: size or biomass spectrum and 1b: derived linear fit parameters) of gridded datasets.
Each table include the following variables:

Taxon or PFT: Taxonomic affiliation according to the World Register of Marine Species database (for Taxon) or to functional groups (for PFT)
Validation_percentage ([0-1]): Percentage of image classification validation
year (yyyy): Year of the temporal bin
month (mm): Month of the temporal bin
ocean: Study area according to the Global Oceans and Seas of the standard georeferenced marine regions website (https://www.marineregions.org/)
latitude (decimal degree): latitude
longitude (decimal degree): longitude
min_depth (meter): Minimum depth
max_depth (meter): Maximum depth
n: Number of individual sub-bins (weekly, half degrees) grouped to compute the average (std) spectrum and parameters
biovolume_size_class (cubic micrometer): Midpoint of the biovolume size bin
normalized_biovolume_mean (cubic micrometer per liter per cubic micrometer): Average of the Normalized biovolume estimates in final bins (month year, 1 degrees)
normalized_biovolume_std (cubic micrometer per liter per cubic micrometer): Standard deviation of the Normalized biovolume estimates in final bins (month year, 1 degrees)
biomass_mid (gram): Midpoint of the biomass bin
range_biomass_bin (gram): Biomass bin width (can be used to calculate total biomass by multiplying with normalized biomass)
normalized_biomass_mean (gram per liter per gram): Average of the Normalized biomass estimates in final bins (month year, 1 degrees)
normalized_biomass_std (gram per liter per gram): Standard deviation of the Normalized biomass estimates in final bins (month year, 1 degrees)
equivalent_circular_diameter_mean (micrometer): Midpoint of the equivalent circular diameter size bin
normalized_abundance_mean (particles per liter per micrometer): Average of the Normalized abundance estimates in final bins (month year, 1 degrees)
normalized_abundance_std (particles per liter per micrometer): Standard deviation of the Normalized abundance estimates in final bins (month year, 1 degrees)
NBSS_slope_mean (liter per cubic micrometer (biovolume) or liter per gram (biomass)): Average of the Normalized biovolume/biomass size spectrum slopes in final bins (month year, 1 degrees)
NBSS_slope_std (liter per cubic micrometer (biovolume) or liter per gram (biomass)): Standard deviation of the Normalized biovolume/biomass size spectrum slopes in final bins (month year, 1 degrees)
NBSS_intercept_mean (per liter): Average of the Normalized biovolume/biomass size spectrum intercepts in final bins (month year, 1 degrees)
NBSS_intercept_std (per liter): Standard deviation of the Normalized biovolume/biomass size spectrum intercepts in final bins (month year, 1 degrees)
NBSS_r2_mean: Average of the linear fit determination coefficients in final bins (month year, 1 degrees)
NBSS_r2_std: Standard deviation of the linear fit determination coefficients in final bins (month year, 1 degrees)
PSD_slope_mean (per liter per micrometer): Average of the particle size distribution slopes in final bins (month year, 1 degrees)
PSD_slope_std (per liter per micrometer): Standard deviation of the particle size distribution slopes in final bins (month year, 1 degrees)
PSD_intercept_mean (particles per liter per micrometer): Average of the particle size distribution intercepts in final bins (month year, 1 degrees)
PSD_intercept_std (particles per liter per micrometer): Standard deviation of the particle size distribution intercepts in final bins (month year, 1 degrees)
PSD_r2_mean: Average of the linear fit determination coefficients in final bins (month year, 1 degrees)
PSD_r2_std: Standard deviation of the linear fit determination coefficients in final bins (month year, 1 degrees)

Contact us: contact@pssdb.net