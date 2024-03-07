README file for project taxa-specific NBSS estimates (First created on December 2023):

This directory contains the taxa-specific products (1a: size or biomass spectrum and 1b: derived linear fit parameters) of gridded datasets.C Biomass and weight products are derived from published allometric relationship linking taxa-specific biovolume to biomass. Each particle is assigned to a biomass class according to their taxa-specific allometric scaling to compute C biomass products, or to the general relationship linking weight to biovolume published in Lehette and Hernández-León (2009). All parameters are summarised in the allometric look-up table (under ancillary/plankton_elemental_quotas.xlsx, first sheet). Prior to allometric scaling assignment, all taxonomic annotations were standardized to the World Register of Marine Species database, following Neeley et al. (2021) guidelines. The annotations look-up table can be found under ancillary/plankton_annotated_taxonomy.xlsx.
Each table include the following variables:

Taxon or PFT: Taxonomic affiliation according to the World Register of Marine Species database (for Taxon) or to functional groups (for PFT)
Validation_percentage ([0-1]): Percentage of image classification validation
year (yyyy): Year of the temporal bin
month (mm): Month of the temporal bin
ocean: Study area according to the Global Oceans and Seas of the standard geo-referenced marine regions website (https://www.marineregions.org/)
latitude (decimal degree): latitude
longitude (decimal degree): longitude
min_depth (meter): Minimum depth
max_depth (meter): Maximum depth
n: Number of individual sub-bins (weekly, half degree) grouped to compute the average (std) size spectrum and parameters
biovolume_size_class (cubic micrometer): Midpoint of the biovolume size bin
normalized_biovolume_mean (cubic micrometer per liter per cubic micrometer): Average of the Normalized biovolume estimates in final bins (month year, 1x1 degree)
normalized_biovolume_std (cubic micrometer per liter per cubic micrometer): Standard deviation of the Normalized biovolume estimates in final bins (month year, 1x1 degree)
biomass_mid (gram): Midpoint of the biomass bin
range_biomass_bin (gram): Biomass bin width (can be used to calculate total biomass by multiplying with normalized biomass)
normalized_biomass_mean (gram per liter per gram): Average of the Normalized biomass estimates in final bins (month year, 1x1 degree)
normalized_biomass_std (gram per liter per gram): Standard deviation of the Normalized biomass estimates in final bins (month year, 1x1 degree)
equivalent_circular_diameter_mean (micrometer): Midpoint of the equivalent circular diameter size bin
normalized_abundance_mean (particles per liter per micrometer): Average of the Normalized abundance estimates in final bins (month year, 1x1 degree)
normalized_abundance_std (particles per liter per micrometer): Standard deviation of the Normalized abundance estimates in final bins (month year, 1x1 degree)
NBSS_slope_mean (liter per cubic micrometer (biovolume) or liter per gram (biomass)): Average of the Normalized biovolume/biomass size spectrum slopes in final bins (month year, 1x1 degree)
NBSS_slope_std (liter per cubic micrometer (biovolume) or liter per gram (biomass)): Standard deviation of the Normalized biovolume/biomass size spectrum slopes in final bins (month year, 1x1 degree)
NBSS_intercept_mean (log_10 of per liter): Average of the Normalized biovolume/biomass size spectrum intercepts in final bins (month year, 1x1 degree)
NBSS_intercept_std (log_10 of per liter): Standard deviation of the Normalized biovolume/biomass size spectrum intercepts in final bins (month year, 1x1 degree)
NBSS_r2_mean: Average of the linear fit determination coefficients in final bins (month year, 1x1 degree)
NBSS_r2_std: Standard deviation of the linear fit determination coefficients in final bins (month year, 1x1 degree)
PSD_slope_mean (log_10 of per liter per square micrometer): Average of the particle size distribution slopes in final bins (month year, 1x1 degree)
PSD_slope_std (log_10 of per liter per square micrometer): Standard deviation of the particle size distribution slopes in final bins (month year, 1x1 degree)
PSD_intercept_mean (log_10 of particles per liter per micrometer): Average of the particle size distribution intercepts in final bins (month year, 1x1 degree)
PSD_intercept_std (log_10 of particles per liter per micrometer): Standard deviation of the particle size distribution intercepts in final bins (month year, 1x1 degree)
PSD_r2_mean: Average of the linear fit determination coefficients in final bins (month year, 1x1 degree)
PSD_r2_std: Standard deviation of the linear fit determination coefficients in final bins (month year, 1x1 degree)
QC_3std_dev (1 if data is flagged, 0 otherwise): Flag identifying spectra whose slopes are outside the mean ± 3 standard deviations
QC_min_n_size_bins (1 if data is flagged, 0 otherwise): Flag identifying spectra which only include 4 or less non-empty size classes
QC_R2 (1 if data is flagged, 0 otherwise): Flag identifying spectra whose linear fits yield a coefficient of determination below 0.8

Contact us: contact@pssdb.net
References:
Lehette, P. & Hernández-León, S. Zooplankton biomass estimation from digitized images: a comparison between subtropical and Antarctic organisms. Limnol. Oceanogr. Methods 7, 304–308 (2009).
Neeley, A. et al. Standards and practices for reporting plankton and other particle observations from images. (2021)