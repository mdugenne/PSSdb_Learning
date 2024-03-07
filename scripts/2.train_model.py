## Objective: This script uses the xgboost package to predict NBSS parameters as a function of environmental variables

## Functions: All functions required to run this step are included in funcs_predict_NBSS.py

## Modules:
from natsort import natsorted
import xarray as xr
import pandas as pd
import numpy as np # Arrays
import datetime
from pathlib import Path # Handling of path object
# Main functions to download environmental datasets and predict NBSS:
try:
    from funcs_predict_NBSS import *
except:
    from scripts.funcs_predict_NBSS import *
import pint_pandas # Use pip install pint-pandas
PA=pint_pandas.PintArray
from pint import UnitRegistry # Use pip install pint to install
ureg=UnitRegistry()
PQ=ureg.Quantity
ureg.load_definitions(path_to_git / 'data'/ 'units_def.txt')  # This text file is used to define custom units from standard units  (e.g. square_pixel etc.)
full_list_units = list(dict.fromkeys(sorted(dir(ureg))))  # list(dict.fromkeys(sorted(list(np.concatenate([dir(getattr(ureg.sys,system)) for system in dir(ureg.sys)]).flat))))


## Workflow starts here:
path_to_datafile=path_to_git / 'data' /'Model'
# Load allometric look-up table to convert size in biomass
df_allometry=pd.read_excel(path_to_git / 'data' /'plankton_elemental_quotas.xlsx')
df_allometry=df_allometry[df_allometry.Size_proxy=='Biovolume']
# Convert all size units to cubic micrometers and mass units to gram
df_allometry['C_Intercept']=df_allometry['C_Intercept']*(df_allometry.Size_unit.apply(lambda x: PQ(1,'cubic_micrometer').to(x).magnitude)).values/(df_allometry.Elemental_mass_unit.apply(lambda x: PQ(1,'milligram').to(x).magnitude)).values
df_allometry['Size_unit']='cubic_micrometer'
df_allometry['Elemental_mass_unit']='milligram'

# Load dataframe generated on step 0: NBSS parameters per bin per taxon
data=pd.read_csv(path_to_datafile /'Model_input.csv',dtype={'year':int,'month':int})
# Load dataframe generated on step 1 for global predictions: Environmental matrix
df_env_all = xr.open_dataset('{}/Model_environment.nc'.format(str(path_to_git / 'data' / 'Model'))).to_dataframe().reset_index().dropna()#df_environement.dropna().reset_index(drop=True)
df_env_all=pd.merge(df_env_all,df_env_all.drop_duplicates(subset=['Longitude', 'Latitude'], ignore_index=True)[['Longitude', 'Latitude']].reset_index().rename( {'index': 'Group_index'}, axis='columns'),how='left',on=['Longitude', 'Latitude'])
gdf = gpd.GeoDataFrame(df_env_all[['Group_index','Longitude', 'Latitude']].drop_duplicates().dropna(), geometry=gpd.points_from_xy( df_env_all[['Group_index', 'Longitude', 'Latitude']].drop_duplicates().dropna().Longitude,df_env_all[['Group_index','Longitude', 'Latitude']].drop_duplicates().dropna().Latitude))
df_env_all['Study_area'] = pd.merge(df_env_all, gpd.tools.sjoin(gdf, oceans, predicate="within", how='left')[['Group_index', 'name']], how='left',on='Group_index')['name'].astype(str)
df_env_all['Longhurst_area'] = pd.merge(df_env_all, gpd.tools.sjoin(gdf, longhurst, predicate="within", how='left')[['Group_index', 'ProvDescr']], how='left',on='Group_index')['ProvDescr'].astype(str)
df_env_all['Longhurst_subarea'] = df_env_all.Longhurst_area.replace(longhurst_dict)
df_env_all['Hemisphere'] = np.where(df_env_all.Latitude.astype(float)>0,'North','South')
df_env_all['Region']=df_env_all.Study_area+'\n'+df_env_all.Longhurst_subarea
df_env_all['month']=df_env_all.Datetime.str[5:]

# Loop through taxonomic groups to train and fit xgboost model
df_predictions_all=pd.DataFrame({})
df_stat=pd.DataFrame({})
path_to_input=Path(path_to_git/ 'data' / 'Model' )
path_to_input.mkdir(exist_ok=True,parents=True)
for taxa in data['PFT'].dropna().unique():
    taxa=[taxa]
    instrument=data.loc[data.PFT.isin(taxa),'Instrument'].value_counts().index[0]
    list(path_to_input.rglob('xgboost*'))
    if Path("{}/xgboost_{}_{}.json".format(path_to_input,'_'.join(taxa),instrument)).exists():
        model=load_xgbmodel(Path("{}/xgboost_{}_{}.json".format(path_to_input,'_'.join(taxa),instrument)))
        model.evals_result_['validation_0']['logloss'][-1]
    else:
        model, df_stat_summary=nbss_predict(df=data,selection={'PFT':taxa,'Instrument':instrument},var_x=['PSD_slope','Picophytoplankton','Nanophytoplankton','Microphytoplankton','AOU','MLD','CHL','Phosphates','Silicates','Nitrates','Temperature','Iron','PAR','Aerosol_thickness','DUWT001','SLA','FSLE'],var_y=['Intercept','Slope','Min_size','Max_size'],figures_path='') #
        # Concatenate statistics to final table
        df_stat=pd.concat([df_stat,df_stat_summary.assign(taxa='_'.join(taxa))],axis=0).reset_index(drop=True)
    # Display feature importance scores
    df_features=pd.DataFrame(sorted(dict(zip(model.feature_names_in_, model.feature_importances_)).items(), key=lambda x: x[1], reverse=True),columns=['Variable','Score'])
    plot = (ggplot(df_features) +
            geom_col(mapping=aes(x='Variable', y='Score'),fill='black', size=0.01) +
            scale_x_discrete(limits=df_features.Variable.values)+
            labs(x=r'', y=r'Feature importance', fill='') +
            theme_paper).draw(show=True)
    plot.set_size_inches(10, 3)
    plot.savefig(fname='{}/figures/PSSdb_paper/Model_scores_{}.svg'.format(str(path_to_git),"-".join(taxa)), dpi=300, bbox_inches='tight')
    # Generate predictions and re-construct size spectrum using predicted intercept, slope, size range
    df_env = df_env_all[[column for column in model.feature_names_in_ if column in df_env_all.columns]]
    predictions = model.predict(df_env)
    df_predictions=pd.concat([df_env_all, pd.DataFrame({'taxa':"_".join(taxa),'Predicted_intercept': predictions[:, [0]].flatten(), 'Predicted_slope': predictions[:, [1]].flatten(), 'Predicted_min_size': 10**predictions[:, [2]].flatten(), 'Predicted_max_size': 10**predictions[:, [3]].flatten()})], axis=1)
    # Use predicted parameters to reconstruct the size spectrum
    df_predictions['size_class_max']= pd.cut(df_predictions.Predicted_max_size.values,pd.IntervalIndex(df_bins.sizeClasses.values))
    df_predictions.loc[df_predictions.size_class_max.isna(),'size_class_max']=pd.IntervalIndex(df_bins.sizeClasses.values)[pd.IntervalIndex(df_bins.sizeClasses.values).array.contains(np.nanquantile(df_predictions.Predicted_max_size,0.25))].values[0]
    df_predictions['size_class_min']= pd.cut(df_predictions.Predicted_min_size.values,pd.IntervalIndex(df_bins.sizeClasses.values))
    mask=(np.greater_equal.outer(df_bins.sizeClasses.to_numpy(), df_predictions['size_class_min'].to_numpy()) & np.less_equal.outer(df_bins.sizeClasses.to_numpy(), df_predictions['size_class_max'].to_numpy()))
    size_array=np.tile(np.log10((1/6)*np.pi*df_bins.size_class_mid**3).values,(len(df_predictions),1)).T*mask
    size_array[size_array==0]=np.nan
    df_predictions['Predicted_abundance']=np.nansum(10**(df_predictions.Predicted_intercept.to_numpy().reshape(len(df_predictions),1)+df_predictions.Predicted_slope.to_numpy().reshape(len(df_predictions),1)*(size_array.T-np.tile(np.nanmin(size_array,axis=0),(len(df_bins),1)).T)),axis=1) # in cubic meters
    df_predictions['Predicted_meanESD']=(6*np.nansum((10**(df_predictions.Predicted_intercept.to_numpy().reshape(len(df_predictions),1)+df_predictions.Predicted_slope.to_numpy().reshape(len(df_predictions),1)*(size_array.T-np.tile(np.nanmin(size_array,axis=0),(len(df_bins),1)).T))/np.tile(df_predictions.Predicted_abundance.to_numpy(),((len(df_bins),1))).T)*(10**size_array.T),axis=1)/np.pi)**(1/3) # in micrometers
    # Convert biovolume to carbon biomass using the allometric look-up table
    C_intercept,C_slope=df_allometry.loc[df_allometry.Taxon.str.contains(taxa[0].split("_")[1]),['C_Intercept','C_Slope']].values[0] if len(taxa[0].split("_"))>1 else df_allometry.loc[df_allometry.Taxon.str.contains(taxa[0]),['C_Intercept','C_Slope']].values[0]
    df_predictions['Predicted_biomass']=np.nansum(C_intercept*np.power(((np.tile(df_bins.range_size_bin.values,(len(df_predictions),1)))*10**(df_predictions.Predicted_intercept.to_numpy().reshape(len(df_predictions),1)+df_predictions.Predicted_slope.to_numpy().reshape(len(df_predictions),1)*(size_array.T-np.tile(np.nanmin(size_array,axis=0),(len(df_bins),1)).T))),C_slope),axis=1) # in milligram C per cubic meter
    df_predictions.loc[df_predictions.Predicted_biomass>=np.nanmean(df_predictions.Predicted_biomass)+np.nanstd(df_predictions.Predicted_biomass),'Predicted_biomass']=np.nan
    # Alternatively, consider the average biovolume to derive C quota
    df_predictions['Predicted_biomass'] = C_intercept * (((1 / 6) * np.pi * np.nanmean(df_predictions['Predicted_min_size']) ** 3) ** C_slope) * df_predictions[ 'Predicted_abundance']
    if ''.join(taxa)=='Mesophytoplankton': # Remove predictions outside habitat range
        df_predictions.loc[df_predictions.Predicted_intercept<=np.nanquantile(df_predictions.Predicted_intercept,0.8),['Predicted_biomass','Predicted_meanESD']] = np.nan # Pick large quantile because of the log-transformation

    # Generate maps of total abundance (particles per cubic meter), total biomass (mg carbon per cubic meter), and average diameter (micrometer)
    var, kwargs = 'abundance', {'mapping': 'y:Latitude,x:Longitude,fill:Predicted_abundance','labs(fill="Total abundance (particles m$^{-3}$)")':'','scale_fill_gradientn(trans="log10",na_value="#ffffff00",colors={},guide=guide_colorbar(direction="vertical"))'.format(sns.color_palette("blend:#800000ff,#800000ff,#7AB,#f4eed760").as_hex()[::-1]): ''}
    #ggplot_raster(dataframe=df_predictions.dropna(subset=['Predicted_abundance']).groupby(['Latitude','Longitude']).mean().reset_index().astype({'Latitude':float,'Longitude':float}),output_path='',**kwargs, limitsize=False, dpi=600, bbox_inches='tight')
    ggplot_raster(dataframe=df_predictions.dropna(subset=['Predicted_abundance']).groupby(['Latitude','Longitude']).mean().reset_index().astype({'Latitude':float,'Longitude':float}),output_path='{}/GIT/PSSdb_Learning/figures/PSSdb_paper/Map_{}_abundance.svg'.format(str(Path.home()),''.join(taxa)),**kwargs, limitsize=False, dpi=600, bbox_inches='tight')
    var,kwargs = 'biomass',{'mapping': 'y:Latitude,x:Longitude,fill:Predicted_biomass','labs(fill="Total biomass (mg C m$^{-3}$)")':'','scale_fill_gradientn(limits=[{}],trans="sqrt",colors={},guide=guide_colorbar(direction="vertical"))'.format(",".join((np.nanquantile(df_predictions.Predicted_biomass,[0.05,1])).astype(str)),sns.color_palette("blend:#800000ff,#800000ff,#7AB,#f4eed760").as_hex()[::-1]): ''}
    #ggplot_raster(dataframe=df_predictions.dropna(subset=['Predicted_biomass']).groupby(['Latitude','Longitude']).mean().reset_index().astype({'Latitude':float,'Longitude':float}),output_path='',**kwargs, limitsize=False, dpi=600, bbox_inches='tight')
    ggplot_raster(dataframe=df_predictions.dropna(subset=['Predicted_biomass']).groupby(['Latitude','Longitude']).mean().reset_index().astype({'Latitude':float,'Longitude':float}),output_path='{}/GIT/PSSdb_Learning/figures/PSSdb_paper/Map_{}_biomass.svg'.format(str(Path.home()),''.join(taxa)),**kwargs, limitsize=False, dpi=600, bbox_inches='tight')
    var, kwargs = 'esd', {'mapping': 'y:Latitude,x:Longitude,fill:Predicted_meanESD','labs(fill="Average diameter ($\mu$m)")':'','scale_fill_gradientn(trans="sqrt",colors={},na_value="#916f7c60",guide=guide_colorbar(direction="vertical"))'.format(sns.color_palette("blend:#916f7c60,#f4eed760,#f4eed760").as_hex()[::-1]): ''}
    #ggplot_raster(dataframe=df_predictions.dropna(subset=['Predicted_meanESD']).groupby(['Latitude','Longitude']).mean().reset_index().astype({'Latitude':float,'Longitude':float}),output_path='',**kwargs, limitsize=False, dpi=600, bbox_inches='tight')
    ggplot_raster(dataframe=df_predictions.dropna(subset=['Predicted_meanESD']).groupby(['Latitude', 'Longitude']).mean().reset_index().astype({'Latitude': float, 'Longitude': float}), output_path='{}/GIT/PSSdb_Learning/figures/PSSdb_paper/Map_{}_diameter.svg'.format(str(Path.home()),''.join(taxa)), **kwargs,limitsize=False, dpi=600, bbox_inches='tight')


    # Check correlations with iron proxies and save maps
    df_cor = df_predictions.groupby(['Longitude', 'Latitude']).apply(lambda x: pd.DataFrame( {'Taxa': "".join(taxa), 'Abundance': np.nanmean(x.Predicted_abundance), 'Variable': 'Dust', 'value': spearmanr(np.log10(x.Predicted_abundance), np.log10(x.DUWT001))[0], 'pvalue': spearman_correlation(np.log10(x.Predicted_abundance), np.log10(x.DUWT001))}, index=[0])).reset_index()
    df_cor = pd.concat([df_cor, df_predictions.groupby(['Longitude', 'Latitude']).apply(lambda x: pd.DataFrame({'Taxa': "".join(taxa), 'Abundance': np.nanmean(x.Predicted_abundance), 'Variable': 'Aerosol','value': spearmanr(np.log10(x.Predicted_abundance), np.log10(x.Aerosol_thickness))[0], 'pvalue': spearman_correlation(np.log10(x.Predicted_abundance), np.log10(x.Aerosol_thickness))},index=[0])).reset_index()], axis=0)
    df_cor = pd.concat([df_cor, df_predictions.groupby(['Longitude', 'Latitude']).apply(lambda x: pd.DataFrame({'Taxa': "".join(taxa), 'Abundance': np.nanmean(x.Predicted_abundance), 'Variable': 'Iron','value': spearmanr(np.log10(x.Predicted_abundance), np.log10(x.Iron))[0],'pvalue': spearman_correlation(np.log10(x.Predicted_abundance), np.log10(x.Iron))}, index=[0])).reset_index()], axis=0)
    df_cor['p-value'] = df_cor.pvalue.apply(lambda x: ''.join(['*' for t in [.05, .01, .001] if x <= t]))
    df_cor['absolute_value'] = np.where((df_cor.pvalue > 0.05) | (np.log10(df_cor.Abundance) <= np.nanquantile(np.log10(df_cor.Abundance), 0.45)), np.nan, df_cor.value)

    plot = (ggplot(df_cor.query('(absolute_value>0)').dropna(subset=['absolute_value']).groupby(['Longitude', 'Latitude','Variable']).absolute_value.mean().reset_index()) +
            facet_wrap('~Variable',ncol=1)+
            geom_raster(mapping=aes(x='Longitude.astype(float)', y='Latitude.astype(float)', fill='(absolute_value)')) +
            scale_fill_gradientn(trans="log10", colors=sns.color_palette("blend:#800000ff,#7AB,#f4eed760").as_hex()[::-1],guide=guide_colorbar(direction="vertical"), na_value='#ffffff00') +
            labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)', fill=r'Spearman correlation coefficients') +
            coord_cartesian(expand=False) + geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), fill='black', color='black') +
            theme_paper+theme(legend_position='bottom',legend_direction='horizontal',legend_text=element_text(angle=0))).draw(show=True)
    plot.set_size_inches(3, 7)
    plot.savefig(fname='{}/GIT/PSSdb_Learning/figures/PSSdb_paper/Maps_{}_Correlations_iron.svg'.format(str(Path.home()), "_".join(df_cor.Taxa.unique())), dpi=600, bbox_inches='tight')

    df_cor_all = pd.concat([df_cor_all, df_cor], axis=0)


    # Generate a netcdf xarray from Pandas DataFrame
    array=df_predictions[['Latitude','Longitude','Datetime','Predicted_abundance']].astype({'Latitude':float,'Longitude':float}).set_index(['Latitude', 'Longitude','Datetime']).to_xarray()
    array['Latitude'].attrs={'units':'decimal_degrees','description':'Midpoint of the 1 degree latitudinal grid'}
    array['Longitude'].attrs={'units':'decimal_degrees','description':'Midpoint of the 1 degree longitudinal grid'}
    array['Datetime'].attrs={'units':'%Y-%m','description':'Monthly average corresponding to the final resolution of PSSdb products'}
    array['Predicted_abundance'].attrs={'units':'particles per cubic meter','long_name':'Predicted total abundance (spectrum integer)','description':'Total abundance of mesoplankton between 0-200m'}
    array.attrs={'creation':'Dataset generated on March 2024 by M. Dugenne', 'title':'Global predictions of mesoplanktan abundance in the epipelagic layer (0-200m) between 2008-2020', 'summary':'This dataset represents the global predictions of mesoplankton spectral biogeography obtained by boosted gradient regression tree (xgboost) as part of the Pelagic Sie Structure database project (https://pssdb.net)'} # add global attribute metadata
    print(array)
    path_to_output=path_to_input / 'Model_output'
    path_to_output.mkdir(exist_ok=True, parents=True)
    array.to_netcdf('{}/Predictions_{}_abundance.nc'.format(str(path_to_output),"".join(taxa)))
    print('Predictions of {} saved as netcdf file'.format(''.join(taxa)))

