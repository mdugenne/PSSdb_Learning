## Objective: This script generates multiple plots associated to model fits and predictions, including maps and climatologies

## Modules
import pandas as pd
import numpy as np # Arrays
import datetime
import xarray as xr
from pathlib import Path # Handling of path object
# Main functions to download environmental datasets and predict NBSS:
try:
    from funcs_predict_NBSS import *
except:
    from scripts.funcs_predict_NBSS import *

## Workflow starts here
path_to_input=path_to_git / 'data' / 'Model' / 'Model_output'
path_to_files=list(path_to_input.glob('Predictions_*_abundance.nc'))
df_predictions_all=pd.concat(map(lambda path:xr.open_dataset(path).to_dataframe().reset_index().assign(taxa=path.name.split("_")[1]),path_to_files))
df_predictions_all=pd.merge(df_predictions_all,df_predictions_all.drop_duplicates(subset=['Longitude', 'Latitude'], ignore_index=True)[['Longitude', 'Latitude']].reset_index().rename( {'index': 'Group_index'}, axis='columns'),how='left',on=['Longitude', 'Latitude'])
gdf = gpd.GeoDataFrame(df_predictions_all[['Group_index','Longitude', 'Latitude']].drop_duplicates().dropna(), geometry=gpd.points_from_xy( df_predictions_all[['Group_index', 'Longitude', 'Latitude']].drop_duplicates().dropna().Longitude,df_predictions_all[['Group_index','Longitude', 'Latitude']].drop_duplicates().dropna().Latitude))
df_predictions_all['Study_area'] = pd.merge(df_predictions_all, gpd.tools.sjoin(gdf, oceans, predicate="within", how='left')[['Group_index', 'name']], how='left',on='Group_index')['name'].astype(str)
df_predictions_all['Longhurst_area'] = pd.merge(df_predictions_all, gpd.tools.sjoin(gdf, longhurst, predicate="within", how='left')[['Group_index', 'ProvDescr']], how='left',on='Group_index')['ProvDescr'].astype(str)
df_predictions_all['Longhurst_subarea'] = df_predictions_all.Longhurst_area.replace(longhurst_dict)
df_predictions_all['Hemisphere'] = np.where(df_predictions_all.Latitude.astype(float)>0,'North','South')
df_predictions_all['Region']=df_predictions_all.Study_area+'\n'+df_predictions_all.Longhurst_subarea
df_predictions_all['month']=df_predictions_all.Datetime.str[5:].astype(float)


# Plot climatologies
df_predictions=df_predictions_all[(df_predictions_all.Region.isin( ['Arctic Ocean\nPolar','North Pacific Ocean\nTemperate-Subpolar', 'North Atlantic Ocean\nSubtropical', 'South Atlantic Ocean\nSubtropical', 'North Pacific Ocean\nSubtropical', 'South Pacific Ocean\nSubtropical','Southern Ocean\nPolar','South Pacific Ocean\nEquatorial-Upwellings']))]
for taxa in df_predictions.taxa.unique():
    plot = (ggplot(df_predictions.query('taxa=="{}"'.format(taxa))) +
                facet_wrap('~Region', scales='free') +
                stat_summary(mapping=aes(x='month', y='Predicted_abundance'), geom='pointrange', fun_data='median_hilow',fun_args={'confidence_interval': 0.5}) +
                scale_y_log10() + scale_x_continuous(breaks=np.arange(1, 13, 1)) +
                labs(x=r'', y=r'Total abundance (particles m$^{-3}$)') +
                theme_paper).draw(show=True)
    plot.set_size_inches(10,7)
    plot.savefig(fname='{}/GIT/PSSdb_Learning/figures/PSSdb_paper/Climatologies_{}.svg'.format(str(Path.home()), taxa),dpi=600, bbox_inches='tight')

# Plot proportions - Use R ggnewscale to superimpose plots directly
colors_dict={'Crustaceans':"blend:#916f7c73,#ffffff99,#ffffff99",'Rhizaria':"blend:#d0072eff,#ffffff99,#ffffff99",'Mesophytoplankton':"blend:#5f6fbcff,#ffffff99,#ffffff99",'Carnivore':"blend:#916f7c73,#ffffff99,#ffffff99",'Omnivore':"blend:#aaaaff73,#ffffff99,#ffffff99",'Malacostraca':"blend:#d0072eff,#ffffff99,#ffffff99","Copepoda":"blend:#d4c62aff,#ffffff99,#ffffff99,#ffffff99",'Ostracoda':"blend:#536c67ff,#ffffff99,#ffffff99",'Acantharia':"blend:#d38d5fff,#ffffff99,#ffffff99",'Collodaria':"blend:#d3145fff,#ffffff99,#ffffff99",'Phaeodaria':"blend:#b3b30aff,#ffffff99,#ffffff99",'Foraminifera':"blend:#5f6fbcff,#ffffff99,#ffffff99",'puff':"blend:#916f7cff,#ffffff99,#ffffff99",'tuff':"blend:#42dbbdff,#ffffff99,#ffffff99"}
df_predictions=df_predictions_all.query('taxa.isin(["Carnivore","Omnivore"])').groupby(['Longitude','Latitude','taxa']).Predicted_abundance.mean().reset_index()
df_predictions=df_predictions_all.query('taxa.isin(["Acantharia","Collodaria","Phaeodaria","Foraminifera"])').groupby(['Longitude','Latitude','taxa']).Predicted_abundance.mean().reset_index()
df_predictions=df_predictions_all.query('taxa.isin(["Copepoda","Malacostraca","Ostracoda"])').groupby(['Longitude','Latitude','taxa']).Predicted_abundance.mean().reset_index()
df_predictions=df_predictions_all.query('taxa.isin(["puff","tuff"])').groupby(['Longitude','Latitude','taxa']).Predicted_abundance.mean().reset_index()
df_predictions.loc[((df_predictions.Predicted_abundance <= np.nanquantile(df_predictions.Predicted_abundance, 0.85)) | (np.abs(df_predictions.Latitude)>35)), 'Predicted_abundance'] = 0  # Pick large quantile because of the log-transformation

groups=df_predictions.taxa.unique()
df_predictions=df_predictions.pivot_table(values='Predicted_abundance',index=['Longitude','Latitude'],columns='taxa').reset_index()

# Convert to proportions
df_predictions[groups]=df_predictions[groups]/np.nansum(df_predictions[groups],axis=1).reshape(len(df_predictions),1)

for i,group in pd.Series(groups).items():
    plot= (ggplot(df_predictions) +
             geom_raster(mapping=aes(x='Longitude', y='Latitude', fill='{}'.format(group))) +
             scale_fill_gradientn(trans="log10", colors=sns.color_palette(colors_dict[group]).as_hex()[::-1], guide=guide_colorbar(direction="vertical"), na_value='#ffffff00') +
             labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)', fill=r'Proportion') +
             coord_cartesian(expand=False)+geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'),fill='black',color='black') +
             theme_paper).draw(show=True)
             plot.set_size_inches(10, 7)
             plot.savefig(fname='{}/GIT/PSSdb_Learning/figures/PSSdb_paper/Map_ratio_{}.svg'.format(str(Path.home()), group),dpi=600, bbox_inches='tight')