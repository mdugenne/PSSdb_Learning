## Objective: This script downloads and matches NBSS observations to common environmental variables that will be used for prediction

## Functions: All functions required to run this step are included in funcs_predict_NBSS.py

## Modules:
from natsort import natsorted
import pandas as pd
import numpy as np # Arrays
import datetime
from pathlib import Path # Handling of path object
import xarray
import calendar

# Main functions to download environmental datasets and predict NBSS:
try:
    from funcs_predict_NBSS import *
except:
    from scripts.funcs_predict_NBSS import *

## Workflow starts here:
path_to_datafile=path_to_git / 'data' /'Model'
# Load dataframe generated on step 0: NBSS parameters per bin per taxon
nbss_summary=pd.read_csv(path_to_datafile /'Model_input.csv',dtype={'year':int,'month':int}).rename(columns={'latitude':'Latitude','longitude':'Longitude'})
nbss_summary['year']=nbss_summary.year.astype(int).astype(str)
nbss_summary['month']=nbss_summary.month.astype(str).str.zfill(2)
nbss_summary['Date_bin']=nbss_summary[['year','month']].agg('-'.join, axis=1)
# Download environmental variables if missing
for year_bin in natsorted(nbss_summary.year.str.split('-').str[0].unique())[0:1]:
    with tqdm(desc='Downloading environmental data. Please wait\n', total=len(nbss_summary.year.str.split('-').str[0].unique()[0:1]), bar_format='{desc}{bar}', position=0, leave=True) as bar:
        try:
            WOA_download(year_bin=year_bin, output_path=path_to_git / 'data' / 'Environment' / 'WOA', url='https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/', grid='1.00',output_format='netcdf')
        except:
            print('\nDownloading World Ocean Atlas data failed. Skipping\n')

        try:
            Copernicus_download(year_bin=year_bin, username=cfg_pw['Copernicus_user'], password=cfg_pw['Copernicus_pass'], output_path=path_to_git / 'data' / 'Environment' / 'Copernicus')
        except:
            print('\nDownloading Copernicus data failed. Skipping\n')

        try:
            NASA_download(year_bin=year_bin, output_path=path_to_git / 'data' / 'Environment' / 'NASA', appkey=cfg_pw['NASA_appkey'])
        except:
            print('\nDownloading NASA Ocean color data failed. Skipping\n')

        try:
            AVISO_download(year_bin, username=cfg_pw['AVISO_user'], password=cfg_pw['AVISO_pass'], output_path=path_to_git / 'data' / 'Environment' / 'AVISO')
        except:
            print('\nDownloading AVISO data failed. Skipping\n')
        try:
            PSD_download(year_bin, url='https://doi.pangaea.de/10.1594/PANGAEA.939863?format=html#download', output_path=path_to_git / 'data' / 'Environment' / 'PSD')
        except:
             print('\nDownloading PSD data failed. Skipping\n')

        percent = 100*(bar.n/len(nbss_summary.year.str.split('-').str[0].unique()))
        bar.set_description("Downloading environmental data for year {} {}%".format(year_bin,percent), refresh=True)
        bar.update(n=1)

# Merge environmental variables to NBSS parameters and generate associated maps
if len(list(Path(path_to_git / 'data' / 'Environment' /'WOA').expanduser().glob('*.parquet')))==0:
    params_woa_dict = {'ocean_mixed_layer_thickness': 'MLD', 'sea_water_electrical_conductivity': 'Conductivity','sea_water_salinity': 'Salinity','mole_concentration_of_dissolved_molecular_oxygen_in_sea_water': 'Oxygen','Apparent_Oxygen_Utilization': 'AOU','moles_concentration_of_phosphate_in_sea_water': 'Phosphates','moles_concentration_of_silicate_in_sea_water': 'Silicates', 'moles_concentration_of_nitrate_in_sea_water': 'Nitrates', 'sea_water_temperature': 'Temperature'}
    df_woa=pd.concat(map (lambda path:(df:=pd.read_parquet(path),df.rename(columns={'latitude':'Latitude','longitude':'Longitude','month':'Month','value_surface':params_woa_dict[df.variable.unique()[0]]}).drop_duplicates(['Latitude','Longitude','Month',params_woa_dict[df.variable.unique()[0]]])[['Latitude','Longitude','Month',params_woa_dict[df.variable.unique()[0]]]].reset_index(drop=True).astype({'Latitude':str,'Longitude':str}))[-1],list(Path(path_to_git / 'data' / 'Environment'/'WOA').expanduser().rglob('*.parquet'))),axis=1).reset_index(drop=True)
    df_woa=df_woa.loc[:,~df_woa.columns.duplicated()]
    df_woa.to_parquet(path_to_git / 'data' / 'Environment' /'WOA' / 'WOA_all.parquet',index=False)
    [file.unlink(missing_ok=True) for file in list(Path(path_to_git / 'data' / 'Environment' /'WOA').expanduser().rglob('*.parquet'))[1:]]  # Delete original files to save storage space
else:
    df_woa=pd.read_parquet(path_to_git / 'data' / 'Environment' / 'WOA' / 'WOA_all.parquet')
df_woa.columns
# var,kwargs = 'Temperature',{'mapping': 'y:Latitude,x:Longitude,fill:Temperature',"scale_fill_gradientn(colors={},guide=guide_colorbar(direction='vertical'))".format(list(sns.color_palette("blend:#7AB,#f4eed760,#e9afafff,#aa0044ff").as_hex())): ''}
# var,kwargs = 'Nitrates',{'mapping': 'y:Latitude,x:Longitude,fill:Nitrates/1000','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("Purples",15).as_hex()))[::-1]): ''}
# var,kwargs = 'Silicates',{'mapping': 'y:Latitude,x:Longitude,fill:Silicates/1000','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("Reds",15).as_hex()))[::-1]): ''}
# var,kwargs = 'Phosphates',{'mapping': 'y:Latitude,x:Longitude,fill:Phosphates/1000','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("Reds",15).as_hex()))[::-1]): ''}
# var,kwargs = 'Oxygen',{'mapping': 'y:Latitude,x:Longitude,fill:Oxygen','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("Blues",15).as_hex()))[::-1]): ''}
# var,kwargs = 'AOU',{'mapping': 'y:Latitude,x:Longitude,fill:(AOU)','scale_fill_gradientn(colors={},limits=[-10,50],guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("blend:#aa0044ff,#5A9,#ffffffff").as_hex()))): ''}
# var,kwargs = 'MLD',{'mapping': 'y:Latitude,x:Longitude,fill:MLD','scale_fill_gradientn(trans="sqrt",colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("OrRd",15).as_hex()))[::-1]): ''}
# ggplot_raster(dataframe=df_woa.astype({'Latitude':float,'Longitude':float}),output_path='',**kwargs)
# ggplot_raster(dataframe=df_woa[df_woa.Month=='January'].astype({'Latitude':float,'Longitude':float}),output_path='{}/figures/PSSdb_paper/Map_{}.svg'.format(str(path_to_git),var),**kwargs)

nbss_summary=nbss_summary.assign(Month=pd.to_datetime(nbss_summary.Date_bin,format='%Y-%m').dt.strftime('%B'),month=pd.to_datetime(nbss_summary.Date_bin,format='%Y-%m').dt.strftime('%m').astype(float))
nbss_summary=pd.merge(nbss_summary,df_woa,how='left',on=['Latitude','Longitude','Month'])

## Read Copernicus summary file and merge to nbss_summary
if len(list(Path(path_to_git / 'data' / 'Environment'/'Copernicus').expanduser().glob('*.parquet')))==1:
    df_copernicus=pd.read_parquet(list(Path(path_to_git / 'data' / 'Environment'/'Copernicus').expanduser().glob('*.parquet'))[0])
else:
    df_copernicus_sla=pd.concat(map (lambda path:(df:=pd.read_parquet(path),df:=df.assign(lon=(df.lon + 180) % 360 - 180),df.rename(columns={'lat':'Latitude','lon':'Longitude','sla':'SLA'})[['Latitude','Longitude','Date_bin','SLA']].reset_index(drop=True).astype({'Latitude':str,'Longitude':str}))[-1],natsorted(list(Path(path_to_git / 'data' / 'Environment'/'Copernicus'/'SLA').expanduser().rglob('*.parquet')))),axis=0).reset_index(drop=True).sort_values(['Latitude','Longitude','Date_bin'])
    df_copernicus_chl=pd.concat(map (lambda path:(df:=pd.read_parquet(path),df:=df.assign(lon=(df.lon + 180) % 360 - 180),df.rename(columns={'lat':'Latitude','lon':'Longitude'})[['Latitude','Longitude','Date_bin','CHL']].reset_index(drop=True).astype({'Latitude':str,'Longitude':str}))[-1],natsorted(list(Path(path_to_git / 'data' / 'Environment'/'Copernicus'/'CHL').expanduser().rglob('*.parquet')))),axis=0).reset_index(drop=True).sort_values(['Latitude','Longitude','Date_bin'])
    df_copernicus_bgc=pd.concat(map (lambda path:(df:=pd.read_parquet(path),df:=df.assign(lon=(df.lon + 180) % 360 - 180),df.rename(columns={'lat':'Latitude','lon':'Longitude','fe_surface':'Iron'})[['Latitude','Longitude','Date_bin','Iron']].reset_index(drop=True).astype({'Latitude':str,'Longitude':str}))[-1],natsorted(list(Path(path_to_git / 'data' / 'Environment'/'Copernicus'/'BGC').expanduser().rglob('*.parquet')))),axis=0).reset_index(drop=True).sort_values(['Latitude','Longitude','Date_bin'])
    df_copernicus = reduce(lambda left,right:pd.merge(left,right,on=['Latitude','Longitude','Date_bin']), [df_copernicus_sla,df_copernicus_chl,df_copernicus_bgc]).reset_index(drop=True).sort_values(['Latitude','Longitude','Date_bin'])
    if len(list(Path(path_to_git / 'data' / 'Environment'/'Copernicus').expanduser().glob('*.parquet')))==0:
        df_copernicus.to_parquet(path_to_git / 'data' / 'Environment'/'Copernicus' / 'Copernicus_all.parquet',index=False)
    [file.unlink(missing_ok=True) for file in list(Path(path_to_git / 'data' / 'Environment' / 'Copernicus').expanduser().rglob('*.parquet'))[ 1:]]  # Delete original files to save storage space

df_copernicus.columns
# var,kwargs = 'CHL',{'mapping': 'y:Latitude,x:Longitude,fill:CHL','scale_fill_gradientn(trans="log10",colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("BuGn",15).as_hex()))[::-1]): ''}
# var,kwargs ='SLA', {'mapping': 'y:Latitude,x:Longitude,fill:SLA-np.mean(SLA)','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("RdBu",15).as_hex()))): ''}
# var,kwargs = 'Iron',{'mapping': 'y:Latitude,x:Longitude,fill:Iron','scale_fill_gradientn(trans="log10",colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("YlOrBr",15).as_hex()))[::-1]): ''}
# ggplot_raster(dataframe=df_copernicus.dropna(subset=[var])[df_copernicus.Date_bin==df_copernicus.Date_bin.unique()[0]].astype({'Latitude':float,'Longitude':float}),output_path='',**kwargs)
# ggplot_raster(dataframe=df_copernicus.dropna(subset=[var])[df_copernicus.Date_bin==df_copernicus.Date_bin.unique()[0]].astype({'Latitude':float,'Longitude':float}),output_path='{}/figures/PSSdb_paper/Map_{}.svg'.format(str(path_to_git),var),**kwargs)

nbss_summary=pd.merge(nbss_summary,df_copernicus.astype({'Latitude':float,'Longitude':float}),how='left',on=['Latitude','Longitude','Date_bin'])

## Read NASA summary file and merge to nbss_summary
if len(list(Path(path_to_git / 'data' / 'Environment'/'NASA').expanduser().glob('*.parquet')))==0:
    df_nasa_par=pd.concat(map(lambda path: (df := pd.read_parquet(path), df.rename(columns={'lat': 'Latitude', 'lon': 'Longitude'}).astype( {'Latitude': str, 'Longitude': str}))[-1], list(Path(path_to_git / 'data' / 'Environment' / 'NASA' / 'par').expanduser().rglob('*.parquet'))), axis=0).reset_index(drop=True).rename(columns={'par': 'PAR'})
    df_nasa_aerosol=pd.concat(map(lambda path: (df := pd.read_parquet(path), df.rename(columns={'lat': 'Latitude', 'lon': 'Longitude'}).astype( {'Latitude': str, 'Longitude': str}))[-1], list(Path(path_to_git / 'data' / 'Environment' / 'NASA' / 'aot_869').expanduser().rglob('*.parquet'))), axis=0).reset_index(drop=True).rename(columns={'aot_869': 'Aerosol_thickness'})
    df_nasa = reduce(lambda left, right: pd.merge(left.drop(columns=['units', 'datetime']),right.drop(columns=['units', 'datetime']), on=['Latitude', 'Longitude', 'Date_bin'], how='outer'),[df_nasa_par, df_nasa_aerosol]).reset_index(drop=True).sort_values(['Date_bin', 'Latitude', 'Longitude'])
    df_nasa_diagnostic = pd.concat(map(lambda path: (df := pd.read_parquet(path),df.rename(columns={'lat': 'Latitude', 'lon': 'Longitude'}).astype( {'Latitude': str, 'Longitude': str}))[-1], list( Path(path_to_git / 'data' / 'Environment' / 'NASA' / 'aerosol_diagnostic').expanduser().rglob('*.parquet'))),axis=0).reset_index(drop=True)
    df_nasa = reduce(lambda left, right: pd.merge(left, right.drop(columns=['units', 'datetime', 'time']), on=['Latitude', 'Longitude', 'Date_bin'], how='outer'),[df_nasa, df_nasa_diagnostic]).reset_index(drop=True).sort_values( ['Date_bin', 'Latitude', 'Longitude'])
    df_nasa = df_nasa.assign(Month=lambda x: pd.to_datetime(x.Date_bin, format='%Y-%m').dt.strftime('%B'), month=lambda x: pd.to_datetime(x.Date_bin, format='%Y-%m').dt.strftime('%m').astype(float))
    df_nasa.to_parquet(path_to_git / 'Ancillary' / 'NASA' / 'NASA_all.parquet', index=False)
    [file.unlink(missing_ok=True) for file in list(Path(path_to_git / 'data' / 'Environment' / 'NASA').expanduser().rglob('*.parquet'))[1:]]  # Delete original files to save storage space
else:
    df_nasa=pd.read_parquet(path_to_git / 'data' / 'Environment' / 'NASA' / 'NASA_all.parquet')
df_nasa.columns
# var,kwargs = 'Aerosol_thickness',{'mapping': 'y:Latitude,x:Longitude,fill:Aerosol_thickness','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("OrRd",15).as_hex()))[::-1]): ''}
# var,kwargs = 'PAR',{'mapping': 'y:Latitude,x:Longitude,fill:PAR','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("YlOrRd",15).as_hex()))[::-1]): ''}
# var,kwargs = 'DUWT001',{'mapping': 'y:Latitude,x:Longitude,fill:1e10*DUWT001','scale_fill_gradientn(trans="sqrt",colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("PuOr",15).as_hex()))[5:]): ''}
# ggplot_raster(dataframe=df_nasa.dropna(subset=[var])[df_nasa.Date_bin==df_nasa.Date_bin.unique()[1]].astype({'Latitude':float,'Longitude':float}),output_path='',**kwargs)
# ggplot_raster(dataframe=df_nasa.dropna(subset=[var])[df_nasa.Date_bin==df_nasa.Date_bin.unique()[1]].astype({'Latitude':float,'Longitude':float}),output_path='{}/figures/PSSdb_paper/Map_{}.svg'.format(str(path_to_git),var),**kwargs)

nbss_summary=pd.merge(nbss_summary,df_nasa.astype({'Latitude':float,'Longitude':float}),how='left',on=['Latitude','Longitude','Date_bin'])

## Read NASA PSD summary file and merge to nbss_summary
if len(list(Path(path_to_git / 'data' / 'Environment' /'PSD').expanduser().glob('*.parquet')))==0:
    path_parquet=list(Path(ppath_to_git / 'data' / 'Environment' /'PSD').expanduser().rglob('*.parquet'))
    path_parquet.remove(list(Path(path_to_git / 'data' / 'Environment' /'PSD').expanduser().glob('*.parquet'))[0])
    df_psd=pd.concat(map (lambda path:(df:=pd.read_parquet(path),df.drop(columns=['Latitude','Longitude']).rename(columns={'lat':'Latitude','lon':'Longitude','Picophytoplankton_carbon_fraction':'Picophytoplankton','Nanophytoplankton_carbon_fraction':'Nanophytoplankton','Microphytoplankton_carbon_fraction':'Microphytoplankton'}).astype({'Latitude':str,'Longitude':str}))[-1],natsorted(path_parquet)),axis=0).reset_index(drop=True)
    df_psd.to_parquet(path_to_git / 'data' / 'Environment' /'PSD' / 'PSD_all.parquet',index=False)
    [file.unlink(missing_ok=True) for file in list(Path(path_to_git / 'data' / 'Environment' /'PSD').expanduser().rglob('*.parquet'))[1:]]  # Delete original files to save storage space

else:
    df_psd=pd.read_parquet(path_to_git / 'data' / 'Environment' / 'PSD' / 'PSD_all.parquet').rename(columns={'Picophytoplankton_carbon_fraction':'Picophytoplankton','Nanophytoplankton_carbon_fraction':'Nanophytoplankton','Microphytoplankton_carbon_fraction':'Microphytoplankton'})
    df_psd=df_psd.astype({'Latitude':'float','Longitude':'float'}).round({'Latitude':1,'Longitude':1}).copy()
df_psd.columns
# var,kwargs = 'PSD_slope',{'mapping': 'y:Latitude,x:Longitude,fill:PSD_slope','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("blend:#7AB,#EDA",15).as_hex()))[::-1]): ''}
# var,kwargs = 'Microphytoplankton',{'mapping': 'y:Latitude,x:Longitude,fill:Microphytoplankton','scale_fill_gradientn(trans="sqrt",colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("YlOrRd",15).as_hex()))[::-1]): ''}
# var,kwargs = 'Nanophytoplankton',{'mapping': 'y:Latitude,x:Longitude,fill:Nanophytoplankton','scale_fill_gradientn(trans="sqrt",colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("PuOr",15).as_hex()))[5:]): ''}
# var,kwargs = 'Picophytoplankton',{'mapping': 'y:Latitude,x:Longitude,fill:Picophytoplankton','scale_fill_gradientn(trans="sqrt",colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("PuOr",15).as_hex()))[5:]): ''}
# ggplot_raster(dataframe=df_psd.dropna(subset=[var])[df_psd.Date_bin==df_psd.Date_bin.unique()[0]].astype({'Latitude':float,'Longitude':float}),output_path='',**kwargs)
# ggplot_raster(dataframe=df_psd.dropna(subset=[var])[df_psd.Date_bin==df_psd.Date_bin.unique()[0]].astype({'Latitude':float,'Longitude':float}),output_path='{}/figures/Map_{}.svg'.format(str(path_to_git),var),**kwargs)
nbss_summary=nbss_summary.assign(Month=pd.to_datetime(nbss_summary.Date_bin,format='%Y-%m').dt.strftime('%B'),month=pd.to_datetime(nbss_summary.Date_bin,format='%Y-%m').dt.strftime('%m').astype(float))
nbss_summary=pd.merge(nbss_summary,df_psd.astype({'Latitude':float,'Longitude':float})[['Latitude','Longitude','Date_bin','PSD_slope','Picophytoplankton','Nanophytoplankton','Microphytoplankton']],how='left',on=['Latitude','Longitude','Date_bin'])

## Read AVISO summary file and merge to nbss_summary
if len(list(Path(path_to_git / 'data' / 'Environment' /'AVISO').expanduser().glob('*.parquet')))==0:
    df_aviso=pd.concat(map (lambda path:(df:=pd.read_parquet(path),df.rename(columns={'lat':'Latitude','lon':'Longitude'}).astype({'Latitude':str,'Longitude':str}))[-1],list(Path(path_to_git / 'data' / 'Environment'/'AVISO'/'FSLE').expanduser().rglob('*.parquet'))),axis=0).reset_index(drop=True).rename(columns={'fsle_max':'FSLE'})
    df_aviso.to_parquet(path_to_git / 'data' / 'Environment' / 'AVISO' / 'AVISO_all.parquet',index=False)
    [file.unlink(missing_ok=True) for file in list(Path(path_to_git / 'data' / 'Environment' /'AVISO').expanduser().rglob('*.parquet'))[1:]]  # Delete original files to save storage space

else:
    df_aviso=pd.read_parquet(path_to_git / 'data' / 'Environment' / 'AVISO' / 'AVISO_all.parquet')
df_aviso.columns
# var,kwargs = 'FSLE',{'mapping': 'y:Latitude,x:Longitude,fill:FSLE','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(list(reversed(sns.color_palette("RdGy",15).as_hex()))[::-1][0:7]+list(reversed(sns.color_palette("RdGy",15).as_hex()))[::-1][10:8:-1]): ''}
# ggplot_raster(dataframe=df_aviso.dropna(subset=[var])[df_aviso.Date_bin==df_aviso.Date_bin.unique()[1]].astype({'Latitude':float,'Longitude':float}),output_path='',**kwargs)
# ggplot_raster(dataframe=df_aviso.dropna(subset=[var])[df_aviso.Date_bin==df_aviso.Date_bin.unique()[1]].astype({'Latitude':float,'Longitude':float}),output_path='{}/figures/PSSdb_paper/Map_{}.svg'.format(str(path_to_git),var),**kwargs)

nbss_summary=pd.merge(nbss_summary,df_aviso.astype({'Latitude':float,'Longitude':float})[['Latitude','Longitude','Date_bin','FSLE','theta_max']].drop_duplicates(),how='left',on=['Latitude','Longitude','Date_bin'])
data=nbss_summary.copy()

# Generate scatter-plots of environmental parameter vs NBSS parameters
data_melt=data[(data.Intercept!=0) | (data.Intercept.isna()==False)].melt(id_vars=['PFT','PSD_slope','Picophytoplankton', 'Nanophytoplankton', 'Microphytoplankton','CHL','Temperature','PAR','AOU','Nitrates','Phosphates','Silicates','Iron','Aerosol_thickness','DUWT001','SLA','FSLE','MLD'],value_vars=['Slope','Intercept','Min_size','Max_size'])
data_melt['variable']=pd.Categorical(data_melt.variable,categories=['Intercept','Slope','Min_size','Max_size'])
data_melt['SLA']=data_melt.SLA-np.nanmean(data_melt.SLA)
for taxa in data_melt.PFT.unique():
    for var in ['PSD_slope', 'Picophytoplankton', 'Nanophytoplankton', 'Microphytoplankton', 'CHL', 'Temperature', 'PAR', 'AOU', 'Nitrates', 'Phosphates', 'Silicates', 'Iron', 'Aerosol_thickness', 'DUWT001', 'SLA', 'FSLE', 'MLD']:
        var = var + '/1000' if var.split('/')[0] in '\t'.join(['Nitrates', 'Phosphates','Silicates']) else var  # Nutrients are expressed in cubic meters to integrate over depth
        axis_trans = 'log10' if var.split('/')[0] in '\t'.join(['CHL', 'Nitrates', 'Phosphates', 'Silicates', 'Iron','DUWT001','AOU']) else 'identity'# 'sqrt' if var in '\t'.join(['Picophytoplankton', 'Nanophytoplankton', 'Microphytoplankton']) else 'identity'
        plot = (ggplot(data_melt.query('PFT=="{}"'.format(taxa))) +
                facet_grid('variable~PFT', scales='free', space='free') +
                geom_point(mapping=aes(x='{}'.format(var), y='value', size='CHL'), fill='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 0),color='black') +
                scale_size_area(max_size=5) +  scale_x_continuous(trans=axis_trans)+ #
                labs(y=r'', x=' '.join(var.split('_')), size='') +
                theme_paper ).draw(show=True)
        plot.set_size_inches(1.2 * 2 * 1, 6)  # len(data_melt['group'].unique())
        plot.savefig(fname='{}/figures/PSSdb_paper/{}_{}_vs_parameters_{}.svg'.format(str(path_to_git),re.compile('[^a-zA-Z]').sub('', var), taxa), dpi=300, bbox_inches='tight')

## Correlations
env_group = ['Instrument', 'PFT', 'Total_abundance', 'Min_size', 'Max_size', 'Intercept', 'Slope', 'PSD_slope', 'Picophytoplankton', 'Nanophytoplankton', 'Microphytoplankton', 'Temperature', 'PAR', 'AOU', 'Nitrates','Phosphates', 'Silicates', 'Iron', 'Aerosol_thickness', 'DUWT001', 'SLA', 'FSLE', 'MLD', 'CHL']
env_group = [column for column in env_group if column in nbss_summary.columns]
data=data.dropna(subset=env_group)
data_pval = data[env_group].groupby(env_group[0:2]).corr(method=lambda x, y: spearman_correlation(x, y)) - np.tile(np.eye(data[env_group].groupby(env_group[0:2]).corr(method=lambda x, y: spearman_correlation(x, y)).shape[1], data[env_group].groupby(env_group[0:2]).corr(method=lambda x, y: spearman_correlation(x, y)).shape[1]), int(data[env_group].groupby(env_group[0:2]).corr(method=lambda x, y: spearman_correlation(x, y)).shape[0] / data[env_group].groupby(env_group[0:2]).corr(method=lambda x, y: spearman_correlation(x, y)).shape[1])).T
data_pval = data_pval.where( np.tile(np.triu(np.ones(data_pval.iloc[0:data_pval.shape[1]].shape)).astype(bool), int(data_pval.shape[0] / data_pval.shape[1])).T)
p = data_pval.applymap(lambda x: ''.join(['*' for t in [.05, .01, .001] if x <= t]))
data_summary = pd.merge(data[env_group].groupby(env_group[0:2]).corr(method='spearman').where(np.tile(np.triu(np.ones(data_pval.iloc[0:data_pval.shape[1]].shape)).astype(bool), int(data_pval.shape[0] / data_pval.shape[1])).T).melt(ignore_index=False).reset_index().set_axis(labels=env_group[0:2] + ['var1', 'var2', 'value'], axis=1), p.melt(ignore_index=False).reset_index().set_axis(labels=env_group[0:2] + ['var1', 'var2', 'pvalue'], axis=1), how='left', on=env_group[0:2] + ['var1', 'var2'])
data_summary = data_summary.assign(lab_text=lambda x: np.round(x.value, 2).astype(str) + x['pvalue'])
data_summary['var1'] = pd.Categorical(data_summary.var1, ordered=True, categories=['Total_abundance', 'Min_size', 'Max_size', 'Intercept', 'Slope'][ ::-1] + ['PSD_slope', 'Picophytoplankton', 'Nanophytoplankton','Microphytoplankton', 'Temperature', 'CHL', 'MLD', 'PAR','AOU', 'Nitrates', 'Phosphates', 'Silicates', 'Iron', 'Aerosol_thickness', 'DUWT001', 'SLA', 'FSLE'][::-1])
data_summary['var2'] = pd.Categorical(data_summary.var2, ordered=True,categories=['Total_abundance', 'Min_size', 'Max_size', 'Intercept', 'Slope'] + ['PSD_slope','Picophytoplankton', 'Nanophytoplankton', 'Microphytoplankton','Temperature','CHL', 'MLD','PAR','AOU','Nitrates','Phosphates','Silicates','Iron','Aerosol_thickness','DUWT001', 'SLA','FSLE'][ ::-1])
data_summary = data_summary[(data_summary.var2.isin(['Total_abundance', 'Min_size', 'Max_size', 'Intercept', 'Slope'])) & (data_summary.var1.isin( ['PSD_slope', 'Picophytoplankton', 'Nanophytoplankton', 'Microphytoplankton', 'Temperature', 'CHL', 'MLD','PAR', 'AOU', 'Nitrates', 'Phosphates', 'Silicates', 'Iron', 'Aerosol_thickness', 'DUWT001', 'SLA','FSLE']))].sort_values(['PFT', 'var1']).reset_index(drop=True)
data_summary['p_value'] = data_summary.groupby(['PFT', 'var1']).apply(lambda x: pd.DataFrame({'p_value': np.where((x['pvalue'].str.len() == max(x['pvalue'].str.len())) & (np.abs(x.value) == max(np.abs(x.value))), x.pvalue, '')},index=x.index))
data_summary['p_value'] = data_summary.groupby(['PFT', 'var1']).apply(lambda x: pd.DataFrame({'p_value': np.where((x['pvalue'].str.len() == max(x['pvalue'].str.len())) & (np.abs(x.value) == max(np.abs(x.value))), x.pvalue, '')[ (x['pvalue'].str.len() == max(x['pvalue'].str.len())) & (np.abs(x.value) == max(np.abs(x.value)))][0]}, index=x.index))

plot = (ggplot(data_summary.query('var2.isin(["Intercept","Slope"])')) +
        facet_wrap('~PFT', ncol=1) +
        geom_col(mapping=aes(x='var1', y='value', group='var2', fill='var2'), stat='identity', position='dodge',show_legend=False, color='black', size=0.01,alpha=.6) +
        geom_text(mapping=aes(x='var1', y='value', group='var2', label='pvalue', nudge_y='-0.1+1*(value>0)'), color='black', family='Times New Roman', size=8, angle=90, position=position_dodge(width=1)) +
        labs(x=r'', y=r'Spearman correlation coefficient', fill='') +
        scale_fill_manual(values={'Intercept':'#{:02x}{:02x}{:02x}'.format(226,219,227),'Slope':'#{:02x}{:02x}{:02x}'.format(183,190,200)}) +
        theme_paper+theme(panel_grid_minor_x=element_line(color='black'),axis_text_x=element_text(angle=-45,ha='left'))).draw(show=True)
plot.set_size_inches(5.5, 5.6)
plot.savefig(fname='{}/figures/PSSdb_paper/Environmental_correlation_NBSS.svg'.format(str(path_to_git)),dpi=300, bbox_inches='tight')

# Save NBSS and environmental parameters for model input
path_to_input=path_to_git / 'data' / 'Model'
path_to_input.mkdir(exist_ok=True,parents=True)
data.to_csv('{}/Model_input.csv'.format(str(path_to_input)),index=False)

# Generate full environmental dataset for global predictions
df_environement = pd.merge(df_woa.astype({'Latitude': str, 'Longitude': str}).loc[:, ['Latitude', 'Longitude', 'Month','Oxygen', 'AOU', 'Salinity', 'MLD', 'Phosphates', 'Silicates', 'Nitrates','Temperature']], ( df_copernicus.astype({'Latitude': str, 'Longitude': str}).assign(Month=pd.to_datetime(df_copernicus.Date_bin,format='%Y-%m').dt.strftime('%B')).loc[ :, ['Latitude', 'Longitude','Date_bin','Month', 'Iron', 'CHL','SLA']]), how='right', on=['Latitude', 'Longitude', 'Month']) #df_copernicus.Date_bin == '2008-01'
df_environement = pd.merge(df_environement, (df_nasa.astype({'Latitude': str, 'Longitude': str}).loc[:, ['Latitude', 'Longitude','Date_bin', 'PAR', 'Aerosol_thickness', 'DUWT001']]), how='left', on=['Latitude', 'Longitude', 'Date_bin'])
df_environement = pd.merge(df_environement, (df_psd.astype({'Latitude': str, 'Longitude': str}).loc[:, ['Latitude', 'Longitude','Date_bin','PSD_slope', 'Picophytoplankton', 'Nanophytoplankton', 'Microphytoplankton']]), how='left', on=['Latitude', 'Longitude', 'Date_bin'])
df_environement = pd.merge(df_environement, (df_aviso.astype({'Latitude': str, 'Longitude': str}).loc[:, ['Latitude', 'Longitude','Date_bin', 'FSLE']]).drop_duplicates(), how='left', on=['Latitude', 'Longitude', 'Date_bin'])
df_env_all = df_environement.dropna().reset_index(drop=True)
df_environement['Group_index']=df_environement.reset_index(drop=True).index
gdf = gpd.GeoDataFrame(df_environement[['Group_index','Longitude', 'Latitude']].drop_duplicates().dropna(), geometry=gpd.points_from_xy( df_environement[['Group_index', 'Longitude', 'Latitude']].drop_duplicates().dropna().Longitude,df_environement[['Group_index','Longitude', 'Latitude']].drop_duplicates().dropna().Latitude))
df_environement['Study_area'] = pd.merge(df_environement, gpd.tools.sjoin(gdf, oceans, predicate="within", how='left')[['Group_index', 'name']], how='left',on='Group_index')['name'].astype(str)
df_environement['Longhurst_area'] = pd.merge(df_environement, gpd.tools.sjoin(gdf, longhurst, predicate="within", how='left')[['Group_index', 'ProvDescr']], how='left',on='Group_index')['ProvDescr'].astype(str)
df_environement['Longhurst_subarea'] = df_environement.Longhurst_area.replace(longhurst_dict)
df_environement['Hemisphere'] = np.where(df_environement.Latitude.astype(float)>0,'North','South')

# Generate a netcdf xarray from Pandas DataFrame
array=df_env_all.rename(columns={'Date_bin':'Datetime'}).drop(columns=['Month','Oxygen','Salinity']).astype({'Latitude':float,'Longitude':float}).set_index(['Latitude', 'Longitude','Datetime']).to_xarray()
array['Latitude'].attrs={'units':'decimal_degrees','description':'Midpoint of the 1 degree latitudinal grid'}
array['Longitude'].attrs={'units':'decimal_degrees','description':'Midpoint of the 1 degree longitudinal grid'}
array['Datetime'].attrs={'units':'%Y-%m','description':'Monthly average corresponding to the final resolution of PSSdb products'}
array['AOU'].attrs={'units':'micromole per cubic meter','long_name':'Apparent Oxygen Utilization','description':'Surface concentration of apparent oxygen utilization from World Ocean Atlas'}
array['MLD'].attrs={'units':'meters','long_name':'Mixed Layer Depth','description':'Mixed Layer Depth from World Ocean Atlas'}
array['Phosphates'].attrs={'units':'micromole per cubic meter','description':'Surface concentration of phosphate from World Ocean Atlas'}
array['Nitrates'].attrs={'units':'micromole per cubic meter','description':'Surface concentration of nitrate from World Ocean Atlas'}
array['Silicates'].attrs={'units':'micromole per cubic meter','description':'Surface concentration of silicate from World Ocean Atlas'}
array['Temperature'].attrs={'units':'degrees_celsius','description':'Sea surface temperature from World Ocean Atlas'}
array['Iron'].attrs={'units':'micromole per cubic meter','description':'Surface concentration of dissolved iron coumpounds from the biogeochemical model PISCES'}
array['CHL'].attrs={'units':'milligram per cubic meters','description':'Surface concentration of chlorophyll-a measured from satellite'}
array['SLA'].attrs={'units':'meters','long_name':'Sea Level Anomaly','description':'Sea Level Anomaly measured from satellite'}
array['PAR'].attrs={'units':'mole einstein per square meter per day','long_name':'Photosynthetically Available Radiance','description':'Surface light intensity measured from satellite'}
array['Aerosol_thickness'].attrs={'units':'unitless','long_name':'Aerosol optical thickness','description':'Proxy for surface emission of Fe-rich aerosols measured from satellite'}
array['DUWT001'].attrs={'units':'kilogram per square meter per second','long_name':'Wet Dust Deposition Bin 001','description':'Proxy for surface deposition of wet dust measured from satellite'}
array['PSD_slope'].attrs={'units':'particles per liter per square micrometers','long_name':'Particle Size Distribution Slope','description':'Slope of the particle size distribution measured from satellite'}
array['Picophytoplankton'].attrs={'units':'unitless','long_name':'Fraction of picophytoplankton','description':'Fraction of picophytoplankton measured from satellite'}
array['Nanophytoplankton'].attrs={'units':'unitless','long_name':'Fraction of nanophytoplankton','description':'Fraction of nanophytoplankton measured from satellite'}
array['Microphytoplankton'].attrs={'units':'unitless','long_name':'Fraction of microphytoplankton','description':'Fraction of microphytoplankton measured from satellite'}
array.attrs={'creation':'Dataset generated on Sept 2023 by M. Dugenne', 'title':'Full environmental dataset used to predict the spectral biogeography of mesoplankton in the Global Ocean', 'summary':'This dataset represents the set of explanatory variables used to predict mesoplankton spectral biogeography as part of the Pelagic Sie Structure database project (https://pssdb.net)'} # add global attribute metadata
print(array)
array.to_netcdf('{}/Model_environment.nc'.format(str(path_to_input)))

# Plot climatologies
df_env_all = xr.open_dataset('{}/Model_environment.nc'.format(str(path_to_input))).to_dataframe().reset_index().dropna()#df_environement.dropna().reset_index(drop=True)
df_env_all=pd.merge(df_env_all,df_env_all.drop_duplicates(subset=['Longitude', 'Latitude'], ignore_index=True)[['Longitude', 'Latitude']].reset_index().rename( {'index': 'Group_index'}, axis='columns'),how='left',on=['Longitude', 'Latitude'])
gdf = gpd.GeoDataFrame(df_env_all[['Group_index','Longitude', 'Latitude']].drop_duplicates().dropna(), geometry=gpd.points_from_xy( df_env_all[['Group_index', 'Longitude', 'Latitude']].drop_duplicates().dropna().Longitude,df_env_all[['Group_index','Longitude', 'Latitude']].drop_duplicates().dropna().Latitude))
df_env_all['Study_area'] = pd.merge(df_env_all, gpd.tools.sjoin(gdf, oceans, predicate="within", how='left')[['Group_index', 'name']], how='left',on='Group_index')['name'].astype(str)
df_env_all['Longhurst_area'] = pd.merge(df_env_all, gpd.tools.sjoin(gdf, longhurst, predicate="within", how='left')[['Group_index', 'ProvDescr']], how='left',on='Group_index')['ProvDescr'].astype(str)
df_env_all['Longhurst_subarea'] = df_env_all.Longhurst_area.replace(longhurst_dict)
df_env_all['Hemisphere'] = np.where(df_env_all.Latitude.astype(float)>0,'North','South')
df_env_all['Region']=df_env_all.Study_area+'\n'+df_env_all.Longhurst_subarea
df_env_all['month']=df_env_all.Datetime.str[5:]
for var in [ 'Iron', 'Aerosol_thickness', 'DUWT001']:
    var = var + '/1000' if var.split('/')[0] in '\t'.join(['Nitrates', 'Phosphates', 'Silicates']) else var  # Nutrients are expressed in cubic meters to integrate over depth
    axis_trans = 'log10' if var.split('/')[0] in '\t'.join(['CHL', 'Nitrates', 'Phosphates', 'Silicates', 'Iron','DUWT001']) else 'identity'  # 'sqrt' if var in '\t'.join(['Picophytoplankton', 'Nanophytoplankton', 'Microphytoplankton']) else 'identity'
    plot = (ggplot(df_env_all[(df_env_all.Longhurst_subarea.isin(['Equatorial-Upwellings']))]) +
            facet_grid('Hemisphere~Region', scales='free_y') +
            geom_boxplot(mapping=aes(x='month', y=var,group='month'),alpha=0.2,fill='gray',color='black',outlier_size=0.0) +
            scale_y_continuous(trans=axis_trans) +
            labs(x=r'', y=r'') +
            theme_paper+theme(strip_text=element_text(size=4))).draw(show=True)
    plot.set_size_inches(5, 3)
    plot.savefig(fname='{}/figures/PSSdb_paper/{}_climatologies_Longhurst.svg'.format(str(path_to_git), re.compile('[^a-zA-Z]').sub('', var)),dpi=300, bbox_inches='tight')