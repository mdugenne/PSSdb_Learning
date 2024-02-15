## Objective: This script contains the set of functions required to analyze and predict PSSdb NBSS as a response of environmental factors

## Modules:
import warnings

import numpy as np
import pandas as pd
import statsmodels

warnings.filterwarnings(action='ignore')
# Modules for webpage handling/scraping:
from pathlib import Path
import urllib3
import requests
from bs4 import BeautifulSoup # Use pip install bs4 and pip install lxml
from funcy import join_with
from requests_html import HTMLSession
import urllib
import copernicusmarine as cm # Use pip install copernicusmarine

# Progress bar modules
#import progressbar # Attention: Use pip install progressbar2 to install
from tqdm import tqdm
from datetime import date

# Module for netcdf:
import xarray as xr # pip install netCDF4
from dateutil.relativedelta import relativedelta
import itertools
from natsort import natsorted
import scipy
import seaborn as sns
palette_oxygen=list(reversed(sns.color_palette("inferno",15).as_hex()))
palette_temp=list((sns.color_palette("BuPu",15).as_hex()))#colorspace.diverging_hcl(name="Purple-Blue").colors()
palette_chl=list(reversed(sns.color_palette("GnBu",15).as_hex()))#colorspace.diverging_hcl(name="GnBu").colors()
palette_taxa={'puff':'#{:02x}{:02x}{:02x}'.format(145,111,124),'tuff':'#{:02x}{:02x}{:02x}'.format(66,219,189),'Acantharia(Class)':'#{:02x}{:02x}{:02x}'.format(211,141,95),'Polycystina(Class)':'#{:02x}{:02x}{:02x}'.format(211,20,95),'Foraminifera(Phylum)(Class)':'#{:02x}{:02x}{:02x}'.format(95,111,188),'Globothalamea(Class)':'#{:02x}{:02x}{:02x}'.format(95,111,188),'Thecofilosea(Class)':'#{:02x}{:02x}{:02x}'.format(244,238,188),'Malacostraca(Class)':'#{:02x}{:02x}{:02x}'.format(208,7,46),'Ostracoda(Class)':'#{:02x}{:02x}{:02x}'.format(83,108,103),'Copepoda(Class)':'#{:02x}{:02x}{:02x}'.format(212,255,42),'Carnivore':'#{:02x}{:02x}{:02x}'.format(145,111,124),'Omnivore':'#{:02x}{:02x}{:02x}'.format(170,170,255)}

# Copernicus API
import motuclient# In terminal, Use python -m pip install motuclient==1.8.4 --no-cache-dir
import os
import cdsapi # In terminal : touch ~/.cdsapirc
# nano ~/.cdsapirc
# url: https://cds.climate.copernicus.eu/api/v2
#  key: {uid}:{api-key}
# pip install cdsapi

import yaml
path_to_git=Path('~/GIT/PSSdb_model').expanduser()
path_to_config = path_to_git / 'scripts' / 'Configuration_masterfile.yaml'
with open(path_to_config, 'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

from itertools import compress
from time import strptime
from plotnine import *  # Python equivalent of ggplot2. Use shell command: pip install plotnine. Do not execute in Pycharm (bug- no fix yet): https://youtrack.jetbrains.com/issue/PY-42390
import plotnine
theme_paper=theme(panel_grid=element_blank(), legend_position='right',
              panel_background=element_rect(fill='#ffffff00'), legend_background=element_rect(fill='#ffffff00'),
              strip_background=element_rect(fill='#ffffff00'),
              panel_border=element_rect(color='#222222'),
              legend_title=element_text(family='Times New Roman', size=8),
              legend_text=element_text(family='Times New Roman', size=8, rotation=90),
              axis_title=element_text(family='Times New Roman', size=8, linespacing=1),
              axis_text_x=element_text(family='Times New Roman', size=8, linespacing=0),
              axis_text_y=element_text(family='Times New Roman', size=8, rotation=90, linespacing=1),
              plot_background=element_rect(fill='#ffffff00'), figure_size=(5, 3))

import geopandas as gpd
from shapely.geometry import Polygon, mapping
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}), pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}), pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
oceans = gpd.read_file(list(Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent.rglob('goas_v01.shp'))[0])
longhurst=gpd.read_file(list(Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent.rglob('Longhurst*.shp'))[0])
longhurst_dict={'Trades - Eastern Tropical Atlantic Province':'Subtropical',
       'Trades - Pacific Equatorial Divergence Province':'Equatorial-Upwellings',
       'Trades - W. Pacific Warm Pool Province':'Subtropical',
       'Trades - Western Tropical Atlantic Province':'Subtropical',
       'Coastal - Guianas Coastal Province':'Equatorial-Upwellings',
       'Coastal - Chile-Peru Current Coastal Province':'Equatorial-Upwellings',
       'Coastal - Sunda-Arafura Shelves Province':'Coastal',
       'Coastal - NW Arabian Upwelling Province':'Equatorial-Upwellings',
       'Trades - Indian Monsoon Gyres Province':'Subtropical',
       'Coastal - Guinea Current Coastal Province':'Equatorial-Upwellings',
       'Coastal - Australia-Indonesia Coastal Province':'Coastal',
       'Trades - South Atlantic Gyral Province (SATG)':'Subtropical',
       'Westerlies - S. Pacific Subtropical Gyre Province':'Subtropical',
       'Coastal - Brazil Current Coastal Province':'Coastal',
       'Trades - Archipelagic Deep Basins Province':'Coastal',
       'Coastal - E. Africa Coastal Province':'Coastal',
       'Trades - Indian S. Subtropical Gyre Province':'Subtropical',
       'Coastal - East Australian Coastal Province':'Coastal',
       'Coastal - Benguela Current Coastal Province':'Equatorial-Upwellings',
       'Westerlies - S. Subtropical Convergence Province':'Temperate-Subpolar',
       'Westerlies - Tasman Sea Province':'Temperate-Subpolar',
       'Coastal - New Zealand Coastal Province':'Coastal',
       'Coastal - SW Atlantic Shelves Province':'Coastal',
       'Westerlies - Subantarctic Province':'Temperate-Subpolar',
       'Polar - Antarctic Province':'Polar',
       'Polar - Austral Polar Province':'Polar',
       'Coastal - Central American Coastal Province':'Coastal',
       'Trades - N. Pacific Equatorial Countercurrent Province':'Equatorial-Upwellings',
       'Trades - N. Atlantic Tropical Gyral Province (TRPG)':'Subtropical',
       'Trades - Caribbean Province':'Coastal',
       'Coastal - W. India Coastal Province':'Coastal',
       'Coastal - E. India Coastal Province':'Coastal',
       'Trades - N. Pacific Tropical Gyre Province':'Subtropical',
       'Westerlies - N. Pacific Subtropical Gyre Province (West)':'Subtropical',
       'Westerlies - Kuroshio Current Province':'Kuroshio',
       'Coastal - Canary Coastal Province (EACB)':'Equatorial-Upwellings',
       'Coastal - Red Sea, Persian Gulf Province':'Coastal',
       'Coastal - California Upwelling Coastal Province':'Equatorial-Upwellings',
       'Coastal - China Sea Coastal Province':'Coastal',
       'Westerlies - N. Atlantic Subtropical Gyral Province (East) (STGE)':'Subtropical',
       'Westerlies - N. Atlantic Subtropical Gyral Province (West) (STGW)':'Subtropical',
       'Coastal - NW Atlantic Shelves Province':'Temperate-Subpolar',
       'Westerlies - Mediterranean Sea, Black Sea Province':'Mediterranean',
       'Westerlies - Gulf Stream Province':'Temperate-Subpolar',
       'Westerlies - N. Pacific Polar Front Province':'Temperate-Subpolar',
       'Polar - N. Pacific Epicontinental Province':'Temperate-Subpolar',
       'Westerlies - Pacific Subarctic Gyres Province (West)':'Temperate-Subpolar',
       'Westerlies - N. Atlantic Drift Province (WWDR)':'Temperate-Subpolar',
       'Coastal - NE Atlantic Shelves Province':'Temperate-Subpolar',
       'Westerlies - Pacific Subarctic Gyres Province (East)':'Temperate-Subpolar',
       'Polar - Atlantic Arctic Province':'Polar',
       'Coastal - Alaska Downwelling Coastal Province':'Coastal',
       'Polar - Boreal Polar Province (POLR)':'Polar',
       'Polar - Atlantic Subarctic Province':'Polar'}

from shapely.geometry import Polygon, mapping
from shapely import wkt

#pip install xgboost
import xgboost as xgb
from xgboost import XGBClassifier
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold,GridSearchCV
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import StandardScaler
from scipy.stats import t
from scipy.stats import spearmanr

import shutil
import graphviz #use conda install graphviz
from numpy import arange
from pandas import read_csv
from sklearn.linear_model import LassoCV,RidgeCV
from matplotlib import pyplot as plt
from funcy import project
import statsmodels.api as sm
import re
from scipy import stats
from statsmodels.stats.outliers_influence import summary_table

from scipy.stats import poisson,norm,lognorm # Estimate uncertainties assuming count detection follows a Poisson distribution and size a normal distribution
from itertools import compress

import warnings
warnings.filterwarnings('ignore', module='urllib3')


## Functions start here:

# Perform linear regression on NBSS datasets to calculate intercepts at the lowest observed size class rather than log_10(1)=0
def NBSS_regression(nbss,selection='indexing',threshold_sd=0.5,threshold_count=0.2,threshold_size=0.2,n_bins=1):
    """
    Objective: This function computes the slope and intercept of the  Normalized Biovolume Size Spectrum (in \u03BCm\u00b3 dm\u207B\u00b3 \u03BCm\u207B\u00b3) from a log-transformed linear regression
        :param nbss: Dataframe containing size class (size_class_mid, in micrometers) and associated Normalized Biovolume Size Spectrum (NBSS, in \u03BCm\u00b3 dm\u207B\u00b3 \u03BCm\u207B\u00b3) columns.
        :param selection: Method to select the linear portion of the NBSS. Options are 'indexing' (i.e masking size smaller than the NBSS peak and larger than classes presenting 1 ROI) or 'std' (i.e. Coefficient of variation superior to the threshold are ignored,NBSS_std column required)
        :param threshold_sd: Coefficient of variation (std/mean) used to determine data selection (True if lower than threshold)
        :param threshold_count: Uncertainty threshold used to discard size classes whose number of observations is insufficient according to the Poisson distribution
        :param threshold_size: Uncertainty threshold used to discard size classes whose number of pixels is close to the camera resolution according to the Gaussian distribution
        :param n_bins: Number of consecutive bins that should be empty in order to discard the next upper size classes of the size spectrum with partimony

    :return: Returns the NBSS linear regression fit parameters
    """
    native_index=nbss.index
    nbss=nbss.astype(dict(zip(['size_class_mid','NBSS'],[float]*2))).sort_values(['size_class_mid']).reset_index().copy()

    nbss['Count_uncertainty']=poisson.pmf(k=nbss['NBSS_count'], mu=nbss['NBSS_count']) if 'NBSS_count' in nbss.columns else 0
    nbss.loc[nbss.Count_uncertainty>threshold_count,'NBSS']=np.nan
    nbss['Size_uncertainty']=norm.pdf((1/6)*np.pi*(2*(nbss.size_class_pixel.astype(float)/np.pi)**0.5)**3,loc=(1/6)*np.pi*(2*(1/np.pi)**0.5)**3,scale=(1/6)*np.pi*(2*(1/np.pi)**0.5)**3) if 'size_class_pixel' in nbss.columns else 0
    nbss.loc[nbss.Size_uncertainty>threshold_size,'NBSS']=np.nan
    nbss['selected'] = nbss.NBSS.isna()==False

    if selection=='indexing': # Select observations between maximum and last size class before observing n-bins consecutive empty size classes
        nbss['selected'] = (nbss.NBSS.index >= nbss.NBSS.index.to_list()[np.argmax(nbss.NBSS)])
        if any(nbss.drop(nbss.NBSS.index.to_list()[0:np.argmax(nbss.NBSS)]).NBSS.isnull().astype(int).groupby(nbss.drop(nbss.NBSS.index.to_list()[0:np.argmax(nbss.NBSS)]).NBSS.notnull().astype(int).cumsum()).sum() >= n_bins + 1):
            nbss.loc[nbss.drop(nbss.NBSS.index.to_list()[0:np.argmax(nbss.NBSS)]).index[nbss.drop(nbss.NBSS.index.to_list()[0:np.argmax(nbss.NBSS)]).NBSS.isnull().astype(int).groupby(nbss.drop(nbss.NBSS.index.to_list()[0:np.argmax(nbss.NBSS)]).NBSS.notnull().astype(int).cumsum()).sum().index[np.argwhere(nbss.drop(nbss.NBSS.index.to_list()[0:np.argmax(nbss.NBSS)]).NBSS.isnull().astype(int).groupby(nbss.drop(nbss.NBSS.index.to_list()[0:np.argmax(nbss.NBSS)]).NBSS.notnull().astype(int).cumsum()).sum().to_numpy() >= n_bins + 1)[0][0]]]:,'selected'] = False


    if selection=='std': # Select observations based on bootstrapped uncertainties
        nbss['selected']=(nbss.NBSS.index >= nbss.NBSS.index.to_list()[np.argmax(nbss.NBSS)]) & ((nbss.NBSS_std/nbss.NBSS)<=threshold_sd) if (np.argmax(nbss.NBSS) != 0) & (len( (nbss.NBSS_std/nbss.NBSS)<=threshold_sd)> 0) else [True] * len(nbss.NBSS) if (np.argmax(nbss.NBSS) == 0) & (len((nbss.NBSS_std/nbss.NBSS<=threshold_sd)==0)) else [False] * len(nbss)

    boolean_selection = nbss.assign(native_index=pd.Categorical(nbss['index'].values,native_index,ordered=True)).sort_values(['native_index']).selected.to_numpy()
    nbss_subset=nbss[(nbss.selected==True) & (nbss.NBSS.isna()==False)]

    nbss_dict={}
    try :
        X = np.log10((1 / 6) * np.pi * ((1e-03 * nbss_subset.size_class_mid) ** 3))  # in cubic millimeters
        Y = np.log10(nbss_subset.NBSS) if not all(nbss_subset.NBSS == 0) else np.log10(nbss_subset.NBSS + 1)
        X = sm.add_constant(X)

        wls_model=sm.WLS(Y, X )
        reg = wls_model.fit()
        st, data, ss2 = summary_table(reg, alpha=0.05)
        predict_mean_ci_low, predict_mean_ci_upp = data[:, 4:6].T
        predict_ci_low, predict_ci_upp = data[:, 6:8].T
        # Get uncertainty at each original x input value
        ci_uncertainty = (predict_ci_upp - reg.predict()) / reg.predict()

        #reg.params,reg.pvalues

        nbss_dict= dict({'Selection':[boolean_selection.tolist()],'Min_size':nbss_subset.size_class_mid.values[np.argmax(nbss_subset.NBSS)],'Max_size':nbss_subset.size_class_mid.values[-1],'Average_size':np.nansum(nbss.size_class_mid*(nbss.NBSS/np.nansum(nbss.NBSS))),'Slope':reg.params['size_class_mid'],'Slope_sd':reg.bse['size_class_mid'],'Total_abundance':nbss_subset.NBSS.sum(),'Max_abundance':nbss_subset.NBSS.max(),'Intercept_ref':reg.params['const'],'Intercept_ref_sd':reg.bse['const'],'Intercept': reg.predict()[0],'Intercept_sd':ci_uncertainty[0],'R2':reg.rsquared}) if reg.rsquared!=1 else dict({'Selection':[boolean_selection.tolist()],'Min_size':np.nan,'Max_size':np.nan,'Average_size':np.nan,'Slope':np.nan,'Slope_sd':np.nan,'Total_abundance':np.nan,'Max_abundance':np.nan,'Intercept_ref':np.nan,'Intercept_ref_sd':np.nan,'Intercept':np.nan,'Intercept_sd':np.nan,'R2':np.nan})
    except:
       nbss_dict= dict({'Selection':[boolean_selection.tolist()],'Min_size':np.nan,'Max_size':np.nan,'Average_size':np.nan,'Slope':np.nan,'Slope_sd':np.nan,'Total_abundance':np.nan,'Max_abundance':np.nan,'Intercept_ref':np.nan,'Intercept_ref_sd':np.nan,'Intercept':np.nan,'Intercept_sd':np.nan,'R2':np.nan})
    return nbss_dict

# Plot raster
def ggplot_raster(dataframe,output_path,**kwargs):
    kwargs = kwargs
    mapping_args = [values.split(',') for key, values in kwargs.items() if key == 'mapping']
    if len(mapping_args) == 1:
        kwargs_mapping = dict((value.split(':')[0], value.split(':')[1]) for value in mapping_args[0])
    else:
        print('Required argument(s) missing:{}'.format(getattr(getattr(plotnine, 'geom_point'), 'REQUIRED_AES')))
    geom_args = [values.split(',') for key, values in kwargs.items() if key == 'aes']
    if len(geom_args) == 1:
        kwargs_geom = dict((value.split(':')[0], float(value.split(':')[1])) if value.split(':')[0] not in ['fill', 'color', 'colour','size'] else (value.split(':')[0], value.split(':')[1]) for value in geom_args[0])
    else:
        kwargs_geom = dict((key, value) for key, value in getattr(getattr(plotnine, 'geom_point'), 'DEFAULT_AES').items() if key not in kwargs_mapping.keys())

    kwargs_eval = [key for key, value in kwargs.items() if value == '']
    if (len(dataframe) > 0 ) :
        plot = (ggplot(dataframe) +
                getattr(plotnine.mapping, 'aes')(**kwargs_mapping) +
                getattr(plotnine, 'geom_raster')(**{'na_rm':True}) +
                #getattr(plotnine, 'geom_point')(**kwargs_geom) +
                geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'),fill='black',color='black') +
                labs( x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)') +
                guides(fill=guide_colourbar(barwidth=10,direction='horizontal',title_position='top'))+
                theme(legend_position='bottom',panel_grid=element_blank(), panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family='serif', size=10), legend_text=element_text(family='serif', size=10),
                      axis_title=element_text(family='serif',size=10), axis_text_x=element_text(family='serif',size=10),
                      axis_text_y=element_text(family='serif',size=10, rotation=90),
                      plot_background=element_rect(fill='white'))+
                coord_cartesian(expand=False))
        if len(kwargs_eval) > 0:
            for i in range(len(kwargs_eval)):
                plot = plot + eval('plotnine.' + kwargs_eval[i])

            # Saving plot to output_path
        if len(output_path) > 0:
            path = Path(output_path).expanduser()
            print('Saving plot to', path, sep=' ')
            plot.save(filename=path, dpi=300, verbose=False, bbox_inches='tight',width=5,height=3)
        else:
            print('Printing plot. Save manually using the bottom icon')
            plot.draw(show=True)


    else:
        return print('Required argument(s) missing')

# Download all environmental variables for observed year(s)
## Functions linked to web scraping
def get_all_forms(url):
    """Returns all form tags found on a web page's `url` """
    with  HTMLSession() as session:
        # GET request
        res = session.get(url)
        # for javascript driven website
        # res.html.render()
        soup = BeautifulSoup(res.html.html, "html.parser")

    return soup.find_all("form")

# Retrieve inputs required for the form
def get_form_details(form):
    """Returns the HTML details of a form,
    including action, method and list of form controls (inputs, etc)"""
    details = {}
    # get the form action (requested URL)
    action = form.attrs.get("action").lower()
    # get the form method (POST, GET, DELETE, etc)
    # if not specified, GET is the default in HTML
    method = form.attrs.get("method", "get").lower()
    # get all form inputs
    inputs = []
    for input_tag in form.find_all("input"):
        # get type of input form control
        input_type = input_tag.attrs.get("type", "text")
        # get name attribute
        input_name = input_tag.attrs.get("name")
        # get the default value of that input tag
        input_value =input_tag.attrs.get("value", "")
        # add everything to that list
        inputs.append({"type": input_type, "name": input_name, "value": input_value})
    details["action"] = action
    details["method"] = method
    details["inputs"] = inputs
    return details

## Particle Size Distribution: slopes and fractions of pico-, nano-, and microphytoplankton retreived from satellite backscattering coefficients
def PSD_download(year_bin,url='https://doi.pangaea.de/10.1594/PANGAEA.939863?format=html#download',output_path=path_to_git / 'data' / 'Environment'/ 'PSD'):
    """
    Objective: This function uses web scraping modules to download satellite PSD data from Kostodinov et al. 2023 (https://doi.pangaea.de/10.1594/PANGAEA.939863?format=html#download).\n
         :param year_bin: The year (str) of the data of interest.
         :param url: The url used to locate the files to download
         :param output_path: The path for downloaded files storage.
         :return: Returns PSD monthly file at output_path/year_bin/*.nc
    """
    with requests.Session() as session:
        # Authenticate with email/password
        rsp = session.get(url, headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
        soup = BeautifulSoup(rsp.content, "lxml")
        href = { item.text:  item.attrs['href'] for item in soup.find_all('a') if len(item.attrs['href'])}
        rsp.close()
        href={key:ref for key,ref in href.items() if ('-fv5.0.nc' in key.split('_')[-1]) & (key.split('_')[-1].replace('-fv5.0.nc','')[0:4]==year_bin) & (len(key.split('_')[-1].replace('-fv5.0.nc',''))==6)}
        for file,file_url in href.items():
            path_to_file = output_path / year_bin / 'PSD_{}.nc'.format(file.split('_')[-1].replace('-fv5.0.nc', ''))
            path_to_file.parent.mkdir(exist_ok=True, parents=True)
            if path_to_file.exists():
                print('\nExisting {} found. Skipping current file'.format(path_to_file))
            else:
                with session.get(file_url , stream=True, headers={'Connection': 'Keep-Alive','User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'}) as rsp:
                    with open(path_to_file, "wb") as f:
                        for a_chunk in rsp.iter_content(chunk_size=131072):  # Loop over content, i.e. eventual HTTP chunks
                                f.write(a_chunk)

                        #shutil.copyfileobj(rsp.raw, f, length=16 * 1024 * 1024)

            print('Data downloaded at {}'.format(path_to_file))
        if len(list(Path(output_path/year_bin).expanduser().glob("*.parquet"))) == 0:
            df = pd.concat(map(lambda path: nc_to_df_PSD(path, grid_res=1), list(Path(output_path/year_bin).expanduser().glob("*.nc"))))
            df.to_parquet(Path(output_path/year_bin).expanduser() / str('PSD_{}.parquet'.format(year_bin)), index=False)
            [file.unlink(missing_ok=True) for file in list(Path(output_path / year_bin).expanduser().glob("*.nc"))] # Delete original files to save storage space
        if not Path(output_path / 'README.xml').exists():  # Write README file
            rsp = session.get(url=url, headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
            soup = BeautifulSoup(rsp.content, "lxml")
            dict_soup={item.text.split(':')[0]:item.text.split(':')[1] for item in soup.find(id='dataset').find_all(attrs={'class': 'row'})}
            df=pd.DataFrame({'shortname':list(dict_soup.keys()),'description':list(dict_soup.values())})
            rsp.close()
            with open(Path(output_path / 'README.xml'), 'w') as f:
                f.write(df.reset_index(drop=True).to_xml(namespaces={"":url}, prefix="", index=False, parser='etree'))

def nc_to_df_PSD(path,grid_res=1):
    df=pd.DataFrame()
    if path.exists():
        ds = xr.open_dataset(path, decode_times=False)
        ds=ds.rename({'Grid_latitudes_pixel_centers': 'lat', 'Grid_longitudes_pixel_centers': 'lon'})
        df = ds.to_dataframe().reset_index()
        df=df.assign(lat_bins=pd.cut(df.lat,bins=np.arange(-90, 91, grid_res),labels=list(pd.IntervalIndex(pd.cut(np.arange(-90, 91, grid_res), np.arange(-90, 91, grid_res))).mid[ 1:])),lon_bins=pd.cut(df.lon,bins=np.arange(-180, 181, grid_res),labels=list(pd.IntervalIndex(pd.cut(np.arange(-180, 181, grid_res), np.arange(-180, 181, grid_res))).mid[ 1:]))).groupby(['lat_bins','lon_bins']).mean().reset_index().astype(dict(zip(['lat_bins','lon_bins'],2*[float])))
        df['Date_bin']=ds.attrs['Year']+'-'+ds.attrs['Month'] if 'Year' in ds.attrs.keys() else ds.attrs['Months']
        # kwargs = {'mapping': 'y:lat_bins,x:lon_bins,fill:PSD_slope','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(palette_temp): ''}
        # ggplot_raster(dataframe=df,output_path='',**kwargs)

    return df

## Satellite mesoscale proxies: Sea Level Anomaly (indicative of cyclonic or anticyclonic eddies), and Finite Size Lyapunov Exponent (indicate the divergence rate between water masses in close proximity, in other words frontal filaments)
def AVISO_download(year_bin,username=cfg_pw['AVISO_user'],password=cfg_pw['AVISO_pass'],output_path=path_to_git / 'data' / 'Environment'/ 'AVISO',access='ftp'):
    """
    Objective: This function uses the motu python module to download data from AVISO (https://oceandata.sci.gsfc.nasa.gov/,https://github.com/clstoulouse/motu-client-python#Usage).\n
         :param year_bin: The year (str) of the data of interest.
         :param output_path: The path for downloaded files storage.
         :param username: Your AVISO username
         :param password: Your AVISO password
         :param output_path: The local path where new downloads are exported
         :param access: Accessibility options ('ftp' or 'motu'). Default is ftp (faster)
         :return: Returns NASA monthly file at output_path/parameter/period/*.nc
    """
    if pd.to_datetime(year_bin + '-01-01') in pd.date_range('1994-01-01', '2022-08-04'):
        print('Downloading Finite Size Lyapunov Exponent. Please wait')
        output_file=Path(output_path).expanduser()/  'FSLE' / year_bin
        output_file.mkdir(parents=True, exist_ok=True)
        if access=='motu':
            for month in np.arange(1,13,1):
                    for week,week_range in {'week-'+str(index):list(compress(pd.date_range(day, periods=7).to_pydatetime().tolist(),pd.date_range(day, periods=7).month==month))[0::len(list(compress(pd.date_range(day, periods=7).to_pydatetime().tolist(),pd.date_range(day, periods=7).month==month)))-1] for index,day in enumerate(pd.date_range(pd.to_datetime("{}-{}-01".format(year_bin,str(month).zfill(2))), periods=31)[::7]) if day.month==month}.items():
                        path_to_file = output_file / "{}_{}_{}_{}.nc".format(year_bin, str(month).zfill(2),week, 'FSLE')
                        if not path_to_file.exists():
                            os.system('python -m motuclient --user {} --pwd {} -m https://motu.aviso.altimetry.fr/motu-web/Motu -s AvisoFSLE -d dataset-duacs-dt-global-allsat-madt-fsle -x -180 -X 180 -y -90 -Y 90 -t "{}" -T "{}" --outputWritten netcdf -v fsle_max -o {} -f {}'.format(username,password,week_range[0].date().isoformat(),week_range[1].date().isoformat(),path_to_file.parent,path_to_file.name))
        else:
            for month in np.arange(1, 13, 1):
                with urllib.request.urlopen('ftp://{}:{}@ftp-access.aviso.altimetry.fr/value-added/lyapunov/delayed-time/global/{}/dt_global_allsat_madt_fsle_{}{}01_20210921.nc'.format(username,password,year_bin,year_bin,str(month).zfill(2))) as remote_file, open(output_file/'dt_global_allsat_madt_fsle_{}{}01_20210921.nc'.format(year_bin,str(month).zfill(2)), 'wb') as local_file:
                    shutil.copyfileobj(remote_file, local_file)


        if not Path(output_file.parent / 'README.xml').exists():  # Write README file
            os.system( 'python -m motuclient --motu https://motu.aviso.altimetry.fr/motu-web/Motu --service-id AvisoFSLE --product-id dataset-duacs-dt-global-allsat-madt-fsle --describe-product -o console --user {} --pwd {} --out-dir {} --out-name {}'.format(  username, password, output_file.parent, 'README.xml'))

        if len(list(Path(output_path /  'FSLE' / year_bin).expanduser().glob("*.parquet")))==0:
               df=pd.concat(map(lambda path: nc_to_df_AVISO(path, grid_res=1),list(Path(output_path /  'FSLE' / year_bin).expanduser().glob("*.nc"))))
               df.to_parquet(output_path /  'FSLE' / year_bin / str('FLSE_{}.parquet'.format(year_bin)),index=False)
               [file.unlink(missing_ok=True) for file in list(Path(output_path /  'FSLE' / year_bin).expanduser().glob("*.nc"))]  # Delete original files to save storage space
    else:
        print('AVISO data not available for {}'.format(year_bin))

def nc_to_df_AVISO(path,grid_res=1):
    df=pd.DataFrame()
    if path.exists():
        ds = xr.open_dataset(path, decode_times=False)
        ds['lon'] = (ds['lon'] + 180) % 360 - 180
        ds = ds.groupby_bins('lon', np.arange(-180, 181, grid_res)).mean().groupby_bins('lat', np.arange(-90, 91, grid_res)).mean()
        ds['lat_bins'] = pd.IntervalIndex(pd.cut(np.arange(-91, 90, grid_res), np.arange(-91, 90, grid_res))).mid[ 1:]  # -((-(1/(grid_res/2))*ds['lat'])//1)/(1/(grid_res/2))
        ds['lon_bins'] = pd.IntervalIndex(pd.cut(np.arange(-180, 181, grid_res), np.arange(-180, 181, grid_res))).mid[ 1:]  # -((-(1/(grid_res/2))*ds['lon'])//1)/(1/(grid_res/2))
        ds = ds.rename({'lat_bins': 'lat', 'lon_bins': 'lon'})
        df = ds.to_dataframe().reset_index()
        df['datetime']=pd.to_datetime('1950-01-01 00:00:00') + relativedelta(days=df.loc[0,'time'])
        df['Date_bin']=df['datetime'].dt.strftime('%Y-%m')
        # kwargs = {'mapping': 'y:lat,x:lon,fill:fsle_max','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(palette_temp[::-1]): ''}
        # ggplot_raster(dataframe=df,output_path='',**kwargs)

    return df

## Satellite ocean color and atmospheric iron proxies: Photosynthetically Available Radiance (indicative of daily light intensity), Aeorsol optical thickness (indicative of all aerosol emissions, including natural and anthropogenic dust or volcanic and wildfire ashes)
def NASA_download(year_bin,output_path=path_to_git / 'data' / 'Environment' / 'NASA' ,parameters=['OC','AD'],parameters_dict={'NASA':'https://oceandata.sci.gsfc.nasa.gov/directdataaccess/Level-3%20Mapped/','MERRA':'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2_MONTHLY/M2TMNXADG.5.12.4/'},appkey=cfg_pw['NASA_appkey']):
    """
    Objective: This function uses python web scrapping module to download data from NASA (https://oceandata.sci.gsfc.nasa.gov/).\n
         :param year_bin: The year (str) of the data of interest.
         :param output_path: The path for downloaded files storage.
         :param parameters: A list of the variable(s) to download. Default is ['OC','AD'] for Ocean color and Aerosol diagnostic
         :param parameters_dict: A dictionary of the url used to download data. Default is level 3 mapped products
         :param appkey: The appkey generated for a specific profile on NASA at https://oceandata.sci.gsfc.nasa.gov/appkey/
         :return: Returns NASA monthly file at output_path/parameter/period/*.nc
    """
    if 'OC' in parameters:
        # NASA: Photosynthetically Available Radiance and Aerosol optical thickness
        with requests.Session() as session:
            # Authenticate with email/password
            rsp = session.get(parameters_dict.get('NASA') + 'SeaWiFS/', headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
            soup = BeautifulSoup(rsp.content, "lxml")
            href = {item.text.replace('/', ''): 'https://oceandata.sci.gsfc.nasa.gov' + item.attrs['href'].replace(' ', '%20') for item in soup.find_all("tbody")[0].find_all("a")}
            rsp.close()
            rsp = session.get(parameters_dict.get('NASA') + 'Aqua-MODIS/', headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
            soup = BeautifulSoup(rsp.content, "lxml")
            href.update({item.text.replace('/', ''): 'https://oceandata.sci.gsfc.nasa.gov' + item.attrs['href'].replace(' ', '%20') for item in soup.find_all("tbody")[0].find_all("a")})
            rsp.close()
            rsp = session.get(href[year_bin], headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
            soup = BeautifulSoup(rsp.content, "lxml")
            rsp.close()
            href_month = { item.text.split(' ')[1]: 'https://oceandata.sci.gsfc.nasa.gov' + item.attrs['href'].replace(' ', '%20') for item in soup.find_all("tbody")[0].find_all("a") if item.text.split(' ')[0] == '01'}
            for month, url in href_month.items():
                rsp = session.get(url, headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
                soup = BeautifulSoup(rsp.content, "lxml")
                rsp.close()
                href_files = {item.text.split('.')[-3]: item.attrs['href'].replace(' ', '%20') for item in soup.find_all("tbody")[0].find_all("a") if len(re.findall(r'MO.PAR.par.9km.nc|MO.RRS.aot_869.9km.nc', item.text)) == 1}
                for variable, file_url in href_files.items():
                    path_to_file = output_path / variable / year_bin / "{}_{}_{}_MO.nc".format(year_bin, str(strptime(month,'%b').tm_mon).zfill(2), variable)
                    path_to_file.parent.mkdir(parents=True, exist_ok=True)
                    #
                    rsp = session.get(file_url + '?appkey={}'.format(appkey), stream=True, headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
                    with open(path_to_file, "wb") as fd:
                        for a_chunk in rsp.iter_content(chunk_size=131072):  # Loop over content, i.e. eventual HTTP chunks
                            fd.write(a_chunk)
                    print('Data downloaded at {}'.format(path_to_file))
                    rsp.close()
            """ Link has been discontinued
            if not Path(output_path/ 'README.xml').exists():  # Write README file
                rsp = session.get(url='https://oceancolor.gsfc.nasa.gov/docs/format/l2oc_modis/', headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
                soup = BeautifulSoup(rsp.content, "lxml")
                df=pd.read_html(rsp.text)#
                #df_data=df[1+df[0].rename(columns=df[0].iloc[0]).drop(df[0].index[0]).query('dimensions=="geophysical_data"').dimensions.index.values[0]]
                df_data=pd.concat(map(lambda ind: df[ind].rename(columns=df[ind].iloc[1]).drop(df[ind].index[0:2]).dropna(subset=['shortname']) if 'shortname' in df[ind].rename(columns=df[ind].iloc[1]).drop(df[ind].index[0:1]).columns else pd.DataFrame(),np.arange(0,len(df))))
                rsp.close()
                with open(Path(output_path/ 'README.xml'), 'w') as f:
                    f.write(df_data.reset_index(drop=True).to_xml(namespaces={"": "https://oceancolor.gsfc.nasa.gov/docs/format/l2oc_modis/"},prefix="",index=False,attr_cols=['shortname', 'datatype', 'units', 'long_name','standard_name','reference'],parser='etree'))
            """
            for variable in href_files.keys():
                if len(list(Path(output_path / variable / year_bin).rglob('*.parquet')))==0:
                    df=pd.concat(map(lambda path:nc_to_df_NASA(path,grid_res=1),list(Path(output_path / variable / year_bin).rglob('*.nc')))).reset_index(drop=True).sort_values(['Date_bin','lat','lon'])
                    df.to_parquet(output_path / variable / year_bin / str('{}_{}.parquet'.format(variable,year_bin)), index=False)
                    [file.unlink(missing_ok=True) for file in list(Path(output_path / variable / year_bin).rglob('*.nc'))]  # Delete original files to save storage space

    if 'AD' in parameters:
        # MERRA: Aerosol diagnostic including dust wet deposition
        with requests.Session() as session:
            rsp = session.get(parameters_dict.get('MERRA'), headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
            soup = BeautifulSoup(rsp.content, "lxml")
            href = {item.text.replace('/', ''): parameters_dict.get('MERRA') + item.attrs['href'].replace(' ', '%20') for item in soup.find_all("a")}
            rsp.close()
            rsp = session.get(href[year_bin], headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
            soup = BeautifulSoup(rsp.content, "lxml")
            rsp.close()
            href_month = { item.text[-11:-5]: (href[year_bin]).rstrip('contents.html') + item.attrs['href'].replace(' ', '%20') for item in soup.find_all("a") if  ('.html' in item.attrs['href']) and ('.nc4' in item.text)} #{ item.text[-10:-4]: (href[year_bin]).rstrip('contents.html') + item.attrs['href'].replace(' ', '%20') for item in soup.find_all("a") if  ('.nc4' in item.attrs['href']) and ('.nc4.xml' not in item.attrs['href'])}

            for month, file_url in href_month.items():
                variable='aerosol_diagnostic'
                path_to_file = output_path / variable / year_bin / "{}_{}_{}_MO.nc".format(year_bin, month[-2:], variable)
                path_to_file.parent.mkdir(parents=True, exist_ok=True)
                os.system("wget --user={} --password={} -O {} --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --content-disposition {}".format(cfg_pw['NASA_user'],cfg_pw['NASA_pass'],path_to_file,file_url.replace('.html','.nc4?DUWT001'))) #

            if not Path(output_path/ 'README.xml').exists():  # Write README file
                rsp = session.get(url='https://goldsmr4.gesdisc.eosdis.nasa.gov/dods/M2TMNXADG.info', headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
                soup = BeautifulSoup(rsp.content, "lxml")
                df=(pd.read_html(rsp.text)[1]).loc[:,1:2]
                df.columns=['shortname','datatype']
                df=df.assign(units='',long_name='',standard_name='',reference='')#
                df_data=pd.concat([pd.read_xml(Path(output_path/ 'README.xml')),df],axis=0)
                rsp.close()
                with open(Path(output_path/ 'README.xml'), 'w') as f:
                    f.write(df_data.reset_index(drop=True).to_xml(namespaces={"": "https://oceancolor.gsfc.nasa.gov/docs/format/l2oc_modis/,https://goldsmr4.gesdisc.eosdis.nasa.gov/dods/M2TMNXADG.info"},prefix="",index=False,attr_cols=['shortname', 'datatype', 'units', 'long_name','standard_name','reference'],parser='etree'))

            for month in href_month.keys():
                if len(list(Path(output_path / variable / year_bin).rglob('*.parquet')))==0:
                    df=pd.concat(map(lambda path:nc_to_df_NASA(path,grid_res=1),list(Path(output_path / variable / year_bin).rglob('*.nc')))).reset_index(drop=True).sort_values(['Date_bin','lat','lon'])
                    df.to_parquet(output_path / variable / year_bin / str('{}_{}.parquet'.format(variable,year_bin)), index=False)
                    [file.unlink(missing_ok=True) for file in list(Path(output_path / variable / year_bin).rglob('*.nc'))]  # Delete original files to save storage space

def nc_to_df_NASA(path,grid_res=1):
    df=pd.DataFrame()
    if path.exists():
        ds =xr.open_dataset(path, decode_times=False)
        if 'DUWT001' in list(ds.variables.keys()):
            ds_unit = ds.variables.get('DUWT001').attrs.get('units')
            ds_time = ds.attrs['RangeBeginningDate'] +'T00:00:00.0Z'
            ds['lat']=np.arange(float(ds.attrs['SouthernmostLatitude']), float(ds.attrs['NorthernmostLatitude'])+float(ds.attrs['LatitudeResolution']), float(ds.attrs['LatitudeResolution']))
            ds['lon']=np.arange(float(ds.attrs['WesternmostLongitude']), float(ds.attrs['EasternmostLongitude'])+float(ds.attrs['LongitudeResolution']), float(ds.attrs['LongitudeResolution']))
            ds = ds['DUWT001'].groupby_bins('lon', np.arange(-180, 181, grid_res)).mean().groupby_bins('lat',np.arange(-90, 91,grid_res)).mean()

        else:
            ds_time=ds.time_coverage_start
            ds_unit=ds.variables.get('par').attrs.get('units') if Path(path).parts[-3]=='par' else ds.variables.get('aot_869').attrs.get('units') if 'aot_869' in list(ds.variables.keys()) else ds.variables.get('pic').attrs.get('units')
            ds=ds['par'].groupby_bins('lon', np.arange(-180,181,grid_res)).mean().groupby_bins('lat', np.arange(-90,91,grid_res)).mean() if Path(path).parts[-3]=='par' else ds['aot_869'].groupby_bins('lon', np.arange(-180,181,grid_res)).mean().groupby_bins('lat', np.arange(-90,91,grid_res)).mean() if 'aot_869' in list(ds.variables.keys()) else ds['pic'].groupby_bins('lon', np.arange(-180,181,grid_res)).mean().groupby_bins('lat', np.arange(-90,91,grid_res)).mean()#(ds['par'].coarsen(boundary='trim',lon=int(1/np.absolute(np.diff(ds['lon'].values.tolist()))[0]),lat=int(1/np.absolute(np.diff(ds['lat'].values.tolist())[0]))).mean()) if Path(path).parts[-3]=='par' else (ds['aot_869'].coarsen(boundary='trim',lat=int(1/np.absolute(np.diff(ds['lat'].values.tolist())[0]))).mean()).coarsen(boundary='trim',lon=int(1/np.absolute(np.diff(ds['lon'].values.tolist()))[0])).mean()
        ds['lat_bins']=pd.IntervalIndex( pd.cut(np.arange(-91,90,grid_res),np.arange(-91,90,grid_res))).mid[1:]#-((-(1/(grid_res/2))*ds['lat'])//1)/(1/(grid_res/2))
        ds['lon_bins'] =pd.IntervalIndex( pd.cut(np.arange(-180,181,grid_res),np.arange(-180,181,grid_res))).mid[1:]#-((-(1/(grid_res/2))*ds['lon'])//1)/(1/(grid_res/2))
        ds=ds.rename({'lat_bins':'lat','lon_bins':'lon'})
        df =ds.to_dataframe().assign(units=ds_unit,datetime=pd.to_datetime(ds_time,format='%Y-%m-%dT%H:%M:%S.%fZ')).reset_index()
        df['Date_bin']=df['datetime'].dt.strftime('%Y-%m')
        #kwargs = {'mapping': 'y:lat,x:lon,fill:par','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(palette_temp): ''}
        # kwargs = {'mapping': 'y:lat,x:lon,fill:DUWT001','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(palette_temp): ''}
        # kwargs = {'mapping': 'y:lat,x:lon,fill:pic','scale_fill_gradientn(trans="log10",colors={},guide=guide_colorbar(direction="vertical"))'.format(palette_temp): ''}
        #ggplot_raster(dataframe=df,output_path='',**kwargs)

    return df

## Satellite chlorophyll-a concentration and biogeochemical model: surface Chlorophyll (proxy for the quantity of primary producers in the surface ocean), PISCES model (including predictions of total iron concentrations comprising both mineral and organic Fe-compounds)
def Copernicus_download(year_bin,username=cfg_pw['Copernicus_user'],password=cfg_pw['Copernicus_pass'],output_path=path_to_git / 'data' / 'Environment' /'Copernicus' ,parameters_dict={'Sea Level Anomalies':'SLA','Biogeochemical forecast':'BGC','Particle backscatter':'BBP','Chlorophyll a':'CHL','Primary production':'PP'},SLA_dict={'SEALEVEL_GLO_PHY_L4_MY_008_047-TDS':'cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m'},BGC_dict={'GLOBAL_MULTIYEAR_BGC_001_029-TDS':'cmems_mod_glo_bgc_my_0.25_P1M-m','GLOBAL_ANALYSIS_FORECAST_BIO_001_028-TDS':'global-analysis-forecast-bio-001-028-monthly'},OC_dict={'OCEANCOLOUR_GLO_BGC_L4_MY_009_104-TDS':['cmems_obs-oc_glo_bgc-optics_my_l4-multi-4km_P1M','cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M','cmems_obs-oc_glo_bgc-pp_my_l4-multi-4km_P1M']}):
    """
    Objective: This function uses the python Copernicus API module motu (https://help.marine.copernicus.eu/en/articles/4796533-what-are-the-motu-apis#h_3d33beaafc) to download data from Copernicus (https://data.marine.copernicus.eu/products).\n
         :param year_bin: The year (str) of the data of interest.
         :param username: Copernicus username.
         :param password: Copernicus password.
         :param output_path: The path for downloaded files storage.
         :param parameters_dict: A dictionary of variables to download and their corresponding code. Default is all variables and their corresponding dictionaries
         :return: Returns Copernicus monthly file at output_path/parameter/period/*.nc
     """
    if 'SLA' in parameters_dict.values():# SLA
        print('Downloading Sea level anomalies. Please wait')
        output_file=Path(output_path).expanduser()/  'SLA' / year_bin
        output_file.mkdir(parents=True, exist_ok=True)
        path_to_file=output_file / "{}_{}_{}.nc".format(year_bin,'01-12', 'SLA')
        if not path_to_file.exists(): # Download netcdf
            os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-01-01 00:00:00" --date-max "{}-12-31 23:59:59" --variable sla --out-dir {} --out-name {} --user {} --pwd {}'.format(list(SLA_dict.keys())[0],list(SLA_dict.values())[0],year_bin,year_bin, path_to_file.parent, path_to_file.name,username,password))
        if not Path(output_file.parent/'README.xml').exists(): # Write README file
            os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --describe-product -o console --user {} --pwd {} --out-dir {} --out-name {}'.format(list(SLA_dict.keys())[0],list(SLA_dict.values())[0],username,password,output_file.parent,'README.xml'))
        if len(list(Path(output_file).expanduser().glob("*.parquet")))==0: # Save as parquet
            df=pd.concat(map(lambda path: nc_to_df_Copernicus(path, grid_res=1),natsorted(list(Path(output_file).expanduser().glob("*.nc")))))
            df.to_parquet(output_file / str('SLA_{}.parquet'.format(year_bin)),index=False)
            [file.unlink(missing_ok=True) for file in list(Path(output_file).expanduser().glob("*.nc"))]  # Delete original files to save storage space

    if 'BGC' in parameters_dict.values():#BGC 1993-2020
        print('Downloading Biogeochemical variables. Please wait')
        output_file = Path(output_path).expanduser()  / 'BGC' / year_bin
        output_file.mkdir(parents=True, exist_ok=True)
        if not Path(output_file.parent / 'README_.xml').exists():  # Write README file
            os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --describe-product -o console --user {} --pwd {} --out-dir {} --out-name {}'.format( list(BGC_dict.keys())[0], list(BGC_dict.values())[0], username, password, output_file.parent,'README.xml'))

        if pd.to_datetime(year_bin+'-01-01') in pd.date_range('1993','2020'):
            for month in np.arange(1,13,1):
                path_to_file = output_file / "{}_{}_{}.nc".format(year_bin, str(month).zfill(2), 'BGC')
                if not path_to_file.exists():
                    #os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-{}-01 00:00:00" --date-max "{}-{}-31 23:59:59" --depth-min 0.5057600140571594 --depth-max 4000 --variable chl --variable fe --variable no3 --variable nppv --variable o2 --variable phyc --variable po4 --variable si --out-dir {} --out-name {} --user {} --pwd {}'.format(list(BGC_dict.keys())[0],list(BGC_dict.values())[0],year_bin,str(month).zfill(2),year_bin,str(month).zfill(2), path_to_file.parent, path_to_file.name,username,password))
                    os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-{}-01 00:00:00" --date-max "{}-{}-31 23:59:59" --depth-min 0.5057600140571594 --depth-max 200 --variable fe --out-dir {} --out-name {} --user {} --pwd {}'.format(list(BGC_dict.keys())[0],list(BGC_dict.values())[0],year_bin,str(month).zfill(2),year_bin,str(month).zfill(2), path_to_file.parent, path_to_file.name,username,password))

            if len(list(Path(output_file).expanduser().glob("*.parquet"))) == 0: # Save as parquet file
                df = pd.concat(map(lambda path: nc_to_df_Copernicus(path, grid_res=1), natsorted(list(Path(output_file).expanduser().glob("*.nc")))))
                df.to_parquet(output_file / str('BGC_{}.parquet'.format(year_bin)), index=False)
                [file.unlink(missing_ok=True) for file in list(Path(output_file).expanduser().glob("*.nc"))]  # Delete original files to save storage space


        else: # BGC 2020-ongoing
            for month in np.arange(1, 13, 1):
                path_to_file = output_file / "{}_{}_{}.nc".format(year_bin, str(month).zfill(2), 'BGC')
                if not path_to_file.exists():
                   # os.system('python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-{}-01 00:00:00" --date-max "{}-{}-31 23:59:59" --depth-min 0.49402499198913574 --depth-max 4000 --variable chl --variable fe --variable no3 --variable nppv --variable o2 --variable phyc --variable po4 --variable si --out-dir {} --out-name {} --user {} --pwd {}'.format(list(BGC_dict.keys())[1],list(BGC_dict.values())[1],year_bin,str(month).zfill(2),year_bin,str(month).zfill(2), path_to_file.parent, path_to_file.name,username,password))
                    os.system('python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-{}-01 00:00:00" --date-max "{}-{}-31 23:59:59" --depth-min 0.49402499198913574 --depth-max 200  --variable fe --out-dir {} --out-name {} --user {} --pwd {}'.format(list(BGC_dict.keys())[1],list(BGC_dict.values())[1],year_bin,str(month).zfill(2),year_bin,str(month).zfill(2), path_to_file.parent, path_to_file.name,username,password))

    if 'BBP' in parameters_dict.values():# Ocean color
        print('Downloading Ocean color data. Please wait')
        ## Particulate backscatter: 09/1997-ongoing
        output_file = Path(output_path).expanduser() / 'BBP' / year_bin
        output_file.mkdir(parents=True, exist_ok=True)
        path_to_file = output_file / "{}_{}_{}.nc".format(year_bin,'01-12', 'BBP')
        if not path_to_file.exists():
            if pd.to_datetime(year_bin + '-01-01') in pd.date_range('1997-09-01', str(date.today())):
                cm.subset(username=username,password=password,force_download=True,minimum_depth=0,maximum_depth=1,dataset_id=list(OC_dict.values())[0][0],variables=["BBP"],start_datetime="{}-01-01T00:00:00".format(year_bin),end_datetime="{}-12-31T23:59:59".format(year_bin),output_filename = path_to_file.name,output_directory = path_to_file.parent)#os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-01-01 00:00:00" --date-max "{}-12-31 23:59:59" --variable BBP --out-dir {} --out-name {} --user {} --pwd {}'.format(list(OC_dict.keys())[0],list(OC_dict.values())[0][0],year_bin,year_bin, path_to_file.parent, path_to_file.name,username,password))
                #cm.subset(motu_api_request='python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-01-01 00:00:00" --date-max "{}-12-31 23:59:59" --variable BBP --out-dir {} --out-name {} --user {} --pwd {}'.format(list(OC_dict.keys())[0],list(OC_dict.values())[0][0],year_bin,year_bin, path_to_file.parent, path_to_file.name,username,password) )
        if not Path(output_file.parent / 'README.xml').exists():  # Write README file
            os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --describe-product -o console --user {} --pwd {} --out-dir {} --out-name {}'.format(list(OC_dict.keys())[0],list(OC_dict.values())[0][0], username, password, output_file.parent,'README.xml'))

    if 'CHL' in parameters_dict.values():    ## Chlorophyll-a and phytoplankton group concentration: 09/1997-ongoing
        output_file = Path(output_path).expanduser() / 'CHL' / year_bin
        output_file.mkdir(parents=True, exist_ok=True)
        if pd.to_datetime(year_bin + '-01-01') in pd.date_range('1997-09-01', str(date.today())):
            for month in np.arange(1, 13, 1):
                path_to_file = output_file / "{}_{}_{}.nc".format(year_bin, str(month).zfill(2), 'CHL')
                if not path_to_file.exists():
                    #os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-{}-01 00:00:00" --date-max "{}-{}-31 23:59:59" --variable CHL --variable DIATO --variable DINO --variable NANO --variable PICO --variable PROCHLO --variable PROKAR --variable MICRO --variable HAPTO --variable GREEN --out-dir {} --out-name {} --user {} --pwd {}'.format(list(OC_dict.keys())[0],list(OC_dict.values())[0][1],year_bin,str(month).zfill(2),year_bin,str(month).zfill(2), path_to_file.parent, path_to_file.name,username,password))
                    cm.subset(username=username,password=password,force_download=True,minimum_depth=0,maximum_depth=1,dataset_id=list(OC_dict.values())[0][1],variables=["CHL"],start_datetime="{}-01-01T00:00:00".format(year_bin),end_datetime="{}-12-31T23:59:59".format(year_bin),output_filename = path_to_file.name,output_directory = path_to_file.parent)#os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-{}-01 00:00:00" --date-max "{}-{}-31 23:59:59" --variable CHL --out-dir {} --out-name {} --user {} --pwd {}'.format(list(OC_dict.keys())[0], list(OC_dict.values())[0][1], year_bin, str(month).zfill(2), year_bin, str(month).zfill(2), path_to_file.parent, path_to_file.name, username, password))

            if not Path(output_file.parent / 'README.xml').exists():  # Write README file
                os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --describe-product -o console --user {} --pwd {} --out-dir {} --out-name {}'.format(list(OC_dict.keys())[0], list(OC_dict.values())[0][1], username, password, output_file.parent,'README.xml'))

            if len(list(Path(output_file).expanduser().glob("*.parquet")))==0:
               df=pd.concat(map(lambda path: nc_to_df_Copernicus(path, grid_res=1),natsorted(list(Path(output_file).expanduser().glob("*.nc")))))
               df.to_parquet(output_file / str('CHL_{}.parquet'.format(year_bin)),index=False)
               [file.unlink(missing_ok=True) for file in list(Path(output_file).expanduser().glob("*.nc"))]  # Delete original files to save storage space

    if 'PP' in parameters_dict.values():    ## Primary productivity: 09/2017-ongoing
        output_file = Path(output_path).expanduser() / 'PP' / year_bin
        output_file.mkdir(parents=True, exist_ok=True)
        path_to_file = output_file / "{}_{}_{}.nc".format(year_bin,'01-12', 'PP')
        if not path_to_file.exists():
            if pd.to_datetime(year_bin + '-01-01') in pd.date_range('1997-09-01', str(date.today())):
                cm.subset(username=username,password=password,force_download=True,minimum_depth=0,maximum_depth=1,dataset_id=list(OC_dict.values())[0][2],variables=["PP"],start_datetime="{}-01-01T00:00:00".format(year_bin),end_datetime="{}-12-31T23:59:59".format(year_bin),output_filename = path_to_file.name,output_directory = path_to_file.parent)#os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-01-01 00:00:00" --date-max "{}-12-31 23:59:59" --variable PP --out-dir {} --out-name {} --user {} --pwd {}'.format(list(OC_dict.keys())[0],list(OC_dict.values())[0][2],year_bin,year_bin, path_to_file.parent, path_to_file.name,username,password))
                #cm.subset(motu_api_request='python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --date-min "{}-01-01 00:00:00" --date-max "{}-12-31 23:59:59" --variable PP --out-dir {} --out-name {} --user {} --pwd {}'.format(list(OC_dict.keys())[0],list(OC_dict.values())[0][2],year_bin,year_bin, path_to_file.parent, path_to_file.name,username,password))
        if not Path(output_file.parent / 'README.xml').exists():  # Write README file
            os.system('python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --describe-product -o console --user {} --pwd {} --out-dir {} --out-name {}'.format(list(OC_dict.keys())[0],list(OC_dict.values())[0][2], username, password, output_file.parent,'README.xml'))

def nc_to_df_Copernicus(path,grid_res=1,integration_depth=[0,200]):
    df=pd.DataFrame()
    if path.exists():
        ds = xr.open_dataset(path, decode_times=False)
        if ('latitude' in ds.variables) & ('longitude' in ds.variables):
            ds=ds.rename({'latitude':'lat','longitude':'lon'})
        if ds.lon.attrs['valid_max']>180: #sla
            ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
        if len(ds.coords['time'])==2: # remove extra days of the next month for monthly nc files that have been extracted to the 31st
            ds=ds.isel(time=0)
        time_units=ds.time.units
        ds=ds.groupby_bins('lon', np.arange(-180,181,grid_res)).mean().groupby_bins('lat', np.arange(-90,91,grid_res)).mean()
        ds['lat_bins']=pd.IntervalIndex( pd.cut(np.arange(-91,90,grid_res),np.arange(-91,90,grid_res))).mid[1:]#-((-(1/(grid_res/2))*ds['lat'])//1)/(1/(grid_res/2))
        ds['lon_bins'] =pd.IntervalIndex( pd.cut(np.arange(-180,181,grid_res),np.arange(-180,181,grid_res))).mid[1:]#-((-(1/(grid_res/2))*ds['lon'])//1)/(1/(grid_res/2))
        ds=ds.rename({'lat_bins':'lat','lon_bins':'lon'})
        variables_dict={var:attr.attrs.get('units') for var,attr in ds.variables.items() if var not in ['crs','depth','lat','lon','time']}

        if (Path(path).parts[-3] == 'BGC') or ('depth' in list(ds.dims.keys())):  # Retrieve surface and integrated values
            #ds=ds.where(ds.depth<=200, drop=True)
            df=pd.concat(map(lambda variable:pd.merge(ds[variable].groupby_bins('depth', integration_depth).apply(lambda x: x.integrate(coord='depth')).to_dataframe().reset_index().rename(columns={'level_3':'integrated_depth',variable:'integrated_'+variable}),ds[variable].where(ds[variable].depth==ds[variable].depth.min(), drop=True).to_dataframe().reset_index().drop(columns=['depth']).rename(columns={variable:variable+'_surface'}),how='left',on=['time','lat','lon']),variables_dict.keys()),axis=1).reset_index(drop=True)
            df['datetime']=df.time.apply(lambda time:pd.to_datetime(time_units.split(' since ')[1]) + eval("relativedelta({}={})".format(time_units.split(' since ')[0],time)))
            df['Date_bin']=df['datetime'].dt.strftime('%Y-%m')
        # kwargs = {'mapping': 'y:lat,x:lon,fill:fe_surface','scale_fill_gradientn(trans="log10",colors={},guide=guide_colorbar(direction="vertical"))'.format(palette_temp): ''}
        # ggplot_raster(dataframe=df,output_path='',**kwargs)

        else:
            df = ds.to_dataframe().reset_index()
            df['datetime']=df.time.apply(lambda time:pd.to_datetime(time_units.split(' since ')[1]) + eval("relativedelta({}={})".format(time_units.split(' since ')[0],time)))
            df['Date_bin']=df['datetime'].dt.strftime('%Y-%m')
    # kwargs = {'mapping': 'y:lat,x:lon,fill:sla','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(palette_temp): ''}
    # kwargs = {'mapping': 'y:lat,x:lon,fill:CHL','scale_fill_gradientn(colors={},guide=guide_colorbar(direction="vertical"))'.format(palette_temp): ''}
    # ggplot_raster(dataframe=df,output_path='',**kwargs)
    return df

## World Ocean Atlas objective analyses: Sea surface temperature, nutrient (nitrates, phosphates, silicates, dissolved oxygen, apparent oxygen utilization) concentrations, mixed layer depth
def WOA_download(year_bin,output_path=path_to_git / 'data' / 'Environment' / 'WOA',url='https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/',grid='1.00',output_format='netcdf',parameters_dict={'Temperature':'t','Salinity':'s','Conductivity':'C','Mixed Layer Depth':'M','Dissolved Oxygen':'o','Apparent Oxygen Utilization':'A','Silicate':'i','Phosphate':'p','Nitrate':'n'}):
    """
    Objective: This function uses a webpage scraping module to download data from the World Ocean Atlas 2018 release (https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/).\n
        :param year_bin: The year (str) of the data of interest.
        :param output_path: The path for downloaded files storage.
        :param url: The source URL to scroll in order to download World Ocean Atlas data. Default is 2018 WOA release.
        :param output_format: The output format. One of ascii, csv, shape (ArcGIS),or netcdf. Default is netcdf
        :param grid: The output grid. One of ['5deg', '1.00', '0.25'] degree(s). Default is 1 to match PSSdb products
        :param parameters_dict: A dictionary of variables to download and their corresponding code. Default is all variables: {'Temperature':'t','Salinity':'s','Conductivity':'C','Mixed Layer Depth':'M','Dissolved Oxygen':'o','Apparent Oxygen Utilization':'A','Silicate':'i','Phosphate':'p','Nitrate':'n'}
    :return: Returns WOA climatology file at output_path/parameter/period/*.nc
    """
    links_url={'t':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/temperature/{}/'.format(output_format),
              'C':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/conductivity/{}/'.format(output_format),
              'D':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/density/{}/'.format(output_format),
              's':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/salinity/{}/'.format(output_format),
              'M':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/mld/{}/'.format(output_format),
              'o':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/oxygen/{}/'.format(output_format),
              'A':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/AOU/{}/'.format(output_format),
              'i':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/silicate/{}/'.format(output_format),
              'p':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/phosphate/{}/'.format(output_format),
              'n':'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/nitrate/{}/'.format(output_format)}
    with  HTMLSession() as session:
        # GET request
        res = session.get(url)
        # for javascript driven website
        # res.html.render()
        soup = BeautifulSoup(res.html.html, "html.parser")
        if not Path(output_path/ 'README.xml').exists():  # Write README file
            rsp = session.get(url='https://www.ncei.noaa.gov/access/world-ocean-atlas-2018f/', headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
            df_soup = (BeautifulSoup(rsp.content, "lxml").find_all('dl')[0]).find_all('dd')
            df_data=pd.concat(map(lambda ind: pd.DataFrame(df_soup[ind].text.replace('\n','').replace(")",'').split(' (')).T.set_axis(['Variable','Unit'],axis='columns') ,np.arange(0,len(df_soup))))
            df_data['Unit']=np.where(df_data.Unit.str.contains('mol'),df_data.Unit +' * 1025 kg/m3 equivalent to '+df_data.Unit.str.rsplit('/').str[0]+'/m3',df_data.Unit)
            rsp.close()
            output_path.mkdir(exist_ok=True)
            with open(Path(output_path/ 'README.xml'), 'w') as f:
                f.write(df_data.reset_index(drop=True).to_xml(namespaces={"": "https://www.ncei.noaa.gov/access/world-ocean-atlas-2018f/"},prefix="",index=False,attr_cols=['Variable', 'Unit'],parser='etree'))

    form=soup.find_all("form")[0]
    action = form.attrs.get("action").lower()
    href={item.text:item.attrs['href'] for item in soup.find_all("a") if item.text in parameters_dict.keys()}
    if len(href):
        for variable,url_parameter in href.items():
            form=get_all_forms(url+url_parameter)
            inputs=list(map(lambda x: get_form_details(x),form))
            params = join_with(list, [{input['name']:input['value']} for input in inputs[1]['inputs']])
            # Replace with function arguments
            params['format']=output_format if output_format in params['format'] else params['format'][-1]
            params['grid']=grid if grid in params['grid'] else params['grid'][-1]
            params['parameter']=url_parameter.split('?parameter=')[1]
            params['climatology'] ='mn'
            params['decades']={'-'.join([option.text.replace('Averaged decades years','all').replace(' years','').split('-')[0],option.text.replace('Averaged decades years','all').replace(' years','').split('-')[1]]) if (len(option.text.replace('Averaged decades years', 'all').replace(' years', '').split('-'))==2) and (len(option.text.replace('Averaged decades years', 'all').replace(' years', '').split('-')[1])== 4) else '-'.join([option.text.replace('Averaged decades years', 'all').replace(' years', '').split('-')[0], option.text.replace('Averaged decades years', 'all').replace(' years', '').split('-')[0][0:2]+option.text.replace('Averaged decades years', 'all').replace(' years', '').split('-')[1]]) if (len(option.text.replace('Averaged decades years', 'all').replace(' years', '').split('-'))==2) and (len(option.text.replace('Averaged decades years', 'all').replace(' years', '').split('-')[1])== 2)  else option.text.replace('Averaged decades years','all').replace(' years','').split('-')[0]:option['value'] for option in form[1].find(id='decades').find_all('option')} if form[1].find(id='decades') is not None else {'all': 'decav'}
            # Select time period:
            period_code=[value for period,value in params['decades'].items() if  (len(period.split('-'))==2) and (year_bin in pd.date_range(period.split('-')[0],period.split('-')[1],freq='Y').year)]
            period=[period for period,value in params['decades'].items() if  (len(period.split('-'))==2) and (year_bin in pd.date_range(period.split('-')[0],period.split('-')[1],freq='Y').year)]

            if len(period)==0:
                period_code=params['decades']['all'] if 'all' in params['decades'].keys() else params['decades'][list(params['decades'].keys())[0]]
                period='All'
            else:
                period_code=period_code[0]
                period=period[0]

            params['decades'] =period_code
            path_to_zip=output_path/ (variable.replace(' ','_')) / period
            path_to_zip.mkdir(exist_ok=True,parents=True)
            if len(list(path_to_zip.expanduser().glob('*_01.nc')))==0:
                with  HTMLSession() as session:
                    rsp = session.post(url+url_parameter.split('?')[0], data=params, headers={'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36'})
                    soup = BeautifulSoup(rsp.html.html, "html.parser")
                    monthly_links={item.text: item.attrs['href'] for item in soup.find_all("a") if item.text in list(map(lambda x: params['parameter']+str(x).zfill(2)+'_01.nc', range(1, 13)))}
                    with tqdm(desc='Downloading WOA monthly climatology for {}'.format(variable), total=len(monthly_links), bar_format='{desc}{bar}', position=0,leave=True) as bar:
                        for month, url_month in monthly_links.items():
                            ftp_url=links_url[params['parameter']]+os.path.join(*Path(url_month.split('catalog.html?dataset=')[1]).parts[3:])
                            file = session.get(ftp_url)
                            bar.set_description('Downloading WOA monthly climatology for {}/{}'.format(variable,Path(url_month.split('catalog.html?dataset=')[1]).parts[-1]), refresh=True)
                            with open(path_to_zip /month, "wb") as fd:
                                for a_chunk in file.iter_content():  # Loop over content, i.e. eventual HTTP chunks
                                    fd.write(a_chunk)
                            ok = bar.update(n=1)

                print('WOA monthly climatology for {} (grid: {}, period:{}) complete. Check at: \n{}'.format(variable,grid,period,path_to_zip ))
            if len(list(path_to_zip.expanduser().glob('*.parquet'))) == 0:
                try:
                    print('Concatenating and saving data for {}'.format(variable))
                    df=pd.concat(map(lambda path:nc_to_df_woa(path),natsorted(list(path_to_zip.expanduser().glob('*.nc'))))).reset_index(drop=True)
                    df.to_parquet(path_to_zip /str( 'WOA_'+path_to_zip.parts[-2].lower()+'.parquet'),index=False)
                    [file.unlink(missing_ok=True) for file in list(path_to_zip.expanduser().glob('*.nc'))]  # Delete original files to save storage space

                except:
                    print('')


            else:
                print('WOA monthly climatology for {} (grid: {}, period:{}) already downloaded. Skipping'.format(variable,grid,period))

def nc_to_df_woa(path):
    df=pd.DataFrame()
    path=Path(path)
    if path.exists():
        ds = xr.open_dataset(path, decode_times=False)
        variables_dict = {var: attr.attrs.get('units') for var, attr in ds.variables.items() if var not in ['depth', 'lat', 'lon', 'time','crs','lat_bnds','lon_bnds','depth_bnds','climatology_bounds']}
        time_units=ds.time.units
        ds = ds.where(ds.depth <= 200, drop=True)
        df = pd.concat(map(lambda variable:( df_integral:=1025*ds[variable].integrate(coord='depth').to_dataframe().reset_index().rename(columns={variable:  'value_integral'}) if 'moles_concentration' in ds[variable].standard_name else ds[variable].mean(['depth']).to_dataframe().reset_index().rename(columns={variable:  'value_integral'}),df_surface:=1025*ds[variable][:, 0, :, :].to_dataframe().reset_index().drop(columns=['depth']).rename(columns={variable:  'value_surface'}).assign(variable=ds[variable].standard_name,unit=ds[variable].units.replace('per_kilogram','per_square_meters')) if 'moles_concentration' in ds[variable].standard_name else ds[variable][:, 0, :, :].to_dataframe().reset_index().drop(columns=['depth']).rename(columns={variable:  'value_surface'}).assign(variable=ds[variable].standard_name,unit=ds[variable].units),pd.merge(df_integral,df_surface, how='left', on=['time', 'lat', 'lon']))[-1],[variable for variable in variables_dict.keys() if '_an' in variable]), axis=1).reset_index(drop=True)
        df=df.rename(columns={df.columns[-1]: 'value', 'lat':'latitude','lon':'longitude'})
        df['variable']=(ds.title.split(': ')[1]).split(' ')[0]
        df['month']=(ds.title.split(': ')[1]).split(' ')[1]
        df['Period']=(ds.title.split(': ')[1]).split(' ')[2]
        # ((pd.to_datetime('1955-01-01 00:00:00')+relativedelta(months=372))-pd.to_datetime('1955-01-01 00:00:00')) / np.timedelta64(1, 'M')
        #df['datetime'] = df.time.astype(int).apply(lambda month: pd.to_datetime('1955-01-01 00:00:00') + relativedelta(months=month))
    return df

## Spearman correlation coefficients
def spearman_correlation(x, y):
    try:
        return spearmanr(x, y)[1]
    except:
        return np.nan



## XGBoost: Boosted Regression Trees model used to predict a set of continuous, correlated variables (multioutput - NBSS parameters) as a function of explanatory variables (all environemental factors downloaded above)

def load_xgbmodel(path):
    model = xgb.XGBRegressor()
    model.load_model(path)
    return model

def nbss_predict(df,selection={'PFT':'Mesophytoplankton'},var_x=['AOU', 'FSLE','SLA','MLD','CHL','Phosphates','Silicates','Nitrates','Temperature','Iron','PAR','Aerosol_thickness','DUWT001'],var_y=['Intercept','Slope','Min_size','Max_size'],figures_path=path_to_git / 'Plots' /'PSSdb_paper'):
    """
    Objective: This function uses the xgboost module to predict NBSS parameters (var_y) as function of environmental factors (var_x).\n
        :param df: The merged dataframe including the response and explanatory variables.
        :param selection: Filter used to select specific taxonomic group and/or instrument in the dataframe.
        :param var_x: List of explanatory variables (unscaled with native units).
        :param var_y: List of response variable(s) (unscaled with native units. A log-transformation will be applied to ensure interpcets and size ranges are always positive). XGboost allows for multiouput response variable.
        :param figures_path: The path for plots storage. Plots automatically saved include scatterplots of observations vs predictions, model tree, or barplot showing var_x feature importance.
    :return: Returns fitted model at output_path/parameter/period/*.nc
    """
    if len(selection):
        df=df.query('&'.join(['({}=="{}")'.format(key,value) if isinstance(value,str) else '({}.isin({})==True)'.format(key,value) for key,value in selection.items()]))
    X,y=df.loc[df.dropna(subset=var_x+var_y).index,var_x],df.loc[df.dropna(subset=var_x+var_y).index,var_y].transform(lambda x: np.log10(x) if x.name in ['Min_size','Max_size'] else x)
    seed=5
    X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=.2,random_state=seed)
    # create model instance
    model=xgb.XGBRegressor(n_estimators=100, max_depth=7, eta=0.1, subsample=0.7, colsample_bytree=0.8)
    eval_set = [(X_test, y_test)]
    print('Fitting Gradient Boosted Regression Trees model to dataframe. Please wait')
    # define model evaluation method
    cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=seed)
    # evaluate model
    #scores = cross_val_score(model, X, y, cv=cv, n_jobs=-1)
    #scores = np.absolute(scores)
    #scores.mean(),scores.std()
    param_n = {'n_estimators': range(20, 100, 10),'max_depth':range(3,10,1)}
    model_grid = GridSearchCV(estimator=xgb.XGBRegressor(learning_rate=0.1,  eta=0.1, subsample=0.7, colsample_bytree=0.8),param_grid=param_n,cv=cv)
    model_grid.fit(X_train, y_train.loc[:,'Intercept'], eval_metric=[ "logloss"], eval_set=[(X_test, y_test.loc[:,'Intercept'])], verbose=False)
    model_grid.best_params_
    model = xgb.XGBRegressor(n_estimators=model_grid.best_params_['n_estimators'], max_depth=model_grid.best_params_['max_depth'], eta=0.1, subsample=0.7, colsample_bytree=0.8)

    model=model.fit(X_train, y_train, eval_metric=["error", "logloss"], eval_set=eval_set, verbose=False) # fit model model.fit(X, y)
    """
    xgb.plot_tree(model)
    fig = plt.gcf()
    #fig.set_size_inches(30,15)
    if len(str(figures_path)):
        figures_path.expanduser().mkdir(exist_ok=True)
        print('Saving {} tree model to {}'.format(''.join(list(selection.values())[0]),str(figures_path)))
        fig.savefig(fname='{}/Model_tree_{}.svg'.format(str(figures_path),'_'.join(['_'.join(selection) if isinstance(selection,list) else selection for selection in list(selection.values())])),limitsize=False, dpi=600, bbox_inches='tight')
    """
    # Save model into JSON format
    Path(path_to_git / 'data' / 'NBSS_model' / 'Predictions_output').mkdir(exist_ok=True)
    print('Saving {} tree model to {}'.format('_'.join(['_'.join(selection) if isinstance(selection,list) else selection for selection in list(selection.values())]), str(Path(path_to_git / 'data' / 'NBSS_model' / 'Predictions_output'))))
    model.save_model("{}/xgboost_{}.json".format(Path(path_to_git/ 'data' / 'NBSS_model' /'Predictions_output'),'_'.join(['_'.join(selection) if isinstance(selection,list) else selection for selection in list(selection.values())])))
    # make predictions for test data
    test_pred = model.predict(X_test)
    sse = np.sum((test_pred - y_test) ** 2)
    sst = np.sum((y_test - np.mean(y_test)) ** 2)
    rsquared = 1 - sse / sst

    df_pred=pd.merge(y_test.transform(lambda x: 10**x if x.name in ['Min_size','Max_size'] else x).melt(value_vars=y.columns,value_name='observed_value'),pd.DataFrame({y_test.columns[0]:test_pred[:,[0]].flatten(),y_test.columns[1]:test_pred[:,[1]].flatten(),y_test.columns[2]:10**test_pred[:,[2]].flatten(),y_test.columns[3]:10**test_pred[:,[3]].flatten()}).melt(value_vars=y_test.columns,value_name='predicted_value'), left_index=True, right_index=True).rename(columns={'variable_x':'variable'}).drop(columns='variable_y').assign(set='test')
    df_pred = pd.concat([df_pred,pd.merge(y.transform(lambda x: 10 ** x if x.name in ['Min_size', 'Max_size'] else x).melt(value_vars=y.columns,value_name='observed_value'),pd.DataFrame({y.columns[0]: y_pred[:, [0]].flatten(), y.columns[1]: y_pred[:, [1]].flatten(),y.columns[2]: 10 ** y_pred[:, [2]].flatten(),y.columns[3]: 10 ** y_pred[:, [3]].flatten()}).melt(value_vars=y.columns,value_name='predicted_value'),left_index=True, right_index=True).rename(columns={'variable_x': 'variable'}).drop(columns='variable_y').assign(set='all')],axis=0)

    df_pred['variable']=pd.Categorical(df_pred.variable,categories=df_pred.variable.unique(),ordered=True)
    df_pred_summary=df_pred.groupby('variable').apply(lambda x: pd.Series({'k':(X.shape[1] ,model_grid.best_params_['max_depth'],model_grid.best_params_['n_estimators']),'r2':'R$^{2}$: '+str(np.round(r2_score(y[x.variable.unique()[0]], y_pred[:,[x.variable.cat.codes.unique()[0]]]),2)),'logloss':model.evals_result()['validation_0']['logloss'][-1],'n':len(x.observed_value),'x':x.observed_value.min(),'y':x.predicted_value.max(),'R2':'R$^{2}$: '+str(np.round(statsmodels.formula.api.ols(formula='predicted_value~observed_value', data=x).fit().rsquared,2)),'AIC':statsmodels.formula.api.ols(formula='predicted_value~observed_value', data=x).fit().aic,'BIC':statsmodels.formula.api.ols(formula='predicted_value~observed_value', data=x).fit().bic,'Loglik':statsmodels.formula.api.ols(formula='predicted_value~observed_value',data=x).loglike(statsmodels.formula.api.ols(formula='predicted_value~observed_value',data=x).fit().params)})).reset_index()
    plot=(ggplot(df_pred[df_pred.observed_value!=0])+facet_wrap('~variable',ncol=1,scales='free')+geom_abline(slope=1,intercept=0)+geom_point(mapping=aes(x='observed_value',y='predicted_value',fill='set'), stat='identity')+labs(y=r'Predictions', x=r'Observations',size='') +
         geom_text(data= df_pred_summary,mapping=aes(x='x',y='y',label='r2'),ha='left',va='top')+
         scale_fill_manual(values={'test':'black','all':'black'})+
         theme( panel_grid=element_blank(),legend_position='bottom',
            panel_background=element_rect(fill='white'),legend_background=element_rect(fill='white'),
            strip_background=element_rect(fill='white'),
            panel_border=element_rect(color='#222222'),
            legend_title=element_text(family='serif', size=10), legend_text=element_text(family='serif', size=10,rotation=0),
            axis_title=element_text(family='serif', size=10,linespacing=1), axis_text_x=element_text(family='serif',ha='right', size=10,linespacing=0,rotation=0),
            axis_text_y=element_text(family='serif', size=10, rotation=90,linespacing=1),
            plot_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)
    plot[0].set_size_inches(3,10)
    if len(str(figures_path)):
        figures_path.expanduser().mkdir(exist_ok=True)
        print('Saving {} test set regression to {}'.format('_'.join(['_'.join(selection) if isinstance(selection,list) else selection for selection in list(selection.values())]),str(Path(path_to_git / 'Data' / 'NBSS' / 'Predictions_output'))))
        plot[0].savefig(fname='{}/Model_predictions_vs_observations_{}.svg'.format(str(figures_path),'_'.join(['_'.join(selection) if isinstance(selection,list) else selection for selection in list(selection.values())])),limitsize=False, dpi=600, bbox_inches='tight')

    # retrieve performance metrics
    metrics = model.evals_result()
    df_metrics=pd.DataFrame({'Epochs':np.arange(len(np.exp(metrics['validation_0']['logloss']))),'Loss':np.exp(metrics['validation_0']['logloss'])})
    perm_importance = permutation_importance(model, X_test, y_test)
    dict_gain=dict(zip(model.feature_names_in_,model.feature_importances_))#
    dict_gain=dict(zip(model.feature_names_in_,np.fmax(0,perm_importance.importances_mean)))
    #dict_gain=model.get_booster().get_score(importance_type='total_gain')
    if any(pd.Series(dict_gain.keys()).isin(['time','coordinates_latitude','coordinates_longitude'])):
        dict_gain=project(dict_gain, list(filter(lambda var: var not in ['time','coordinates_latitude','coordinates_longitude'], var_x)))
    features_importance=dict(zip(list(dict_gain.keys()),np.array(list(dict_gain.values()))/np.array(list(dict_gain.values())).sum()))
    df_predictors=pd.DataFrame.from_dict(model.get_booster().get_score(importance_type='gain'),orient='index').reset_index().rename(columns={'index':'variable',0:'scores'}).sort_values(['scores'],ascending=False)
    df_predictors=pd.merge(pd.DataFrame.from_dict(features_importance,orient='index').reset_index().rename(columns={'index':'variable',0:'scores_permutation'}),pd.DataFrame.from_dict(dict(zip(var_x,model.feature_importances_)),orient='index').reset_index().rename(columns={'index':'variable',0:'scores_gain'}),how='inner',on='variable').sort_values(['scores_permutation'],ascending=False)
    plot=(ggplot(df_predictors)+scale_x_discrete(limits=df_predictors.variable.values)+geom_col(mapping=aes(x='variable',y='scores_permutation+1e-02'),fill='black', stat='identity')+labs(y=r'Feature importances', x=r'',size='') +
         theme( panel_grid=element_blank(),legend_position='bottom',
            panel_background=element_rect(fill='white'),legend_background=element_rect(fill='white'),
            strip_background=element_rect(fill='white'),
            panel_border=element_rect(color='#222222'),
            legend_title=element_text(family='serif', size=10), legend_text=element_text(family='serif', size=10,rotation=0),
            axis_title=element_text(family='serif', size=10,linespacing=1), axis_text_x=element_text(family='serif',ha='right', size=10,linespacing=0,rotation=45),
            axis_text_y=element_text(family='serif', size=10, rotation=0,linespacing=1),
            plot_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)
    plot[0].set_size_inches(6,3)
    if len(str(figures_path)):
        figures_path.expanduser().mkdir(exist_ok=True)
        print('Saving {} features importance to {}'.format('_'.join(['_'.join(selection) if isinstance(selection,list) else selection for selection in list(selection.values())]),str(Path(path_to_git / 'Data' / 'NBSS' / 'Predictions_output'))))
        plot[0].savefig(fname='{}/Model_scores_{}.svg'.format(str(figures_path),'_'.join(['_'.join(selection) if isinstance(selection,list) else selection for selection in list(selection.values())])),limitsize=False, dpi=600, bbox_inches='tight')

    return model, df_pred_summary

## Execute this script to start downloading specific datasets - main
if __name__ == '__main__':
    print('Downloading environmental climatologies assistant. \nPlease enter the following information (attention: method sensitive to empty space and lower/upper case. Do not use quotation)')
    funcs_args=input('Type:\n 1 for downloading data from Copernicus\n 2 for downloading data from World Ocean Atlas\n 3 for downloading data from NASA Ocean color\n 4 for downloading data from AVISO\n 5 for downloading PSD data from Kostadinov et al. 2023 (https://doi.pangaea.de/10.1594/PANGAEA.939863?format=html#download)\n')
    if funcs_args == '1':
        year_bin_args=input('Year for data coverage (%yyyy, use semi-colon for multiple years):')
        username_args = input('Copernicus username:')
        password_args = input('Copernicus password:')
        params_args=input('Copernicus data: all or SLA (Sea Level Anomalies) BGC (Biogeochemistry forecast) BBP (Particle backscatter), CHL (Chlorophyll a concentration), PP (Primary productivity)\nUse semi-colon to separate multiple parameters')

        parameters_dict = {key:item for key,item in {'Sea Level Anomalies': 'SLA', 'Biogeochemical forecast': 'BGC','Particle backscatter': 'BBP', 'Chlorophyll a': 'CHL', 'Primary production': 'PP'}.items() if item in params_args.split(';')} if params_args!='all' else {'Sea Level Anomalies': 'SLA', 'Biogeochemical forecast': 'BGC','Particle backscatter': 'BBP', 'Chlorophyll a': 'CHL', 'Primary production': 'PP'}
        if len(year_bin_args) :
            if len(username_args) and len(password_args):
                for years in year_bin_args.split(';'):
                    Copernicus_download(year_bin=years, username=username_args, password=password_args,output_path=path_to_git / 'data' / 'Environment' / 'Copernicus',parameters_dict=parameters_dict)
            else:
                for years in year_bin_args.split(';'):
                    Copernicus_download(year_bin=years,output_path=path_to_git / 'data' / 'Environment' / 'Copernicus',parameters_dict=parameters_dict)

        quit()
    if funcs_args == '2':
        year_bin_args = input('Year for data coverage (%yyyy):')
        if len(year_bin_args):
            WOA_download(year_bin=year_bin_args, output_path=path_to_git / 'data' / 'Environment'/ 'WOA',url='https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/', grid='1.00',output_format='netcdf')
            quit()
    if funcs_args == '3':
        year_bin_args = input('Year for data coverage (%yyyy, use semi-colon for multiple years):')
        params_args = input('Variable of interest: OC for Ocean color and AD for Aerosol diagnostic')
        appkeys_args = input('NASA App key:')
        for years in year_bin_args.split(';'):
            if len(appkeys_args):
                NASA_download(year_bin=years, output_path=path_to_git / 'data' / 'Environment' / 'NASA',parameters=[params_args] if len(params_args) else ['OC','AD'] , appkey=appkeys_args)
            else:
                NASA_download(year_bin=years, output_path=path_to_git / 'data' / 'Environment' / 'NASA',parameters=[params_args] if len(params_args) else ['OC','AD'])
        quit()
    if funcs_args == '4':
        year_bin_args = input('Year for data coverage (%yyyy, use semi-colon for multiple years):')
        username_args = input('AVISO username:')
        password_args = input('AVISO password:')
        for years in year_bin_args.split(';'):
            if len(username_args) and len(password_args):
                AVISO_download(years, username=username_args, password=password_args, output_path=path_to_git / 'data' / 'Environment' / 'AVISO')
            else:
                AVISO_download(years, output_path=path_to_git / 'data' / 'Environment' / 'AVISO')
        quit()
    if funcs_args == '5':
        year_bin_args = input('Year for data coverage (unique string, %yyyy, use semi-colon for multiple years):')
        if len(year_bin_args):
            for years in year_bin_args.split(';'):
                PSD_download(year_bin=years, output_path=path_to_git / 'data' / 'Environment' / 'PSD')
            quit()