## Objective: This script generates a series of plots to check and select PSSdb taxa-specific products for habitat modelling

## Functions: All functions required to run this step are included in funcs_predict_NBSS.py

## Modules:
from natsort import natsorted
import pandas as pd
import numpy as np # Arrays
import datetime
from pathlib import Path # Handling of path object
# Main functions to download environmental datasets and predict NBSS:
try:
    from funcs_predict_NBSS import *
except:
    from scripts.funcs_predict_NBSS import *

## Workflow starts here:
path_to_datafile=path_to_git / 'data' /'NBSS_ver_12_2023'
nbss=pd.concat(map(lambda path:pd.read_csv(path).assign(Instrument=path.name.split('_')[0]),list(Path(path_to_datafile /'PFT').expanduser().glob('*_1a_Size-distribution*.csv')))).reset_index(drop=True)
nbss_class=pd.concat(map(lambda path:pd.read_csv(path).assign(Instrument=path.name.split('_')[0]),list(Path(path_to_datafile /'Class').expanduser().glob('*_1a_Size-distribution*.csv')))).reset_index(drop=True)

# Filter out observations that do not span at least 80% of the euphotic layer & where mesophytoplankton is dominated by chains of neritic diatoms
group = ['Instrument','latitude','longitude','ocean','year','month']
nbss_stat=nbss_class.groupby(group+['Taxon']).normalized_biovolume_mean.sum().reset_index()
nbss_stat=nbss_stat.query('(Instrument=="UVP") & (Taxon.isin(["Bacillariophyceae","Cyanophyceae"]))').groupby(['latitude','longitude','Taxon']).mean().reset_index()#.groupby(['latitude','longitude']).apply(lambda x:x.loc[x.normalized_biovolume_mean.idxmax(),'Taxon']=='Bacillariophyceae')
nbss_stat['selection']=nbss_stat.groupby(['latitude','longitude']).apply(lambda x: pd.DataFrame({'selected':x.loc[x.normalized_biovolume_mean.idxmax(), 'Taxon'] == 'Cyanophyceae'},index=x.index)).reset_index()[['selected']]

plot = (ggplot(nbss_stat) +
        geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), color='black', fill='black') +
        geom_point(mapping=aes(x='longitude', y='latitude', size='normalized_biovolume_mean',fill='Taxon',color='selection'), alpha=0.1) +
        labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)', fill='', title='') +
        scale_size_area(trans='log10',max_size=3, na_value=1) +scale_color_manual(values={True:'black',False:'blue'})+
        theme_paper +
        coord_cartesian(xlim=[-180, 180], ylim=[-90, 90], expand=False)).draw(show=True)
nbss=nbss.drop(index=nbss[((nbss.Instrument=='UVP') & (nbss.PFT=='Mesophytoplankton') & (nbss[['latitude','longitude']].apply(tuple, axis=1).isin(nbss_stat.query('selection==False')[['latitude','longitude']].apply(tuple, axis=1))))].index).reset_index(drop=True)
nbss=nbss[(nbss.groupby(group).apply(lambda x:x[['min_depth','max_depth']].diff(axis=1)[['max_depth']]>=200*0.8)['max_depth'])].reset_index(drop=True)
nbss=pd.merge(nbss,nbss.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename( {'index': 'Group_index'}, axis='columns'),how='left',on=group)

# Append observations without records: null abundances should be provided to fit habitat models
multiindex = pd.MultiIndex.from_product([list(nbss.astype({column: 'category' for column in ['Group_index', 'biovolume_size_class', 'PFT']})[column].cat.categories) for column in ['Group_index', 'biovolume_size_class', 'PFT']],names=['Group_index', 'biovolume_size_class', 'PFT'])
nbss=nbss.set_index(['Group_index', 'biovolume_size_class', 'PFT']).drop(index=nbss.set_index(['Group_index', 'biovolume_size_class', 'PFT']).index[nbss.set_index(['Group_index', 'biovolume_size_class', 'PFT']).index.duplicated()==True]).reset_index()
nbss = pd.merge(pd.merge( nbss.drop(columns=['Instrument','latitude','longitude','ocean','year','month','min_depth', 'max_depth', 'n','equivalent_circular_diameter_mean']).set_index(['Group_index', 'biovolume_size_class', 'PFT']).reindex(multiindex, fill_value=pd.NA).reset_index(),nbss.drop_duplicates(['Group_index'])[['Instrument','latitude','longitude','ocean','year','month','Group_index']], how='left', on=['Group_index']),nbss[['biovolume_size_class','equivalent_circular_diameter_mean']].drop_duplicates(['biovolume_size_class']).reset_index(drop=True), how='left', on=['biovolume_size_class']).sort_values(['Group_index', 'PFT', 'biovolume_size_class']).reset_index(drop=True)

# Select only colonial N2-fixers, Crustacea, or Rhizaria
group_taxa=['Mesophytoplankton','Rhizaria','Crustaceans']
nbss_subset=nbss[(nbss.PFT.isin(group_taxa))].reset_index(drop=True)

# Use ggplot to plot group-specific Normalized Biovolume Size Spectrum and NBSS integral histograms:
plot = (ggplot(nbss_subset.dropna(subset=['normalized_biovolume_mean'])) +
        facet_grid('PFT~Instrument',scales='free_x') +
        geom_line(mapping=aes(x='equivalent_circular_diameter_mean', y='normalized_biovolume_mean*1e+03', group='Group_index'), size=0.17, alpha=1) + #
        geom_point(mapping=aes(x='equivalent_circular_diameter_mean', y='normalized_biovolume_mean*1e+03', group='Group_index'), size=0.07, alpha=1) +  #
        labs(x='Equivalent circular diameter (µm)',y='Normalized Biovolume Size Spectrum (mm$^{3}$ m$^{-3}$ mm$^{-3}$)', title='',colour='') +
        scale_x_log10(limits=[1e+02, 5e+04], expand=(0, 0), breaks=np.multiply( 10 ** np.arange(np.floor(np.log10(1e+02)), np.ceil(np.log10(5e+04)), step=1).reshape( int((np.ceil(np.log10(1e+02)) - np.floor(np.log10(5e+04)))), 1), np.arange(1, 10, step=1).reshape(1, 9)).flatten(), labels=lambda l: ['x10$^%s$' % round(np.floor(np.log10(v))) if ((v / (10 ** np.floor(np.log10(v)))) == 1) else '%s' % round( v / (10 ** np.floor(np.log10(v)))) for v in l]) +
        scale_y_log10(limits=[1e-01, 1e+03], breaks=10 ** np.arange(np.floor(np.log10(1e-05)) - 1, np.ceil(np.log10(1e+03)), step=1), labels=lambda l: ['10$^{%s}$' % int(np.log10(v)) if (np.log10(v)) / int(np.log10(v)) == 1 else '10$^{0}$' if v==1 else '' for v in l]) +
        theme_paper).draw(show=True)
plot.set_size_inches(4.5,5)
plot.savefig(fname='{}/figures/PSSdb_paper/NBSS_products_1a_PFT.svg'.format(str(path_to_git)), dpi=300, bbox_inches='tight')
plot = (ggplot(nbss_subset.groupby(group+['Group_index','PFT']).apply(lambda x:pd.Series({'Total_abundance':(1e+03)*x.normalized_biovolume_mean.sum(min_count=0)})).reset_index()) +
        facet_grid('PFT~Instrument',scales='free_x') +
        geom_histogram(mapping=aes(x='Total_abundance+1',y=after_stat('count')), alpha=1,binwidth=0.5,fill='black') + #
        stat_bin(mapping=aes(x='Total_abundance+1',y=after_stat('count'),label=after_stat('count'),group=1),geom='text', alpha=1,binwidth=0.5,angle=-90, size=4, format_string='{:.1f}')+ #
        coord_flip()+
        labs(y='Number of bins',x='NBSS integral (mm$^{3}$ m$^{-3}$ mm$^{-3}$)', title='',colour='') +
         #scale_x_log10(limits=[5e+02,5e+04],expand=(0,0),breaks=np.multiply(10 ** np.arange(np.floor(np.log10(3e+02)), np.ceil(np.log10(1e+04)),step=1).reshape(int( (np.ceil(np.log10(3e+02)) - np.floor(np.log10(1e+04)))), 1),np.arange(1, 10, step=1).reshape(1, 9)).flatten(),labels=lambda l: ['%s' % round(v) if ((v/(10**np.floor(np.log10(v))))==1) | ((v/(10**np.floor(np.log10(v))))==5) else '' for v in l]) +
        scale_x_log10(labels=lambda l: ['10$^{%s}$' % int(np.log10(v-1)) if (v-1)>0 else 'ND' if v-1==0 else '' for v in l]) +
        theme_paper).draw(show=True)
plot.set_size_inches(1.5,5)
plot.savefig(fname='{}/figures/PSSdb_paper/Histogram_NBSS_integral.svg'.format(str(path_to_git)), dpi=300, bbox_inches='tight')


# Extract NBSS parameters - no thresholding is performed since the thresholding was already applied for products generation
uncertainty_threshold=1#0.2
bins=1#3
nbss_summary=nbss_subset.assign(NBSS=lambda x:1e+03*x.normalized_biovolume_mean).rename(columns={'equivalent_circular_diameter_mean':'size_class_mid'}).groupby(group+['Group_index','PFT'],dropna=False).apply(lambda x: pd.DataFrame(NBSS_regression(nbss=x,selection='indexing',threshold_count=uncertainty_threshold,threshold_size=uncertainty_threshold,n_bins=bins),index=[0])).reset_index().drop(columns='level_{}'.format(len(group+['Group_index','PFT'])))
nbss_summary=nbss_summary.assign(Absolute_intercept=np.where(nbss_summary.Intercept!=0,10**nbss_summary.Intercept,0))
data=pd.concat([nbss_summary[(nbss_summary.PFT=='Crustaceans') & (nbss_summary.Instrument=='Scanner')],nbss_summary[(nbss_summary.PFT.isin(['Rhizaria','Mesophytoplankton'])) & (nbss_summary.Instrument=='UVP')]],axis=0)

# Use ggplot to map group-specific NBSS integral
plot = (ggplot(data) +
        facet_wrap('~PFT', ncol=1) +
        geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), color='black', fill='black') +
        geom_point(mapping=aes(x='longitude', y='latitude', size='Absolute_intercept',shape='Absolute_intercept.isna()',fill='Absolute_intercept.isna()'), alpha=1,color='black') +
        labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)', fill='', title='') +
        scale_fill_manual(values={True: "black", False: '#78214458'})+
        scale_size_area(trans='log10',max_size=3, na_value=1) + scale_shape_manual(limits=[True, False], values={True: "x", False: "o"},guide=None) +
        theme_paper +
        coord_cartesian(xlim=[-180, 180], ylim=[-90, 90], expand=False)).draw(show=True)
plot.set_size_inches(6,10)
plot.savefig(fname='{}/figures/PSSdb_paper/NBSS_maps.svg'.format(str(path_to_git)),dpi=300, bbox_inches='tight')

# Check NBSS parameter distributions and correlations
data_melt=data[(data.Max_abundance>0)].melt(id_vars=['Instrument','latitude','longitude','year','month','PFT','Group_index'],value_vars=['Min_size','Max_size','Slope','Absolute_intercept'],var_name='Parameter',value_name='Value')
data_melt['Parameter']=pd.Categorical(data_melt.Parameter,ordered=True,categories=['Min_size','Max_size','Slope','Absolute_intercept'])
data_melt['PFT']=pd.Categorical(data_melt.PFT,ordered=True,categories=['Mesophytoplankton','Rhizaria','Crustaceans'][::-1])
plot= (ggplot(data_melt[((data_melt.Parameter=='Slope')) ].dropna()) +
         facet_grid('Parameter~PFT',scales='free_x') +
         geom_violin(mapping=aes(y='-1*Value', group='PFT', x='PFT'), alpha=0.6,draw_quantiles=[0.05,0.5,0.95],width=0.1) +
         stat_summary(mapping=aes(y='-1*Value', group='PFT', x='PFT',label='stat(y)'),fun_y=np.median,size=4,fun_data='mean_cl_normal',geom='text',angle=-90, format_string='{:.1f}')+
         scale_y_sqrt(limits=[4,0])+scale_y_reverse(labels=lambda l: [-1*v for v in l])+
         labs(y=r'Slope (m$^{-3}$ $\mu$m$^{-3}$)', x=r'',size='') +
         theme_paper ).draw(show=True)
plot.set_size_inches(5,1.5)
plot.savefig(fname='{}/figures/PSSdb_paper/Distribution_slope.svg'.format(str(path_to_git)), dpi=300, bbox_inches='tight')
plot= (ggplot(data_melt[((data_melt.Parameter=='Absolute_intercept')) ].dropna()) +
         facet_grid('Parameter~PFT',scales='free_x') +
         geom_violin(mapping=aes(y='1*Value', group='PFT', x='PFT'), alpha=0.6,draw_quantiles=[0.05,0.5,0.95],width=0.1) +
         stat_summary(mapping=aes(y='1*Value', group='PFT', x='PFT',label='np.round(10**stat(y))'),fun_y=np.median,size=4,fun_data='mean_cl_normal',geom='text',angle=-90)+
         scale_y_log10()+
         labs(y=r'Intercept (particles m$^{-3}$)', x=r'',size='') +
         theme_paper ).draw(show=True)
plot.set_size_inches(5,1.5)
plot.savefig(fname='{}/figures/PSSdb_paper/Distribution_intercept.svg'.format(str(path_to_git)), dpi=300, bbox_inches='tight')
plot= (ggplot(data_melt[((data_melt.Parameter.isin(['Min_size','Max_size']))) ].dropna()) +
         facet_grid('Parameter~PFT',scales='free_x') +
         geom_violin(mapping=aes(y='Value', group='PFT', x='PFT'), alpha=0.6,draw_quantiles=[0.05,0.5,0.95],width=0.1,trim=True,adjust=0.5) +
         stat_summary(mapping=aes(y='Value', group='PFT', x='PFT',label='np.round(10**stat(y))'),fun_y=np.median,size=4,fun_data='mean_cl_normal',geom='text',angle=-90)+
         scale_y_log10()+
         labs(y=r'Size range ($\mu$m)', x=r'',size='') +
         theme_paper ).draw(show=True)
plot.set_size_inches(5,3)
plot.savefig(fname='{}/figures/PSSdb_paper/Distribution_size.svg'.format(str(path_to_git)),dpi=300, bbox_inches='tight')

data_transformed=data[data.Slope!=0].dropna(subset=['Intercept']).assign(transformed_slope=lambda x:np.log10(np.abs(x.Slope)),transformed_min_size=lambda x:np.log10(np.abs(x.Min_size)),transformed_max_size=lambda x:np.log10(np.abs(x.Max_size))).assign(untransformed_slope=lambda x:10**(x.transformed_slope),untransformed_min_size=lambda x:10**x.transformed_min_size,untransformed_max_size=lambda x:10**x.transformed_max_size)
data_transformed=data_transformed.mask(10**data_transformed.transformed_min_size==np.inf)
plot_data = pd.concat([ (data_transformed[['Group_index','PFT',col1, col2]].rename(columns={col1: 'values1', col2: 'values2'}).assign(col1=col1, col2=col2))   for col1, col2 in itertools.permutations(['transformed_slope','Intercept','transformed_min_size','transformed_max_size'], 2)])
plot_data['col1']=pd.Categorical(plot_data.col1,categories=['transformed_slope','Intercept','transformed_min_size','transformed_max_size'],ordered=True)
for var,lab in dict(zip(plot_data.col1.cat.categories,[r'Slope (m$^{-3}$ mm$^{-3}$)',r'Intercept (mm$^{3}$ m$^{-3}$ mm$^{-3}$)',r'Min size ($\mu$m)',r'Max size ($\mu$m)'])).items():
    for taxa in plot_data['PFT'].unique():
        YLIM=10**plot_data.loc[(plot_data.col1==var),'values1'].describe().loc[['min','max']].values if var!='transformed_slope' else -1*10**plot_data.loc[(plot_data.col1==var),'values1'].describe().loc[['min','max']].values
        TRANS='log10' if var!='transformed_slope' else 'identity'
        y='10**values1' if  var!='transformed_slope' else '-10**values1'
        plot= (ggplot(plot_data[(plot_data.col1==var) & (plot_data['PFT']==taxa)]) +
         facet_wrap('~ col2',scales='free_x') + # Scales free does not work with facet grid
         geom_point(mapping=aes(y=y, x='10**values2'),fill='#{:02x}{:02x}{:02x}{:02x}'.format(0,0, 0 , 0)) + #,size='CHL'
         geom_point(mapping=aes(y=y, x='10**values2'), fill='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 0)) + #,size='CHL'
         scale_y_continuous(trans=TRANS,limits=YLIM)+
         scale_x_log10()+scale_size_area()+
         labs(y=lab, x=r'',size='') +
         theme_paper).draw(show=True)
        plot.set_size_inches(8.5,2)
        plot.savefig(fname='{}/figures/PSSdb_paper/Correlation_{}_vs_parameters_{}.svg'.format(str(path_to_git),var,taxa),  dpi=300, bbox_inches='tight')

# Save NBSS parameters for model input
path_to_input=path_to_git / 'data' / 'Model'
path_to_input.mkdir(exist_ok=True,parents=True)
data.to_csv('{}/Model_input.csv'.format(str(path_to_input)),index=False)