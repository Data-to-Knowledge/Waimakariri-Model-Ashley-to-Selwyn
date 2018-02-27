# -*- coding: utf-8 -*-
"""
Created on Tue Sep 05 13:37:29 2017

@author: briochh
"""

"""
Created on Mon Jul 17 14:45:46 2017

@author: briochh
"""
#%%
#import mf
import os
#import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import reiextract as reiex
from users.MH.Waimak_modeling.models.extended_boundry.nsmc_exploration_results.combine_nsmc_results import \
    rrfextract as rrf
# %%
# import mf
import os

# import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import reiextract as reiex
from users.MH.Waimak_modeling.models.extended_boundry.nsmc_exploration_results.combine_nsmc_results import \
    rrfextract as rrf

#from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.close('all')
save=False
mfdir=r'I:\Groundwater\commercial\ECan\Waimak\Models\AshWaiExtend\results\AWMC20171026_postopt' # directory to MC results #######!!!!!
prefix='aw_ex_mc'  # prefix of MC files

rrffile=os.path.join(mfdir,'{}.rrf'.format(prefix))  # rrf
pstfile=os.path.join(mfdir,'{}.pst'.format(prefix))  # pest conrtol file
recfile=os.path.join(mfdir,'{}.rec'.format(prefix))  # group phi contribution results for each realisation

optmoddir=r'I:\Groundwater\commercial\ECan\Waimak\Models\AshWaiExtend\results\AWJCO20171024' # optimised model directory (it is not run as part of MC) #####!!!!!
optmodpref='aw_ex_jco' # optimised model prefix
optrei=os.path.join(optmoddir,'{}.rei'.format(optmodpref)) # optimised model misfit file (all observations, results and residuals)
optpst=os.path.join(optmoddir,'{}.pst'.format(optmodpref)) # optimised model pest control file
optrec=os.path.join(optmoddir,'{}.rec'.format(optmodpref)) # optimised model group phi contributions 
#%% read residuals from optimised model
modvobs=rrf.reiextract(optpst,optrei) # read the optimised model rei. returns an ordered dictionary of dictionaries
# keys are observation groups, {obsgroup:{obsname:['obs','modelled','residual','wgt']}} # DataFrame would be nicer!!!!!


#if save:
#    if not os.path.exists(os.path.join(mfdir,'plots')):
#        os.mkdir(os.path.join(mfdir,'plots'))


#%% # GET INDICIES of each observations (makes plotting marginally quiker as we can use the .iloc pandas access)
count={}
runningcount=0
obsgroups={}
for group,values in modvobs.iteritems():
    count[group]=len(values)
    obsgroups[group]=range(runningcount,runningcount+count[group])
    runningcount+=count[group]


#%% 
# extract realisation outpouts and parameter values from pest_hp MC rrf  
outvals_df,parvals_df=rrf.extractrrf(rrffile) 
# extract actual observation values
obs=rrf.extractobs(pstfile, NOBS=1409) # there are alternative to this (e.g. already in modvobs dict)
# calculate the residual and the phi contribution    
phicont_df,residual_df=rrf.calcphires(outvals_df, obs)
# extract the group hpi contributions from the MC rec file
gpphicont_df=rrf.extractphisummary(recfile)
# extract the groupwise phi contributions from the optimised model for comparison
optphicont_df=rrf.extractoptphi(optrec)


#%% FILTER 1
# Filter on total phi meeting a threshold.
acceptsets=rrf.filterphi(optphicont_df,gpphicont_df,thresh=15484) # this threshold of acceptance was the phi of the other optimiesd model (i1)

# refine our DataFrame of models to just those the can be considered acceptable based on filter 1 
phicont_df=phicont_df.loc[:,acceptsets]
outvals_df=outvals_df.loc[:,acceptsets]
residual_df=residual_df.loc[:,acceptsets]
# total and groupwise phi
gpphicont_df=gpphicont_df.loc[:,acceptsets]

# reset index means we can use quicker iloc  
outvals_df=outvals_df.reset_index(drop=False)
residual_df=residual_df.reset_index(drop=False)
phicont_df=phicont_df.reset_index(drop=False)

# merge calib and results for plot
groupplot_df=optphicont_df.merge(gpphicont_df,left_index=True,right_index=True)
# plot boxplot of groupwise and total phi for all realiastions meeting threshold 1
rrf.plotgroupphi(groupplot_df)

#%%
# plot boxplots for each observation for realisations meeting Filter 1
print('plotting boxes')
for group in obsgroups.keys():
    nobs=len(outvals_df.iloc[obsgroups[group]])
    nperplot=20
    plotperfig=5
    boxperfig=nperplot*plotperfig
    nplots=(nobs/boxperfig)+1 # integer
    optres=pd.DataFrame.from_dict(modvobs[group],orient='index')
    optres.index.name='index'
    optres.columns=['obs','mod','res','wgt']
    optres['phi']=(optres.res*optres.wgt)**2
    optres=optres.reindex(phicont_df.iloc[obsgroups[group]].set_index('index').index)
    print('plotting {} output boxes'.format(group)) 
    for i,idx in enumerate(zip(range(0,nplots*boxperfig,boxperfig),\
                                   range(boxperfig,nobs,boxperfig)+[nobs])):    
        print('{} figure {}/{}'.format(group,i+1,nplots))
        fighout,axhout = rrf.boxplot(outvals_df.iloc[obsgroups[group]][idx[0]:idx[1]],
                             obsdata=obs.iloc[obsgroups[group]][idx[0]:idx[1]],
                             logy=False,
                             ylabel='{} model output'.format(group),box_per_plot=nperplot)    
        print('plotting {} phi contribution boxes'.format(group))  
        fighphi,axhphi = rrf.boxplot(phicont_df.iloc[obsgroups[group]][idx[0]:idx[1]],
                         optdata=optres['phi'][idx[0]:idx[1]],
                         ylabel='{} phi contrib.'.format(group),box_per_plot=nperplot)    


plt.show()
