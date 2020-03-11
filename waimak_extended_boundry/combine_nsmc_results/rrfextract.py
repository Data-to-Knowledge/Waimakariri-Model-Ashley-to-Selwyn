# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:12:55 2017

Functions to extract results and pars from pest_hp rrf file.


@author: briochh
"""
#import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import OrderedDict
import re
#import reiextract as reiex

#%%

def reiextract(pst,rei):
    # read in parameter groups from pest control file
    with open(pst,'rb') as fp:
        for line in fp: # loop over lines in control file
            # print line
            # find the number of parameter groups
            if '* control data' in line: # NPARGP in contol data secton
                for i,line in enumerate(fp):
                    if i==1: # secont line of control data section
                        #print line
                        NOBS=int(line.split()[1]) # NPARGP is 2 element in line
                        NOBSGP=int(line.split()[4]) # NOBSGP is 5 element in line
                        NPARGP=int(line.split()[2]) # NPARGP is 3 element in line
                        #print NPARGP
                        break # break to outer loop
                continue
            if '* observation groups' in line:
                obsgroups=[]
                modvobs=OrderedDict()
                # read NOBSGP number of lines after
                for i,line in enumerate(fp):
                    if i < NOBSGP:
                        group=line.split()[0]
                        obsgroups.append(group)
                        if 'regul' not in group:
                            modvobs[group]={}
                        # print line
                    else:
                        break
                break
    with open(rei,'rb') as fp:
        obscount=0
        readnext=False
        for i,line in enumerate(fp):
            if readnext:
                group=line.split()[1]
                obsnam=line.split()[0]
                modvobs[group][obsnam]=[np.float(num) for num in line.split()[2:]]
                obscount+=1
            if i==1:
                itnum=int(re.findall(r'\d+',line)[0])
            elif "Name" in line and line.split()[0] == "Name":   
                readnext=True
            if obscount==NOBS:
                readnext=False
                break
    return modvobs

def extractrrf(rrffile):
    """
    Extracts realisation parameters and run results from pest_hp run results file (rrf)
    
    Returns: Pandas dataframes of run outputs and parameter values.
    
    rrffile : string, run results file (incl path if required)
    """
    initflag='* case dimensions'
    parnamflag='* parameter names'
    obsnamflag='* observation names' 
    parsetflag='* parameter set index'
    
    with open(rrffile,'rb') as fp:
        print('Reading rrf file\n{}\n'.format(rrffile))
        rrf=fp.read().splitlines()
#    if pstfile==None:
#        case=os.path.splitext(os.path.basename(rrffile))
#        commondir=os.path.dirname(rrffile)
#        pstfile=os.path.join(commondir,'{}.pst'.format(case))
    # find indicies in rrf relating to each realisation set
    parset_idx=[]
    step=100 # for infor printing control
    print('Scanning rrf for realisation information start points')
    for i, s in enumerate(rrf):
        if initflag in s:
            casedim_idx = i+1
        if parnamflag in s:
            parnam_idx = i+1
        if obsnamflag in s:
            obsnam_idx=i+1
        if parsetflag in s:
            parset_idx.append(i+1)


    # use indices to get par names and obs names and numbers of obs and pars  
    NPAR,NOBS=[int(e) for e in rrf[casedim_idx].split()]
    parnams=rrf[parnam_idx:parnam_idx+NPAR]
    obsnams=rrf[obsnam_idx:obsnam_idx+NOBS]
    
    # EXTRACT DATA IN RRF TO PANDAS DATFRAMES
    parvals_df=pd.DataFrame(index=parnams)
    outvals_df=pd.DataFrame(index=obsnams)
    parvals_ls=[] # empty list for parameter values
    outvals_ls=[] # empty list for output values
    parset_ls=[] # empty list for parset name
    print('Extracting data from rrf')
    for i,rrfidx in enumerate(parset_idx): # loop over results sets in rrf
        if np.mod(i,step)==0:
            if i+step>len(parset_idx):
                to=len(parset_idx)
            else:
                to=i+step
            print('Extracting realisation {}->{} of {}'.format(i,to,len(parset_idx)))
        parset=int(rrf[rrfidx]) # paramter set name
        #parfile=rrf[rrfidx+2].split()[-1] # par set filename
        # paramter values start 4 rows down
        parvals=pd.to_numeric(
                pd.Series(rrf[rrfidx+4:rrfidx+4+NPAR],
                          index=parnams))
        # model outputs start after params 
        outvals=pd.to_numeric(
                pd.Series(rrf[rrfidx+4+NPAR+1:rrfidx+4+NPAR+1+NOBS],
                          index=obsnams))
        if not outvals.le(-1.11e+35).all(): # if failre pest_hp writes a dummy number
            parvals_ls.append(parvals) # list append in loop is quicker
            outvals_ls.append(outvals) 
            parset_ls.append(parset) # currently not used - param set names not preserved - just sequential
        else:
            parvals_ls.append(parvals)
            outvals_ls.append(outvals * np.nan)
            parset_ls.append(parset) # currently not used - param set names not preserved - just sequential
    print('Compiling DataFrames')
    parvals_df=pd.concat(parvals_ls,axis=1) # cat together dataframe
    parvals_df.columns=parset_ls
    outvals_df=pd.concat(outvals_ls,axis=1)
    outvals_df.columns=parset_ls
    return outvals_df,parvals_df

def calcphires(outvals_df, obs):
    """ Caluculate residual and phi contribution from each observation and 
    each realisation.
    
    outvals_df: pandas Dataframe,  model outputs for all realisations. 
    Indicies should be observation names, columns are realisation parameter sets
    obs: observation data as dataframe from pest control file. 
    Indicies consitent with outvals_df,observation weights also required for calculating contribution to phi.
    
    """
    # calculate associated residuals and phi contributions
    #obs=extractobs(pstfile) # pull observation values from pst
    residual_df=outvals_df.sub(obs['obs'],axis='index')
    phicont_df=residual_df.mul(obs['wgt'],axis='index')**2
    return phicont_df,residual_df
    
def extractobs(pstfile, NOBS=None):
    """ pstfile = path to file location
    """    
    # READ .PST FOR OBS DATA
    obsdataflag='* observation data'
    controldataflag='* control data'
    with open(pstfile,'rb') as fp:
        pst=fp.read().splitlines()
    for i,s in enumerate(pst):
        if NOBS==None:
            if controldataflag in s:
                NOBS=pst[i+2].split()[1]
        if obsdataflag in s:
            obsdata_idx=i+1
            break
    #obsdata_idx=[i+1 for i, s in enumerate(pst) if obsdataflag in s][0]
    obs=[line.split() for line in pst[obsdata_idx:NOBS+obsdata_idx]]
    obs=pd.DataFrame.from_records(obs).set_index(0)
    obs[[1,2]]=obs[[1,2]].astype('float64')
    obs.index=obs.index.str.lower()
    obs.index.name='index'
    obs.columns=['obs', 'wgt', 'group']
    return obs



def extractphisummary(recfile):
    """
    Read realisation results summaries from MC rec file
    
    Returns a pandas DataFrame of the total and groupwise phi contributions for each realisation
    """
    #READ REC with group and total phi contributions
    psetflag='File = '
    totalphiflag='Sum of squared weighted residuals (ie phi)'
    headphiflag='Contribution to phi from observation group "head"'
    vertphiflag='Contribution to phi from observation group "vert"'
    coastphiflag='Contribution to phi from observation group "coast"'
    drnphiflag='Contribution to phi from observation group "drn"'
    sfxphiflag='Contribution to phi from observation group "sfx"'
    sfophiflag='Contribution to phi from observation group "sfo"'


    with open(recfile,'rb') as fp:
        rec=fp.read().splitlines()
        
    pset=[]
    totalphi=[]
    headphi=[]
    vertphi=[]
    coastphi=[]
    drnphi=[]
    sfxphi=[]
    sfophi=[]
    readphi=False
    for s in rec:
        if psetflag in s:
            if readphi:
                totalphi.append(np.nan)
                headphi.append(np.nan)
                vertphi.append(np.nan)
                coastphi.append(np.nan)
                drnphi.append(np.nan)
                sfxphi.append(np.nan)
                sfophi.append(np.nan)
            pset.append(int(filter(str.isdigit,s)))
            readphi=True
        if totalphiflag in s:
            totalphi.append(float(s.split('=')[-1]))
            readphi=False
        if headphiflag in s:
            headphi.append(float(s.split('=')[-1]))
            readphi=False
        if vertphiflag in s:
            vertphi.append(float(s.split('=')[-1]))
            readphi=False
        if coastphiflag in s:
            coastphi.append(float(s.split('=')[-1]))
            readphi=False
        if drnphiflag in s:
            drnphi.append(float(s.split('=')[-1]))
            readphi=False
        if sfxphiflag in s:
            sfxphi.append(float(s.split('=')[-1]))
            readphi=False
        if sfophiflag in s:
            sfophi.append(float(s.split('=')[-1]))
            readphi=False
            
    if readphi: # should catch last parset if it is a failure
        totalphi.append(np.nan)
        vertphi.append(np.nan)
        headphi.append(np.nan)
        coastphi.append(np.nan)
        drnphi.append(np.nan)
        sfxphi.append(np.nan)
        sfophi.append(np.nan)
    gphidict={'total':totalphi,'head':headphi,'vert':vertphi,
              'coast':coastphi,'drn':drnphi,'sfx':sfxphi,
              'sfo':sfophi}    
    gpphicont_df=pd.DataFrame.from_dict(gphidict,orient='index')        
    gpphicont_df.columns=gpphicont_df.columns+1
    return gpphicont_df


def extractoptphi(optrec):
    """
    Extract the total and groupwise phi for the optimised model.
    
    Returns: Pandas DataFrame of group contributions for just the optimiesd model (1 column)
    """
    totalphiflag='Sum of squared weighted residuals (ie phi)'
    headphiflag='Contribution to phi from observation group "head"'
    vertphiflag='Contribution to phi from observation group "vert"'
    coastphiflag='Contribution to phi from observation group "coast"'
    drnphiflag='Contribution to phi from observation group "drn"'
    sfxphiflag='Contribution to phi from observation group "sfx"'
    sfophiflag='Contribution to phi from observation group "sfo"'
    #READ OPTIMISED GROUP PHI CONTRIBUTIONS        
    totalphi=[]
    headphi=[]
    vertphi=[]
    coastphi=[]
    drnphi=[]
    sfxphi=[]
    sfophi=[]
    
    with open (optrec,'rb') as fp:
        rec=fp.read().splitlines()
    for s in rec:
        if totalphiflag in s:
            totalphi.append(float(s.split('=')[-1]))
        if headphiflag in s:
            headphi.append(float(s.split('=')[-1]))
        if vertphiflag in s:
            vertphi.append(float(s.split('=')[-1]))
        if coastphiflag in s:
            coastphi.append(float(s.split('=')[-1]))
        if drnphiflag in s:
            drnphi.append(float(s.split('=')[-1]))
        if sfxphiflag in s:
            sfxphi.append(float(s.split('=')[-1]))
        if sfophiflag in s:
            sfophi.append(float(s.split('=')[-1]))
            break
    optphidict={'total':totalphi,'head':headphi,'vert':vertphi,
              'coast':coastphi,'drn':drnphi,'sfx':sfxphi,
              'sfo':sfophi}    
    optphicont_df=pd.DataFrame.from_dict(optphidict,orient='index')
    optphicont_df.columns=['opt']
    return optphicont_df
#%%
def filterphi(optphicont_df,gpphicont_df,thresh=None,group='total'):
    """
    Filter dataframe on phi contribution and return a list of indicies
    """
    #CUT DOWN TO JUST MODELS WITH TOTAL PHI BELOW THRESHOLD
    acceptable_mul=1.1
    optphi_tot=optphicont_df.loc['total'].values[0]
    if thresh==None:
        thresh=acceptable_mul*optphi_tot
    acceptsets=gpphicont_df.columns[gpphicont_df.loc[group]<thresh]
    return acceptsets

    
def plotgroupphi(groupplot_df):
    fig,ax = plt.subplots(figsize=(12,8))
    groupplot_df.reset_index().plot(x=groupplot_df.reset_index().index+1,
                         y='opt',ax=ax,kind='line',label='opt',color='grey')
    groupplot_df.T[1:].boxplot(ax=ax)

    ax.set_yscale('log')
    plt.show()
        

#%% Plot boxes for each observation
def boxplot(df,logy=True,obsdata=None,optdata=None,ylabel=None,box_per_plot=5.0):
    numsub=int(np.ceil(df.shape[0]/float(box_per_plot)))
    per_plot=int(box_per_plot)
    df=df.set_index('index')
    fig,ax = plt.subplots(numsub,1,figsize=(12,8))
    if obsdata is not None:
        df['obs']=obsdata['obs']
    if optdata is not None:
        df['optphi']=optdata
    if numsub==1: ax = [ax]
    for i,subax in enumerate(ax):
        print(i)
        plotdf=df[i*per_plot:(i+1)*per_plot]
        if obsdata is not None:
            #plotdf['obs']=obsdata[i*per_plot:(i+1)*per_plot]['obs']
            plotdf.reset_index().plot(x=plotdf.reset_index().index+1,
                               y='obs',kind='line',style='*',color='y',
                          ax=subax,label='obs',zorder=5,use_index=True)
            #subax.legend()
        if optdata is not None:
            #plotdf['optphi']=optdata[i*per_plot:(i+1)*per_plot]
            plotdf.reset_index().plot.line(x=plotdf.reset_index().index+1,
                               y='optphi',color='grey',
                          ax=subax,label='cal',zorder=5,use_index=True)
        plotdf.drop(['obs','optphi'],axis=1,errors='ignore').T.plot(kind='box',ax=subax)
        xticklabels=subax.get_xticklabels()
        #xlims=subax.get_xlim()
        #subax.set_xlim(xlims)
        subax.set_xticklabels(xticklabels,y=0.15,va='bottom',ha='right',x=0.1,color='grey',rotation='vertical')
        if logy:
            subax.set_yscale('log')
        subax.grid('on')
        subax.set_xlabel('')
        if ylabel is not None:
            subax.set_ylabel(ylabel)
    fig.subplots_adjust(hspace=0.001)
    fig.tight_layout(h_pad=0.05)
    return fig,ax
    
if __name__ == '__main__':
    obs, param = extractrrf(r"C:\Users\MattH\Desktop\nsmc_stuffs\nsmc\aw_ex_mc\aw_ex_mc.rrf")
    print 'done'
