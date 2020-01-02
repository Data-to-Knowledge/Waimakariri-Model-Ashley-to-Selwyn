"""
Author: matth
Date Created: 23/03/2017 2:20 PM
"""
# this could be made better to retun the UNC(//gisdata/Projects/SCI/) as default or mapped path (P://) if I ever have time.
# or could simple make a function to convert UNC to mapped for supporting the odd thing that needs it (e.g. DOS)

def gisdata(path):
    return '{}/{}'.format('//GISDataFS/GISData',path)

def sci(path):
    return '{}/{}'.format(r'D:\ecan_data',path)

def data(path):
    return '{}/{}'.format('//fileservices02/managedshares/data', path)

def transfers(path):
    return '{}/{}'.format(r'C:\matt_modelling_unbackedup',path)

def temp(path):
    return transfers(r"temp/{}".format(path))

def gw_met_data(path):
    return r'D:\ecan_data/{}'.format(path)

