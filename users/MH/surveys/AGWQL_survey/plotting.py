"""
Author: matth
Date Created: 13/03/2017 9:59 AM
"""

from __future__ import division
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import rasterio as rast
import matplotlib.image as mpimg
import matplotlib.patches as patches
from users.MH.Well_classes import mav, gv
import os


# plot background
def plot_background(gis_dir):
    background = rast.open(os.path.join(gis_dir, "background.tif"))
    ll = background.transform * (0, background.height)
    ur = background.transform * (background.width, 0)

    fig, ax = plt.subplots(1, figsize=(10,12.25))
    image = mpimg.imread(os.path.join(gis_dir,"background.tif"))
    ax.imshow(image, extent=[ll[0], ur[0], ll[1], ur[1]])  # make this a function at somepoint
    lakes = gpd.read_file(os.path.join(gis_dir,"lakes.shp"))
    seds = gpd.read_file(os.path.join(gis_dir,"cen_seds.shp"))
    cwms = gpd.read_file(os.path.join(gis_dir,"cwms_zones.shp"))

    seds_plot = seds.plot(ax=ax, alpha=1, color='khaki', edgecolor='none')
    lakes_plot = lakes.plot(ax=ax, alpha=0.25, color='navy', edgecolor='none')
    cwms_plot = cwms.plot(ax=ax, facecolor='none', linewidth=1.5)
    # add CWMS labels
    labels = [(1584829.848, 5310914.210, 'Kaikoura'),
              (1629701.336, 5318952.893, 'Kaikoura'),
              (1526278.030, 5273534.192, 'Hurunui-Waiau'),
              (1525064.257, 5226076.800, 'Waimak-\nariri'),
              (1582640.995, 5187074.497, 'Christchurch-West Melton'),
              (1614379.364, 5155883.341, 'Banks\nPeninsula'),
              (1458970.779, 5227568.278, 'Selwyn-Waihora'),
              (1431610.136, 5191999.416, 'Ashburton'),
              (1414646.525, 5107728.575, 'Orari-Temuka-\nOpihi-Pareora'),
              (1374152.743, 5165733.180, 'Upper\nWaitaki'),
              (1457876.372, 5059026.594, 'Lower Waitaki-\nSouth Coastal\nCanterbury')
              ]
    for label in labels:
        ax.text(*label, size=8, horizontalalignment='left',
                verticalalignment='top', fontweight='bold')

    # add arrow
    length = 30000
    x=1.58032e6
    y=5.02824e6
    plt.arrow(x,y,0,length,linewidth=1.5, head_length=0.2*length,head_width=0.2*length,facecolor='k')
    plt.text(x,y+length/3,'N', size=20, fontweight='bold', ha='center', va='center')

    # add scale bar
    start = (1.46044e6, 5.00845e6)
    ax.add_patch(patches.Rectangle((start[0],start[1]),25000,5000,facecolor='k',linewidth=1))
    ax.add_patch(patches.Rectangle((start[0]+25000,start[1]),25000,5000,facecolor='w',linewidth=1))
    ax.add_patch(patches.Rectangle((start[0]+50000,start[1]),25000,5000,facecolor='k',linewidth=1))
    ax.add_patch(patches.Rectangle((start[0]+75000,start[1]),25000,5000,facecolor='w',linewidth=1))
    ax.add_patch(patches.Rectangle((start[0]+100000,start[1]),25000,5000,facecolor='k',linewidth=1))
    text = [0,25,50,75,100]
    for val in text:
        str_val = str(val)
        if val==100:
            str_val = str(val)+' km'
        plt.text(start[0]+val*1000,start[1]+10000,str_val,size=10)

    return fig, ax


if __name__ == '__main__':
    plotdir = "P:/Groundwater/Annual groundwater quality survey 2016/figures/"
    gis_dir = "P:/Groundwater/Annual groundwater quality survey 2016/GIS_data"
    # plotting for the AGWQL survey
    data = pd.read_csv("P:/Groundwater/Annual groundwater quality survey 2016/data/all_output.csv")
    data = data[pd.notnull(data['NO3_N_trend'])]

    # well depth
    fig, ax = plot_background(gis_dir)
    tempdata = data[data['depth'] <= 20]
    p1 = plt.plot(tempdata.lon, tempdata.lat, color='b', marker='o', label=u'\u2264 20 m', linestyle='None')

    tempdata = data[(data['depth'] > 20) & (data['depth'] <= 50)]
    plt.plot(tempdata.lon, tempdata.lat, color='g', marker='v', label='20-50 m', linestyle='None')

    tempdata = data[(data['depth'] > 50) & (data['depth'] <= 100)]
    plt.plot(tempdata.lon, tempdata.lat, color='darkorange', marker='s', label='50-100 m', linestyle='None')

    tempdata = data[data['depth'] > 100]
    plt.plot(tempdata.lon, tempdata.lat, color='r', marker='p', label='> 100 m', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    # add seds and zones to legend
    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])
    plt.legend(handels, labels, title='Depth', loc='lower right')
    plt.savefig(plotdir+'depth.png',dpi=300, bbox_inches='tight')
    plt.close()

    # Nitrate N
    fig, ax = plot_background(gis_dir)

    tempdata = data[data['NO3_N'] <= 3]
    plt.plot(tempdata.lon, tempdata.lat, color='b', marker='o', label=u'\u2264 3 mg/L', linestyle='None')

    tempdata = data[(data['NO3_N'] > 3) & (data['NO3_N'] <= 5.6)]
    plt.plot(tempdata.lon, tempdata.lat, color='c', marker='o', label='3 to 5.6 mg/L', linestyle='None')

    tempdata = data[(data['NO3_N'] > 5.6) & (data['NO3_N'] <= 6.9)]
    plt.plot(tempdata.lon, tempdata.lat, color='orange', marker='o', label='5.6 to 6.9 mg/L', linestyle='None')

    tempdata = data[(data['NO3_N'] > 6.9) & (data['NO3_N'] <= 11.3)]
    plt.plot(tempdata.lon, tempdata.lat, color='lightcoral', marker='o', label='6.9 to 11.3 mg/L', linestyle='None')

    tempdata = data[data['NO3_N'] > 11.3]
    plt.plot(tempdata.lon, tempdata.lat, color='r', marker='o', label='Over 11.3 mg/L', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])

    plt.legend(handels, labels, title='Nitrate Nitrogen', loc='lower right')
    plt.savefig(plotdir+'NO3_N.png',dpi=300, bbox_inches='tight')
    plt.close()

    # Nitrate Trends
    fig, ax = plot_background(gis_dir)

    tempdata = data[data['NO3_N_trend'] == 'increasing']
    plt.plot(tempdata.lon, tempdata.lat, color='r', marker='^', label='Increasing Trend', linestyle='None')

    tempdata = data[data['NO3_N_trend'] == 'decreasing']
    plt.plot(tempdata.lon, tempdata.lat, color='b', marker='v', label='Decreasing Trend', linestyle='None')

    tempdata = data[data['NO3_N_trend'] == 'no trend']
    plt.plot(tempdata.lon, tempdata.lat, color='grey', marker='o', label='No Trend', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])
    plt.legend(handels, labels, title='Nitrate Nitrogen Trend', loc='lower right')
    plt.savefig(plotdir+'NO3_N_trends.png',dpi=300, bbox_inches='tight')
    plt.close()

    # DRP
    fig, ax = plot_background(gis_dir)

    tempdata = data[data['drp'] <= 0.004]
    plt.plot(tempdata.lon, tempdata.lat, color='b', marker='o', label='< 0.004 mg/L', linestyle='None')

    tempdata = data[(data['drp'] > 0.004) & (data['drp'] <= 0.009)]
    plt.plot(tempdata.lon, tempdata.lat, color='c', marker='o', label='0.004 to 0.009 mg/L', linestyle='None')

    tempdata = data[(data['drp'] > 0.009) & (data['drp'] <= 0.030)]
    plt.plot(tempdata.lon, tempdata.lat, color='orange', marker='o', label='0.009 to 0.030 mg/L', linestyle='None')

    tempdata = data[data['drp'] > 0.030]
    plt.plot(tempdata.lon, tempdata.lat, color='r', marker='o', label='Over 0.030 mg/L', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])

    plt.legend(handels, labels, title='Dissolved Reactive\nPhosphorous', loc='lower right')
    plt.savefig(plotdir+'drp.png',dpi=300, bbox_inches='tight')
    plt.close()

    # DRP trends
    fig, ax = plot_background(gis_dir)

    tempdata = data[data['drp_trend'] == 'increasing']
    plt.plot(tempdata.lon, tempdata.lat, color='r', marker='^', label='Increasing Trend', linestyle='None')

    tempdata = data[data['drp_trend'] == 'decreasing']
    plt.plot(tempdata.lon, tempdata.lat, color='b', marker='v', label='Decreasing Trend', linestyle='None')

    tempdata = data[data['drp_trend'] == 'no trend']
    plt.plot(tempdata.lon, tempdata.lat, color='grey', marker='o', label='No Trend', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])

    plt.legend(handels, labels, title='Dissolved Reactive\nPhosphorous Trend', loc='lower right')
    plt.savefig(plotdir+'drp_trend.png',dpi=300, bbox_inches='tight')
    plt.close()

    # Ecoli
    fig, ax = plot_background(gis_dir)

    tempdata = data[data['ecoli'] < 0]
    plt.plot(tempdata.lon, tempdata.lat, color='b', marker='o', label='Not Detected', linestyle='None')

    tempdata = data[data['ecoli'] > 0]
    plt.plot(tempdata.lon, tempdata.lat, color='r', marker='o', label='Detected', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])

    plt.legend(handels, labels, title='E. coli Detection', loc='lower right')
    plt.savefig(plotdir+'ecoli.png',dpi=300, bbox_inches='tight')
    plt.close()

    # Mn
    fig, ax = plot_background(gis_dir)

    tempdata = data[data['Mn'] <= 0]
    plt.plot(tempdata.lon, tempdata.lat, color='b', marker='o', label='Not Detected', linestyle='None')

    tempdata = data[(data['Mn'] > 0) & (data['Mn']<=gv['Mn'])]
    plt.plot(tempdata.lon, tempdata.lat, color='c', marker='o', label=u'Mn \u2264 GV', linestyle='None')

    tempdata = data[data['Mn'] > gv['Mn']]
    plt.plot(tempdata.lon, tempdata.lat, color='orange', marker='o', label='Mn > GV', linestyle='None')

    tempdata = data[data['Mn'] > mav['Mn']]
    plt.plot(tempdata.lon, tempdata.lat, color='r', marker='o', label='Mn > MAV', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])

    plt.legend(handels, labels, loc='lower right')
    plt.savefig(plotdir+'Mn.png',dpi=300, bbox_inches='tight')
    plt.close()

    # Fe
    fig, ax = plot_background(gis_dir)

    tempdata = data[data['Fe'] <= 0]
    plt.plot(tempdata.lon, tempdata.lat, color='b', marker='o', label='Not Detected', linestyle='None')

    tempdata = data[(data['Fe'] > 0) & (data['Mn'] <= gv['Fe'])]
    plt.plot(tempdata.lon, tempdata.lat, color='c', marker='o', label=u'Fe \u2264 GV', linestyle='None')

    tempdata = data[data['Fe'] > gv['Fe']]
    plt.plot(tempdata.lon, tempdata.lat, color='orange', marker='o', label='Fe > GV', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])

    plt.legend(handels, labels, loc='lower right')
    plt.savefig(plotdir+'Fe.png',dpi=300, bbox_inches='tight')
    plt.close()

