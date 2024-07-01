import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from shapely.geometry import Point
from shapely.geometry import Polygon as shapelyPolygon
from matplotlib.collections import PatchCollection
import numpy as np
import xarray as xr
import rasterio as rio
import rioxarray as rxr
from datetime import datetime
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns
import re, os, sys

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                                                            'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf



# plotting params
axlabel_fontdict={'fontsize': 14}
ticklabel_fontsize=10

# main behavioral params
seed = 1
np.random.seed(seed)

# make fig
fig = plt.figure(figsize=(13,9))
gs = fig.add_gridspec(nrows=90,ncols=100)

# load and plot flux-tower evaluation results
for i, rs_var in enumerate(['NIRv', 'SIF']):
    if i == 0:
        biome_start_col = 0
        biome_end_col = 32
        scat_start_col = 0
        scat_end_col = 32
    else:
        biome_start_col = 58
        biome_end_col = 90
        scat_start_col = 58
        scat_end_col = 90
    ax1 = fig.add_subplot(gs[:40, biome_start_col:biome_end_col])
    results_df = pd.read_csv('./FLUXNET_evaluation_results_%s.csv' % rs_var)

    whittaker = pd.read_csv('./whittaker_biomes.csv', sep=';')
    whittaker['temp_c'] = whittaker['temp_c'].apply(lambda x:
                                                float(x.replace(',', '.')))
    whittaker['precp_mm'] = whittaker['precp_cm'].apply(lambda x:
                                                float(x.replace(',', '.'))*10)


    def add_label_newline(label):
        if ' ' in label and label not in ['Woodland & shrubland',
                                          'Temperate grassland & desert']:
            label = '\n'.join(label.split(' '))
        return label


    biomes = []
    biome_labels = []
    centroids = []
    polygons = []
    patches = []
    for biome in whittaker['biome'].unique():
        subwhit = whittaker[whittaker.biome == biome].loc[:, ['temp_c', 'precp_mm']].values
        centroids.append(np.mean(subwhit, axis=0))
        shapely_poly = shapelyPolygon(subwhit)
        polygons.append(shapely_poly)
        poly = Polygon(np.array(shapely_poly.buffer(-0.12).exterior.coords.xy).T, True)
        patches.append(poly)
        biomes.append(biome)
        biome_labels.append(biome.replace('/', ' & '))

    # make whittaker biomes GeoDataFrame
    gdf = gpd.GeoDataFrame(pd.DataFrame({'biome': biomes,
                                         'geometry': polygons,
                                        }), geometry='geometry', crs=4236)

    biome_cmap = plt.get_cmap('tab20b', 17)
    # choose 9 colors that perform best for colorblindness
    color_indices = [0, 2, 3, 12, 13, 15, 4, 5, 7]
    colors = []
    hex_colors = []
    for i in color_indices:
        rgba = biome_cmap(i)
        hex_color = mpl.colors.rgb2hex(rgba)
        hex_colors.append(hex_color)

    # NOTE: colors manually pulled from the tab20b palette display, and a little
    # from other palettes on the matplotlib reference page,
    # after cross-checking with a colorblindness simulator
    hex_colors = [
                  '#393B79', # dark blue
                  '#9C9EDE', # light blue
                  '#FFD92F', # yellow
                  '#E7969C', # light pink
                  '#B5CF6B', # light green
                  '#637939', # drab olive green
                  '#31A354', # jungle green
                  '#D6616B', # rose-magenta
                  '#843C39', # dark magenta
                 ]
    rgba_colors = [mpl.colors.hex2color(h) for h in hex_colors]

    #if rs_var == 'SIF':
    if False:
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.1)
    p = PatchCollection(patches, alpha=1, color='white')
    p.set_edgecolor(rgba_colors)
    p.set_linewidth(1.5)
    ax1.add_collection(p)
    scat = ax1.scatter(results_df['mat'],
                       results_df['map'],
                       c = results_df['r2'],
                       vmin=0,
                       vmax=1,
                       edgecolor='black',
                       linewidths=0.6,
                       s=10,
                       alpha=1,
                       cmap='gray')
    ax1.set_xlabel('MAT ($^{\circ}C$)',
                  fontdict=axlabel_fontdict)
    ax1.set_ylabel('MAP ($mm$)', fontdict=axlabel_fontdict)
    ax1.tick_params(labelsize=ticklabel_fontsize)
    if rs_var == 'NIRv':
        title = '$NIR_{V}$'
    else:
        title = '$SIF$'
    ax1.set_title(title, fontdict={'fontsize': 20})
    # add colorbar at far right for the R^2 scatter colors
    if rs_var == 'SIF':
        cbar_ax = fig.add_subplot(gs[:40, 94:97])
        cbar = mpl.colorbar.ColorbarBase(cbar_ax,
                                         cmap='gray',
                                         orientation='vertical',
                                        )
        cbar_ax.set_ylabel('$GPP_{FLUXNET}\sim LSP\ R^2$', fontdict={'fontsize': 11})

    ax2 = fig.add_subplot(gs[50:, scat_start_col:scat_end_col])
    pt_biomes = []
    for i, row in results_df.iterrows():
        if pd.notnull(row['mat']) and pd.notnull(row['map']):
            pt = Point(row['mat'], row['map'])
            pt_biome = gdf[[pt.within(geom) for geom in gdf.geometry]].biome.values
            try:
                assert len(pt_biome) == 1
                pt_biomes.append(pt_biome[0])
            except Exception as e:
                pt_biomes.append(np.nan)
        else:
            pt_biomes.append(np.nan)

    # tweak formatting to thin the legend
    new_pt_biomes = []
    for pb in pt_biomes:
        if pd.notnull(pb):
            space_idxs = []
            #for i, char in enumerate(pb):
            #    if char == ' ':
            #        space_idxs.append(i)
            #idx_closest_to_20 = space_idxs[np.argmin(np.abs(
            #                                            np.array(space_idxs)-20))]
            #new_pb = pb[:idx_closest_to_20] + '\n' + pb[idx_closest_to_20+1:]
            new_pb = re.sub('\s', '\n', pb)
            new_pb = re.sub('/', ';\n', new_pb)
            new_pt_biomes.append(new_pb)
        else:
            new_pt_biomes.append(pb)

    results_df['biome'] = new_pt_biomes

    # do same thing to biome-order list
    new_biomes = []
    for b in biomes:
        if pd.notnull(b):
            space_idxs = []
            #for i, char in enumerate(b):
            #    if char == ' ':
            #        space_idxs.append(i)
            #idx_closest_to_20 = space_idxs[np.argmin(np.abs(
            #                                            np.array(space_idxs)-20))]
            #new_b = b[:idx_closest_to_20] + '\n' + b[idx_closest_to_20+1:]
            new_b = re.sub('\s', '\n', b)
            new_b = re.sub('/', ';\n', new_b)
            new_biomes.append(new_b)
        else:
            new_biomes.append(b)

    # add dotted vertical lines for each year,
    # and a solid line indicating length of the remote sensing ts, for comparison
    for yrs in np.arange(1, np.max(results_df['gpp_ts_len'])+1, 1):
        ax2.axvline(x=yrs, ymin=-1, ymax=2, linewidth=0.5, alpha=0.6,
                    color='black', linestyle=':', zorder=0)
    if rs_var == 'NIRv':
        rs_dataset_len = 10
    else:
        rs_dataset_len = 4.33
    ax2.axvline(x=rs_dataset_len, ymin=-1, ymax=2, linewidth=1.3, alpha=0.8,
                color='black', linestyle='-', zorder=0)

    sns.set_palette(sns.color_palette(hex_colors))

    scat = sns.scatterplot(x=results_df['gpp_ts_len'] + np.random.uniform(-0.25,
                                                                          0.25,
                                                                len(results_df)),
                           y='r2',
                           data=results_df.sort_values('r2', ascending=False),
                           hue='biome',
                           hue_order=new_biomes,
                           #style='match_veg',
                           #markers=['X', 'o', '*'],
                           style_order=[0, 1, 2],
                           color='black',
                           edgecolor='none',
                           alpha=0.75,
                           s = 30,
                           legend=rs_var=='NIRv',
                           #ax=ax2,
                          )
    if rs_var == 'NIRv':
        leg = scat.legend(loc='upper center',
                    bbox_to_anchor=(1.32, 1.8),
                    ncol=1,
                    fontsize=11,
                    labelspacing=0.8,
                    markerscale=1.5,
                   )
    ax2.set_xlim(0, np.max(results_df['gpp_ts_len'])+1)
    ax2.set_ylim(0,1)
    ax2.set_xlabel('GPP time series length (years)', fontdict=axlabel_fontdict)
    if rs_var == 'NIRv':
        ylab_str = 'NIR_{V}'
    else:
        ylab_str = rs_var
    ax2.set_ylabel('$R^2$', fontdict=axlabel_fontdict)
    ax2.tick_params(labelsize=ticklabel_fontsize)

    # add part labels
    if rs_var == 'NIRv':
        ax1_lab = 'A.'
        ax2_lab = 'B.'
    else:
        ax1_lab = 'C.'
        ax2_lab = 'D.'
    ax1.text(ax1.get_xlim()[0]-(0.24*np.diff(ax1.get_xlim())),
             1.005*ax1.get_ylim()[1],
             ax1_lab,
             weight='bold',
             size=24,
            )
    ax2.text(-5.1,
             1,
             ax2_lab,
             weight='bold',
             size=24,
            )


fig.subplots_adjust(top=0.94, bottom=0.08, left=0.08, right=0.98)

fig.savefig(os.path.join(phf.FIGS_DIR, 'FIG_SUPP_flux_evaluation_results.png'), dpi=700)

