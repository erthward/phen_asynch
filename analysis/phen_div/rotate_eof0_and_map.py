import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import geopandas as gpd
import palettable
import xycmap
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.colors import LinearSegmentedColormap
from colorsys import hls_to_rgb



# TODO:
    # produce scree plot for EOF?
    # need to look into rotated EOFs as well?
    # need to center data (i.e., calc 'anomaly') prior to EOF, to make sure means = 0?
    # THERE'S STILL AN ARTEFACT IN EQUATORIAL REGION

# load the EOFs
eofs = rxr.open_rasterio('../../../results/maps/global_4_EOFs_coswts.tif')

# rescale each layer 0-1
for i in range(eofs.shape[0]):
    eofs[i] = (eofs[i]-eofs[i].min())/(eofs[i].max()-eofs[i].min())

# create a vertical vector that's 0 in the far north, 1 in the far south,
# and progresses from 0 to 1 as a sigmoid function across the tropics
# back, by cosine, within tropics
minmax_scale = lambda vals: (vals-np.min(vals))/(np.max(vals)-np.min(vals))
ys = eofs.y.values
wts = eofs.y*0
wts[ys<0] = 1
lats_in_tropics = eofs.sel(y=slice(23.4934, -23.4934)).y.values
# create sigmoid weighting
# NOTE: MINMAX VAL CHOSEN HEURISTICALLY, TO MINIMIZE NOTICEABLE COLOR-WARPING
#       ARTEFACTS IN EQUATORIAL REGION 
minmax_val = 10
wts_in_tropics = minmax_scale(1/(1 + np.exp(-np.linspace(-minmax_val,
                                                         minmax_val,
                                                         len(lats_in_tropics)))))
wts[(ys >= -23.4934) * (ys < 23.4934)] = wts_in_tropics

# use that to weighted-sum the regular eof0 array and the inverted eof0 array
eofs_wt_sum = deepcopy(eofs)
eofs_wt_sum[0] = (wts*(1-eofs[0])) + ((1-wts)*eofs[0])

#fig, axs = plt.subplots(1,2)
#eofs[:3].plot.imshow(ax=axs[0])
#axs[0].set_title('RGB EOFs: no transform')
#eofs_wt_sum[:3].plot.imshow(ax=axs[1])
#eofs_wt_sum[[1,0,2]].plot.imshow(ax=axs[1])
#eofs_wt_sum[[1,2,0]].plot.imshow(ax=axs[1])
#axs[1].set_title('RGB EOFs: tropical sigmoid-weighted inversion')
#fig.show()

# write to raster file
#eofs_wt_sum = eofs_wt_sum.rio.write_crs(4326)
#eofs_wt_sum.rio.to_raster('./transformed_scaled_EOFs.tif')



# create bivariate colormap within luminance-difference on each of the two
# axes proportional to the percent variation explained by the two axes' EOFs
pct_var_eof0 = 0.779
pct_var_eof1 = 0.115
# define colors
eof0_hue = 0.8
eof0_saturation = 0.7
eof1_hue = 0.3
eof1_saturation = 0.7
eof0_col0 = hls_to_rgb(eof0_hue, 0.2, eof0_saturation)
eof0_col1 = hls_to_rgb(eof0_hue, 0.9, eof0_saturation)
eof1_col0 = hls_to_rgb(eof1_hue, 0.2, eof1_saturation)
eof1_col1 = hls_to_rgb(eof1_hue, 0.2+((pct_var_eof1/pct_var_eof0)*0.7), eof1_saturation)
f, a = plt.subplots(1,1)
a.scatter([0,1,0,1], [0,0,1,1], c=[eof0_col0, eof0_col1, eof1_col0, eof1_col1])
# NOTE: chosen manually
ur_corner_col = '#a367c2'


# bivariate colormap
#fig2, ax = plt.subplots(1,1, figsize=(16,8))
#n = (50, 50)
#xcolors = ['#fae034', '#dcfa34', '#abfa34', '#52fa34', '#34fa9a', '#34fae0',
#           '#34cffa', '#3480fa', '#3445fa', '#6f34fa']
#xcmap = LinearSegmentedColormap.from_list('custom', xcolors)
#ycolors = ['#000000', '#ffffff']
#ycmap = LinearSegmentedColormap.from_list('custom', ycolors)
#cmap = xycmap.mean_xycmap(xcmap=xcmap, ycmap=ycmap, n=n)

# either pretty yellow, blue, and maroon colormap...
corner_colors = ("#d1d1d1", "#f5e102", "#0098d9", "#87054d")

# or a data-driven colormap based on the code above...
#corner_colors = (eof0_col0, eof0_col1, eof1_col1, ur_corner_col)

# or a data-driven colormap hand-chosen with Google's hex color picker
# (with logic displayed in HSL values and arithmetic on luminosity value,
#  in comments to right)
#corner_colors = ("#014a46", # HSL = 177deg, 97%, 15%
#                 "#9df1fa", # HSL = 186deg, 90%, 80%
#                 "#7d7102", # HSL = 54deg, 96%, 25% (25% = 15% + (.115/.779)*65%)
#                 "#9ff9c4", # HSL = 145deg, 88%, 80%,
#                )
cmap = xycmap.custom_xycmap(corner_colors=corner_colors, n=n)

eof0_vals = eofs_wt_sum[0].values.ravel()
eof1_vals = eofs_wt_sum[1].values.ravel()
non_nans = np.where(np.invert(np.isnan(eof0_vals)) *
                    np.invert(np.isnan(eof1_vals)))[0]
sx=eof0_vals[np.invert(np.isnan(eof0_vals))]
sy=eof1_vals[np.invert(np.isnan(eof1_vals))]
colors = xycmap.bivariate_color(sx=sx, sy=sy, cmap=cmap)
colors_arr = np.ones([len(eof0_vals), 3]) * np.nan
colors_arr[non_nans] = np.array([c[:3] for c in colors])
colors_arr = colors_arr.reshape([*eofs_wt_sum[0].shape, 3])
colors_xr = deepcopy(eofs_wt_sum[:3])
colors_xr[:,:,:] = np.swapaxes(np.swapaxes(colors_arr, 1, 2), 0, 1)

countries = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
countries.to_crs(4326).plot(color='none',
                            linewidth=0.5,
                            edgecolor='#4a4a4a',
                            ax=ax,
                            zorder=0,
                           )
colors_xr.plot.imshow(ax=ax)
ax.set_xticks(())
ax.set_yticks(())
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('')

# TODO: HOW TO ROTATE BIVAR LEGEND 45 DEG?
#plot_extents = -160, -110, -55, -20
#transform = Affine2D().rotate_deg(45)
#helper = floating_axes.GridHelperCurveLinear(transform, plot_extents)
#cax = floating_axes.FloatingSubplot(fig2, 111, grid_helper=helper)

n_lgd = (8,8)
cax = fig2.add_axes([0.07, 0.07, 0.12, 0.24])
lgd_vals = np.meshgrid(*[np.linspace(0, 1, i) for i in n_lgd])
lgd_cols = xycmap.bivariate_color(*[grid.ravel() for grid in lgd_vals], cmap=cmap)
lgd_cols_arr = np.array([c[:3] for c in lgd_cols]).reshape([*n_lgd, 3])
cax.imshow(lgd_cols_arr)
cax.invert_yaxis()
cax.set_xticks(np.linspace(0, n_lgd[0]-1, 5), np.round(np.linspace(0, 1, 5), 1))
cax.set_yticks(np.linspace(0, n_lgd[1]-1, 5), np.round(np.linspace(0, 1, 5), 1))
cax.set_xlabel('EOF 1', fontdict={'fontsize': 16})
cax.set_ylabel('EOF 2', fontdict={'fontsize': 16})
cax.tick_params(labelsize=9)
fig2.subplots_adjust(left=0, bottom=0, right=1, top=1)
fig2.show()

fig2.savefig('phen_bivar_plot.png', dpi=700)
