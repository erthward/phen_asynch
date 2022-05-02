import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import geopandas as gpd
import palettable
import xycmap
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes


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
minmax_val = np.e
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



# bivariate colormap
fig2, ax = plt.subplots(1,1, figsize=(16,8))
n = (50, 50)
#xcmap = palettable.cartocolors.sequential.Peach_3.mpl_colormap
#xcmap = palettable.colorbrewer.sequential.Oranges_8.mpl_colormap
#ycmap = palettable.colorbrewer.sequential.PuRd_8.mpl_colormap
#cmap = xycmap.mean_xycmap(xcmap=xcmap, ycmap=ycmap, n=n)
corner_colors = ("#d1d1d1", "#f5e102", "#0098d9", "#87054d")
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
