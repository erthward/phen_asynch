"""
produce a global mp4 of annual,
minmax-scaled, fitted LSP patterns
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
from palettable.cmocean.sequential import Speed_20
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import imageio
import os
import sys
import re
import cv2


sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf


# data dir for saving intermediate images and final video
data_dir = '/media/deth/SLAB/diss/3-phn/phen_video'

# set crs for mapping
crs = 8857

# read in coefficients map 
coeffs = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                        '%s%s_coeffs.tif') % ('NIRv', ''))
# NOTE: swap coeffs axis from first to last
# (for numpy broadcast mult against individual day's design matrix vals)
coeffs = coeffs.transpose("y", "x", "band")

# create the harmonic regression design matrix
# (for recreating the fitted LSP curves)
dm = phf.make_design_matrix()


def plot_calendar_bar(ax, doy):
    '''
    plot a calendar timebar for a given day of the year
    NOTE: code adapted from:
        https://stackoverflow.com/questions/32111705/overlay-a-graph-over-a-video
    '''
    # plot the bar
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xrng = xmax - xmin
    yrng = ymax - ymin
    # plot timeline
    ax.plot([xmin+0.125*xrng, xmin+0.875*xrng],
            [ymin-0.045*yrng]*2,
            '-k',
            linewidth=3,
            clip_on=False,
           )
    # plot the months
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    tick_locs = np.linspace(xmin+0.125*xrng,
                            xmin+0.875*xrng,
                            len(months)+1,
                           )
    lab_locs = np.linspace(xmin+0.121*xrng,
                           xmin+0.871*xrng,
                           len(months)+1,
                          )
    for month, tick_loc, lab_loc in zip(months+[months[0]], tick_locs, lab_locs):
        ax.plot([tick_loc]*2,
                [ymin-0.06*yrng,
                 ymin-0.03*yrng],
                '-k',
                linewidth=3,
                clip_on=False,
               )
        ax.text(lab_loc,
                ymin-0.018*yrng,
                month,
                color='black',
                fontsize=5,
                clip_on=False,
               )
    # plot the current timepoint
    timepoint_frac = doy/365
    loc = xmin + 0.125*xrng + ((0.875-0.125)*timepoint_frac)*xrng
    ax.plot([loc]*2,
            [ymin-0.06*yrng, ymin-0.03*yrng],
            '-',
            color='#6d00a3',
            linewidth=8,
            alpha=0.5,
            clip_on=False,
           )


# get max and min fitted value maps
# (for minmax scaling of individual days' images)
if (not os.path.isfile('LSP_fitted_mins.txt') or
    not os.path.isfile('LSP_fitted_maxs.txt')):
    mins = np.ones(coeffs[:, :, 0].shape)*np.nan
    maxs = np.ones(coeffs[:, :, 0].shape)*np.nan
    I, J = np.where(pd.notnull(coeffs[:, :, 0]))
    for i, j in zip(I, J):
        # get the cell's fitted coefficients
        print(f"\n\tprocessing cell [{i}, {j}]...")
        coeffs_vec = coeffs[i, j, :].values
        # calc ts
        ts = phf.calc_time_series(coeffs_vec, dm)
        # get the min and max values
        mins[i, j] = np.min(ts)
        maxs[i, j] = np.max(ts)
    np.savetxt('LSP_fitted_mins.txt', mins)
    np.savetxt('LSP_fitted_maxs.txt', maxs)
else:
    mins = np.loadtxt('LSP_fitted_mins.txt')
    maxs = np.loadtxt('LSP_fitted_maxs.txt')

# create and save an image for each day of the year
cmap = Speed_20.mpl_colormap
for doy in range(365):
    fn = f'map_img_doy{doy}.png'
    if not os.path.isfile(os.path.join(data_dir, fn)):
        print(f'\n\tmapping day number {doy}...\n\n')
        # get the fitted image daily image
        fit = np.sum(coeffs * dm[doy, :], axis=2)
        # minmax-scale it
        fit_scaled = (fit - mins)/(maxs - mins)
        # reproject
        # (and mask out erroneous huge values after reprojection)
        fit_proj = fit_scaled.rio.reproject(crs)
        fit_proj = fit_proj.where(fit_proj<1e4, np.nan)
        assert np.nanmin(fit_proj) >= 0
        assert np.nanmax(fit_proj) <= 1
        # mask values outside global Equal Area projection bounds
        fit_proj = mask_xarr_to_other_xarr_bbox(fit_proj, fit)
        # plot it
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
        phf.plot_juris_bounds(lev1=False,
                              lev0_color='black',
                              lev0_alpha=1,
                              lev0_zorder=0,
                              strip_axes=False,
                             )
        fit_proj.plot.imshow(ax=ax,
                            vmin=0,
                            vmax=1,
                            cmap=cmap,
                            add_colorbar=False,
                            zorder=1,
                           )
        phf.plot_juris_bounds(lev1_linecolor='gray',
                              lev1_linewidth=0.1,
                              lev1_alpha=0.8,
                              lev1_zorder=2,
                              lev0_linecolor='gray',
                              lev0_linewidth=0.2,
                              lev0_alpha=0.8,
                              lev0_zorder=3,
                              strip_axes=True,
                             )
        fig.subplots_adjust(left=0,
                            right=1,
                            bottom=0,
                            top=1,
                           )
        # add the timeline
        plot_calendar_bar(ax, doy)
        fig.savefig(os.path.join(data_dir, f'map_img_doy{doy}.png'),
                    dpi=700,
                    pad_inches=0,
                    orientation='landscape',
                   )
        plt.close('all')
    else:
        print(f'\nday number {doy} already mapped.\n')

# turn into movie
print('\n\ncompiling video file...\n\n')
avifile = os.path.join(data_dir, 'VID_SUPP_normalized_NIRv_LSP.avi')
frame_filename_patt = 'map_img_doy\d{1,3}\.png'
pngs = [f for f in os.listdir(data_dir) if re.search(frame_filename_patt, f)]
# sort in day order
pngs.sort(key=lambda f: int(re.sub('\D', '', f)))
# load 0th image, to get dims
RGB_img = imageio.imread(os.path.join(data_dir, pngs[0]))
height, width, layers = RGB_img.shape
# about one month per second
frame_rate = int(365/12)
video = cv2.VideoWriter(avifile, 0, frame_rate, (width,height))
for png in pngs:
    video.write(cv2.imread(os.path.join(data_dir, png)))
cv2.destroyAllWindows()
video.release()

# compress the full-size video (otherwise it's ~32GB!)
cmd = "ffmpeg -i VID_SUPP_normalized_NIRv_LSP.avi -vcodec libx265 -crf 28 VID_SUPP_normalized_NIRv_LSP_COMPRESSED.mp4"
cwd = os.getcwd()
os.chdir(data_dir)
os.system(cmd)
os.chdir(cwd)

# knit into GIF
#print('\n\ncompiling GIF...\n\n')
#giffile = os.path.join(data_dir, "LSP_vid.gif")
#cmd = "convert -delay 10 -loop 0 $(ls -1 %s/map_img_doy*png | sort -V) %s" % (data_dir, giffile)
#os.system(cmd)

print('Yeehaw!')

