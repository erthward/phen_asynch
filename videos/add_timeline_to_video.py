import cv2
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
import imageio
import os
import re

# from mp4 or from PNGs?
#from_fmt = 'mp4'
from_fmt = 'png'

# delete intermediate files?
delete_intermediate = False


# NOTE: code adapted from:
    # https://stackoverflow.com/questions/32111705/overlay-a-graph-over-a-video

# function to make timebar for a given day
def plot_calendar_bar(ax, frame_num, frame_ct):
    # plot the bar
    if from_fmt == 'mp4':
        p = Polygon(np.array([[100,380],[100,350], [700,350], [700,380], [100,380]]))
    else:
        p = Polygon(np.array([[260,3500],[260,3200], [6940,3200], [6940,3500],
                              [260,3500]]))
        p2 = Polygon(np.array([[360,3300],[360,3220], [6840,3220], [6840,3300],
                              [360,3300]]))
    pc = PatchCollection([p], alpha=1, edgecolor='gray', facecolors='#e5e3e6')
    ax.add_collection(pc)
    if from_fmt == 'png':
        pc2 = PatchCollection([p2], alpha=1, edgecolor='black',
                              facecolors='white')
        ax.add_collection(pc2)
    # plot the months
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    if from_fmt == 'mp4':
        locs = np.linspace(105, 685, len(months))
        for month, loc in zip(months, locs):
            ax.plot([loc+5]*2, [350, 365], '-k', linewidth=0.5)
            ax.text(loc-3, 377, month, color='black', fontsize=12,)
                    #fontdict={'fontweight': 'bold'})
    else:
        tick_locs = np.linspace(360, 6840, len(months)+1)
        lab_locs = np.linspace(315, 6795, len(months)+1)
        for month, tick_loc, lab_loc in zip(months+[months[0]], tick_locs, lab_locs):
            ax.plot([tick_loc]*2, [3225, 3300], '-k', linewidth=1.5)
            ax.text(lab_loc, 3430, month, color='black', fontsize=12,)
                    #fontdict={'fontweight': 'bold'})
    # plot the current timepoint
    timepoint_frac = frame_num/frame_ct
    if from_fmt == 'mp4':
        loc = 100 + ((700-100)*timepoint_frac)
        ax.plot([loc]*2, [376, 354], '-', color='#6d00a3', linewidth=9, alpha=0.5)
    else:
        loc = 360 + ((6840-360)*timepoint_frac)
        ax.plot([loc]*2, [3275, 3245], '-', color='#6d00a3', linewidth=9, alpha=0.5)


subplots_adj_left=0
subplots_adj_bottom=0
subplots_adj_right=1
subplots_adj_top=1
subplots_adj_wspace=0
subplots_adj_hspace=0

if from_fmt == 'mp4':
    data_dir = '/home/deth/Desktop/CAL/research/projects/seasonality/results/videos'
else:
    data_dir = '/media/deth/SLAB/seasonality/other/video_images/'

if from_fmt == 'mp4':
    for vidfile in ['globalNIRvseasonalityVideo.mp4',
                    'globalSIFseasonalityVideo.mp4']:
        print('\n\nNOW PROCESSING %s...\n\n' % vidfile)
        cap = cv2.VideoCapture(os.path.join(data_dir,vidfile))

        try:
            frame_ct = int(cap.get(cv2.cv.CV_CAP_PROP_FRAME_COUNT))
            width  = int(cap.get(cv2.cv.CV_CAP_PROP_FRAME_WIDTH))
            height = int(cap.get(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT))
        except AttributeError:
            frame_ct = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
            width  = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
            height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

        for frame_num in range(frame_ct):
            print('\n\tprocessing frame %i of %i...\n\n' % (frame_num, frame_ct))
            flag, frame = cap.read()
            # reorder from BGR (OpenCV's convention) to RGB
            RGB_img = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches(16,8)
            ax.imshow(RGB_img)
            plot_calendar_bar(ax, frame_num, frame_ct)
            plt.subplots_adjust(left=subplots_adj_left,
                                bottom=subplots_adj_bottom,
                                right=subplots_adj_right,
                                top=subplots_adj_top,
                                wspace=subplots_adj_wspace,
                                hspace=subplots_adj_hspace)
            fig.savefig(os.path.join(data_dir, 'out_frame_%i.png' % frame_num),
                        dpi=100,
                        pad_inches=0,
                        orientation='landscape')
            plt.close('all')
    # knit into GIF
    giffile = os.path.join(data_dir, os.path.splitext(vidfile)[0] + ".gif")
    cmd = "convert -delay 20 -loop 0 $(ls -1 %s/out_frame*png | sort -V) %s" % (data_dir, giffile)
    os.system(cmd)

    # delete intermediate pngs
    if delete_intermediate:
        for tmp_file in os.listdir(data_dir):
            if tmp_file.startswith('out_frame') and tmp_file.endswith('.png'):
                os.remove(os.path.join(data_dir, tmp_file))

elif from_fmt == 'png':
    pngs = [f for f in os.listdir(data_dir) if (f.startswith('global') and
                                                f.endswith('png'))]
    frame_ct = len(pngs)
    for frame_num in range(frame_ct):
        next_png_file = [png for png in pngs if re.search(
                            'ImgDay%i.png' % frame_num, png)]
        assert len(next_png_file) == 1
        next_png_file = next_png_file[0]
        print('\n\tprocessing file number %i: %s...\n\n' % (frame_num, next_png_file))
        RGB_img = imageio.imread(os.path.join(data_dir, next_png_file))
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(16,8)
        ax.imshow(RGB_img)
        plot_calendar_bar(ax, frame_num+1, frame_ct)
        plt.subplots_adjust(left=subplots_adj_left,
                            bottom=subplots_adj_bottom,
                            right=subplots_adj_right,
                            top=subplots_adj_top,
                            wspace=subplots_adj_wspace,
                            hspace=subplots_adj_hspace)
        fig.savefig(os.path.join(data_dir, 'out_frame_%i.png' % frame_num),
                    dpi=700,
                    pad_inches=0,
                    orientation='landscape')
        plt.close('all')

    # turn into movie
    file_basename = re.split('ImgDay\d+\.png', pngs[0])[0]
    avifile = os.path.join(data_dir, file_basename + "_NORMALIZED.avi")
    out_ims = [f for f in os.listdir(data_dir) if
                re.search('out_frame_\d+\.png', f)]
    # sort in day order
    out_ims.sort(key=lambda f: int(re.sub('\D', '', f)))
    height, width, layers = RGB_img.shape
    video = cv2.VideoWriter(avifile, 0, 1, (width,height))
    for im in out_ims:
        video.write(cv2.imread(os.path.join(data_dir, im)))
    cv2.destroyAllWindows()
    video.release()


    # knit into GIF
    giffile = os.path.join(data_dir, file_basename + "_NORMALIZED.gif")
    cmd = "convert -delay 10 -loop 0 $(ls -1 %s/out_frame*png | sort -V) %s" % (data_dir, giffile)
    os.system(cmd)

    # delete intermediate pngs
    if delete_intermediate:
        for tmp_file in os.listdir(data_dir):
            if tmp_file.startswith('out_frame') and tmp_file.endswith('.png'):
                os.remove(os.path.join(data_dir, tmp_file))

print('Yeehaw!')
