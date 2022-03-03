import cv2
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
import os

# NOTE: code adapted from:
    # https://stackoverflow.com/questions/32111705/overlay-a-graph-over-a-video

# function to make timebar for a given day
def plot_calendar_bar(ax, frame_num, frame_ct):

    # plot the bar
    p = Polygon(np.array([[100,380],[100,350], [700,350], [700,380], [100,380]]))
    pc = PatchCollection([p], alpha=1, edgecolor='gray', facecolors='#e5e3e6')
    ax.add_collection(pc)

    # plot the months
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    locs = np.linspace(105, 685, len(months))
    for month, loc in zip(months, locs):
        ax.plot([loc+5]*2, [350, 365], '-k', linewidth=0.5)
        ax.text(loc-3, 377, month, color='black', fontsize=12,)
                #fontdict={'fontweight': 'bold'})

    # plot the current timepoint
    timepoint_frac = frame_num/frame_ct
    loc = 100 + ((700-100)*timepoint_frac)
    ax.plot([loc]*2, [376, 354], '-', color='#6d00a3', linewidth=9, alpha=0.5)

subplots_adj_left=0
subplots_adj_bottom=0
subplots_adj_right=1
subplots_adj_top=1
subplots_adj_wspace=0
subplots_adj_hspace=0

data_dir = '/home/deth/Desktop/CAL/research/projects/seasonality/results/videos'

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
    for tmp_file in os.listdir(data_dir):
        if tmp_file.startswith('out_frame') and tmp_file.endswith('.png'):
            os.remove(os.path.join(data_dir, tmp_file))

print('Yeehaw!')
