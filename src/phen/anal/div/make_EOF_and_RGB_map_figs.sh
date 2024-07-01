# NOTE: save the scaled, folded, reprojected EOF raster on the first call only
python ./plot_EOF_and_RGB_results.py main_rgb_map save_rast
python ./plot_EOF_and_RGB_results.py reg_figs
python ./plot_EOF_and_RGB_results.py eof_summ_fig
python ./plot_EOF_and_RGB_results.py raw_rgb_maps

