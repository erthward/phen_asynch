import asynch_fns as af
DESIGN_MAT = af.make_design_matrix()
DATA_DIR = "../../GEE_output"
mix = af.read_mixer_file(DATA_DIR)
(DIMS, CRS, XMIN, YMIN, XRES, YRES, patches_per_row,
                                    tot_patches) = af.get_mixer_info(mix)
INDICES = af.make_indices_array(DIMS)
PATT_B4_FILENUM = "SIF-"
FILES_DICT = af.get_row_col_patch_ns_allfiles(DATA_DIR, PATT_B4_FILENUM)
INBANDS = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']
OUTBANDS = ['asynch', 'asynch_R2', 'asynch_euc', 'asynch_euc_R2', 'asynch_n']
ip, op = af.get_inpatches_outpatches([*FILES_DICT][74], INBANDS, DIMS)
fi = [*FILES_DICT.items()][74]
row_is = fi[1]['row_is']
col_js = fi[1]['col_js']
patch_ns = fi[1]["patch_ns"]
ops = af.calc_asynch(ip, op, row_is, col_js, patch_ns, DIMS, XMIN, YMIN, XRES, YRES, DESIGN_MAT, INDICES, 30)
