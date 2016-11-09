"""
Python script "diff_vi.py"
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

example script for demonstration on UW HTCondor
"""


import sys
import datetime
import numpy as np
import h5py as hdf
from plot import plot_grid


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print(char_string)
    sys.stdout.flush()
    return


message(' ')
message('diff_vi.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 4:
    message('input error: need variable name')
    sys.exit(1)
else:
    varname = sys.argv[3]
#
if len(sys.argv) < 3:
    message('input error: need second date in pair')
    sys.exit(1)
else:
    date2str = sys.argv[2]
#
if len(sys.argv) < 2:
    message('input error: need first date in pair')
    sys.exit(1)
else:
    date1str = sys.argv[1]
#
date1fname = '%s_%s.h5' % (date1str, varname)
with hdf.File(date1fname, 'r') as h5infile:
    orig_file_1 = np.copy(h5infile['orig_file'])
    projection_tags = np.copy(h5infile['meta/projection_tags'])
    projection_1 = np.copy(h5infile['meta/projection'])
    clip_bounds_tags = np.copy(h5infile['meta/clip_bounds_tags'])
    clip_bounds_1 = np.copy(h5infile['meta/clip_bounds'])
    var_grid_1 = np.copy(h5infile[varname])
grid_1_mask = np.where(var_grid_1 == -9999.0, 0, 1)
#
date2fname = '%s_%s.h5' % (date2str, varname)
with hdf.File(date2fname, 'r') as h5infile:
    orig_file_2 = np.copy(h5infile['orig_file'])
    projection_2 = np.copy(h5infile['meta/projection'])
    clip_bounds_2 = np.copy(h5infile['meta/clip_bounds'])
    var_grid_2 = np.copy(h5infile[varname])
grid_2_mask = np.where(var_grid_2 == -9999.0, 0, 1)
#
grid_combined_mask = grid_1_mask * grid_2_mask
#
if varname == 'nbr':
    diff_grid = np.where(grid_combined_mask == 1,
                         (var_grid_1 - var_grid_2), -9999.0)
    denom_grid = np.where(grid_combined_mask == 1,
                          np.sqrt(abs(var_grid_1 / 1000.0)), 0.0)
    rel_diff_grid = np.where(grid_combined_mask == 1,
                             (diff_grid / denom_grid), -9999.0)
else:
    diff_grid = np.where(grid_combined_mask == 1,
                         (var_grid_2 - var_grid_1), -9999.0)
#
h5outfname = '%s_%s_d%s.h5' % (date1str, date2str, varname)
with hdf.File(h5outfname, 'w') as h5outfile:
    h5outfile.create_dataset('orig_file_1', data=orig_file_1)
    h5outfile.create_dataset('vi_file_1', data=date1fname)
    h5outfile.create_dataset('orig_file_2', data=orig_file_2)
    h5outfile.create_dataset('vi_file_2', data=date2fname)
    h5outfile.create_dataset('meta/projection_tags', data=projection_tags)
    h5outfile.create_dataset('meta/projection', data=projection_1)
    h5outfile.create_dataset('meta/clip_bounds_tags', data=clip_bounds_tags)
    h5outfile.create_dataset('meta/clip_bounds', data=clip_bounds_1)
    dvarname = 'd%s' % varname
    h5outfile.create_dataset(dvarname, data=diff_grid, dtype=np.float32,
                             compression='gzip')
    if varname == 'nbr':
        rdvarname = 'rd%s' % varname
        h5outfile.create_dataset(rdvarname, data=rel_diff_grid,
                                 dtype=np.float32, compression='gzip')
message('wrote %s' % h5outfname)
message(' ')
#
UTM_zone = int(projection_1[1].tolist())
UTM_bounds = clip_bounds_1[0:4]
titlestr = '%s to %s d%s' % (date1str, date2str, varname.upper())
figfname = '%s_%s_d%s.png' % (date1str, date2str, varname)
plot_grid(diff_grid, UTM_zone, UTM_bounds, titlestr, figfname,
          cmin=-0.6, cmax=0.6)
message('- saved figure as %s' % figfname)
#
if varname == 'nbr':
    titlestr = '%s to %s Rd%s' % (date1str, date2str, varname.upper())
    figfname = '%s_%s_rd%s.png' % (date1str, date2str, varname)
    plot_grid(rel_diff_grid, UTM_zone, UTM_bounds, titlestr, figfname,
              cmin=0.0, cmax=40.0)
    message('- saved figure as %s' % figfname)
#
message(' ')
message('diff_vi.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end diff_vi.py
