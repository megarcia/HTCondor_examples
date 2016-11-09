"""
Python script "stitch_vi.py"
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

example script for demonstration on UW HTCondor
"""


import sys
import datetime
import glob
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
message('stitch_vi.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 3:
    message('input error: need variable name')
    sys.exit(1)
else:
    varname = sys.argv[2]
#
if len(sys.argv) < 2:
    message('input error: need date')
    sys.exit(1)
else:
    datestr = sys.argv[1]
#
chunkfnames = sorted(glob.glob('%s_%s_*.h5' % (datestr, varname)))
#
with hdf.File(chunkfnames[0], 'r') as h5infile:
    orig_file = np.copy(h5infile['orig_file'])
    projection_tags = np.copy(h5infile['meta/projection_tags'])
    projection = np.copy(h5infile['meta/projection'])
    clip_bounds_tags = np.copy(h5infile['meta/clip_bounds_tags'])
    clip_bounds = np.copy(h5infile['meta/clip_bounds'])
nrows = int(clip_bounds[9])
ncols = int(clip_bounds[8])
complete_check = np.zeros((nrows))
var_grid = np.zeros((nrows, ncols))
#
for h5infname in chunkfnames:
    with hdf.File(h5infname, 'r') as h5infile:
        start_row = np.copy(h5infile['start_row'])
        end_row = np.copy(h5infile['end_row'])
        var_chunk = np.copy(h5infile[varname])
    complete_check[start_row:end_row + 1] = 1
    var_grid[start_row:end_row + 1, :] = var_chunk[:, :]
#
if sum(complete_check) < nrows:
    message('at least one chunk file covering %d rows is missing' %
            (nrows - sum(complete_check)))
    message('- stopping')
    sys.exit(1)
#
h5outfname = '%s_%s.h5' % (datestr, varname)
with hdf.File(h5outfname, 'w') as h5outfile:
    h5outfile.create_dataset('orig_file', data=orig_file)
    h5outfile.create_dataset('meta/projection_tags', data=projection_tags)
    h5outfile.create_dataset('meta/projection', data=projection)
    h5outfile.create_dataset('meta/clip_bounds_tags', data=clip_bounds_tags)
    h5outfile.create_dataset('meta/clip_bounds', data=clip_bounds)
    h5outfile.create_dataset(varname, data=var_grid, dtype=np.float32,
                             compression='gzip')
message('wrote %s' % h5outfname)
message(' ')
#
UTM_zone = int(projection[1].tolist())
UTM_bounds = clip_bounds[0:4]
titlestr = '%s %s' % (datestr, varname.upper())
figfname = '%s_%s.png' % (datestr, varname)
plot_grid(var_grid, UTM_zone, UTM_bounds, titlestr, figfname, cmax=0.6)
message('- saved figure as %s' % figfname)
#
message(' ')
message('stitch_vi.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end stitch_vi.py
