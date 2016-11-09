"""
Python script "chunk_L57_bands.py"
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


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print(char_string)
    sys.stdout.flush()
    return


message(' ')
message('chunk_L57_bands.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 3:
    message('input error: need number of chunks')
    sys.exit(1)
else:
    nchunks = int(sys.argv[2])
#
if len(sys.argv) < 2:
    message('input error: need file name')
    sys.exit(1)
else:
    h5infname = sys.argv[1]
#
datestr = h5infname[:8]
#
with hdf.File(h5infname, 'r') as h5infile:
    projection_tags = np.copy(h5infile['meta/projection_tags'])
    projection = np.copy(h5infile['meta/projection'])
    clip_bounds_tags = np.copy(h5infile['meta/clip_bounds_tags'])
    clip_bounds = np.copy(h5infile['meta/clip_bounds'])
    b1 = np.copy(h5infile['b1_refl_scswmask'])
    b2 = np.copy(h5infile['b2_refl_scswmask'])
    b3 = np.copy(h5infile['b3_refl_scswmask'])
    b4 = np.copy(h5infile['b4_refl_scswmask'])
    b5 = np.copy(h5infile['b5_refl_scswmask'])
    b7 = np.copy(h5infile['b7_refl_scswmask'])
#
nrows, ncols = np.shape(b1)
chunkrows = nrows // nchunks
#
for i in range(nchunks):
    start_row = chunkrows * i
    end_row = chunkrows * (i + 1) - 1
    if i == nchunks - 1:  # last chunk, make sure to catch all that's left
        end_row = nrows - 1
    h5outfname = '%s_L57_rows_%s-%s.h5' % \
        (datestr, str(start_row).zfill(4), str(end_row).zfill(4))
    with hdf.File(h5outfname, 'w') as h5outfile:
        h5outfile.create_dataset('orig_file', data=h5infname)
        h5outfile.create_dataset('start_row', data=start_row)
        h5outfile.create_dataset('end_row', data=end_row)
        h5outfile.create_dataset('meta/projection_tags',
                                 data=projection_tags)
        h5outfile.create_dataset('meta/projection', data=projection)
        h5outfile.create_dataset('meta/clip_bounds_tags',
                                 data=clip_bounds_tags)
        h5outfile.create_dataset('meta/clip_bounds', data=clip_bounds)
        h5outfile.create_dataset('b1', data=b1[start_row:end_row + 1, :],
                                 dtype=np.float32, compression='gzip')
        h5outfile.create_dataset('b2', data=b2[start_row:end_row + 1, :],
                                 dtype=np.float32, compression='gzip')
        h5outfile.create_dataset('b3', data=b3[start_row:end_row + 1, :],
                                 dtype=np.float32, compression='gzip')
        h5outfile.create_dataset('b4', data=b4[start_row:end_row + 1, :],
                                 dtype=np.float32, compression='gzip')
        h5outfile.create_dataset('b5', data=b5[start_row:end_row + 1, :],
                                 dtype=np.float32, compression='gzip')
        h5outfile.create_dataset('b7', data=b7[start_row:end_row + 1, :],
                                 dtype=np.float32, compression='gzip')
    message('wrote %s' % h5outfname)
#
message(' ')
message('chunk_L57_bands.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end chunk_L57_bands.py
