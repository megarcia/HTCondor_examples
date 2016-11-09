"""
Python script "calc_L57_chunk_vi.py"
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


def calc_ndxi(b4, bx):
    b4_mask = np.where(b4 > 0.0, 1, 0)
    bx_mask = np.where(bx > 0.0, 1, 0)
    combined_mask = b4_mask * bx_mask
    ndxi_num = b4 - bx
    ndxi_den = b4 + bx
    ndxi = np.where(combined_mask != 0, ndxi_num / ndxi_den, -9999.0)
    return ndxi


message(' ')
message('calc_L57_chunk_vi.py started at %s' %
        datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 2:
    message('input error: need chunk file name')
    sys.exit(1)
else:
    h5infname = sys.argv[1]
#
datestr = h5infname[:8]
#
with hdf.File(h5infname, 'r') as h5infile:
    orig_file = np.copy(h5infile['orig_file'])
    start_row = np.copy(h5infile['start_row'])
    end_row = np.copy(h5infile['end_row'])
    projection_tags = np.copy(h5infile['meta/projection_tags'])
    projection = np.copy(h5infile['meta/projection'])
    clip_bounds_tags = np.copy(h5infile['meta/clip_bounds_tags'])
    clip_bounds = np.copy(h5infile['meta/clip_bounds'])
    band_4 = np.copy(h5infile['b4'])
    band_5 = np.copy(h5infile['b5'])
    band_7 = np.copy(h5infile['b7'])
#
ndii = calc_ndxi(band_4, band_5)
message('calcuated normalized difference infrared index (NDII)')
#
h5outfname = '%s_ndii_rows_%s-%s.h5' % \
    (datestr, str(start_row).zfill(4), str(end_row).zfill(4))
with hdf.File(h5outfname, 'w') as h5outfile:
    h5outfile.create_dataset('orig_file', data=orig_file)
    h5outfile.create_dataset('start_row', data=start_row)
    h5outfile.create_dataset('end_row', data=end_row)
    h5outfile.create_dataset('meta/projection_tags', data=projection_tags)
    h5outfile.create_dataset('meta/projection', data=projection)
    h5outfile.create_dataset('meta/clip_bounds_tags', data=clip_bounds_tags)
    h5outfile.create_dataset('meta/clip_bounds', data=clip_bounds)
    h5outfile.create_dataset('ndii', data=ndii, dtype=np.float32,
                             compression='gzip')
message('wrote %s' % h5outfname)
message(' ')
#
nbr = calc_ndxi(band_4, band_7)
message('calculated normalized burn ratio (NBR)')
#
h5outfname = '%s_nbr_rows_%s-%s.h5' % \
    (datestr, str(start_row).zfill(4), str(end_row).zfill(4))
with hdf.File(h5outfname, 'w') as h5outfile:
    h5outfile.create_dataset('orig_file', data=h5infname)
    h5outfile.create_dataset('start_row', data=start_row)
    h5outfile.create_dataset('end_row', data=end_row)
    h5outfile.create_dataset('meta/projection_tags', data=projection_tags)
    h5outfile.create_dataset('meta/projection', data=projection)
    h5outfile.create_dataset('meta/clip_bounds_tags', data=clip_bounds_tags)
    h5outfile.create_dataset('meta/clip_bounds', data=clip_bounds)
    h5outfile.create_dataset('nbr', data=nbr, dtype=np.float32,
                             compression='gzip')
message('wrote %s' % h5outfname)
message(' ')
#
message(' ')
message('calc_L57_chunk_vi.py completed at %s' %
        datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end calc_L57_chunk_vi.py
