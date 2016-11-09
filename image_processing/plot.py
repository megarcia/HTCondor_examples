"""
Python module "plot.py"
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

example script for demonstration on UW HTCondor
"""


import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def plot_grid(img_grid, UTM_zone, UTM_bounds, titlestr, filename,
              cmin=0, cmax=0, zoom=0, xmin=0, xmax=0, ymin=0, ymax=0):
    UTM_W = UTM_bounds[0] / 1000.0
    UTM_N = UTM_bounds[1] / 1000.0
    UTM_E = UTM_bounds[2] / 1000.0
    UTM_S = UTM_bounds[3] / 1000.0
    if zoom:
        img_xmin = xmin
        img_xmax = xmax
        img_ymin = ymin
        img_ymax = ymax
    else:
        img_xmin = UTM_W
        img_xmax = UTM_E
        img_ymin = UTM_S
        img_ymax = UTM_N
    stretch = (img_xmax - img_xmin) / (img_ymax - img_ymin)
    img_mod = np.where(img_grid == -9999.0, np.nan, img_grid)
    plt.imshow(img_mod, extent=(UTM_W, UTM_E, UTM_S, UTM_N), aspect=stretch,
               interpolation='nearest', cmap=plt.get_cmap('jet'))
    if (cmin != 0) or (cmax != 0):
        plt.clim(cmin, cmax)
    plt.colorbar()
    plt.xlim([img_xmin, img_xmax])
    plt.ylim([img_ymin, img_ymax])
    plt.xlabel('UTM %dN easting (km)' % UTM_zone)
    plt.yticks(rotation='vertical')
    plt.ylabel('UTM %dN northing (km)' % UTM_zone)
    plt.title(titlestr, fontsize=12)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.clf()
    return

# end plot.py
