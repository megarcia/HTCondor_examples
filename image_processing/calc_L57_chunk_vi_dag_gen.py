"""
Python script "calc_L57_chunk_vi_dag_gen.py"
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

example script for demonstration on UW HTCondor
"""


import sys
import glob


datestr = sys.argv[1]
#
chunklist = sorted(glob.glob('%s_L57_rows_*.h5' % datestr))
#
dag_fname = 'calc_L57_chunk_vi_%s_dag.sub' % datestr
with open(dag_fname, 'w') as dag_f:
    for i, fname in enumerate(chunklist):
        jobname = 'A_%s_%s' % (datestr, str(i).zfill(2))
        dag_f.write('JOB %s calc_L57_chunk_vi.sub \n' % jobname)
        dag_f.write('VARS %s date="%s" h5name="%s" chunk="%d" \n' %
                    (jobname, datestr, fname, i))
        dag_f.write('\n')
#
sys.exit(0)

# end calc_L57_chunk_vi_dag_gen.py
