"""
Python script 'clean_GHCND.py'
by Matthew Garcia, PhD student
Dept. of Forest and Wildlife Ecology
University of Wisconsin - Madison
matt.e.garcia@gmail.com

Copyright (C) 2015-2016 by Matthew Garcia
Licensed Gnu GPL v3; see 'LICENSE_GnuGPLv3.txt' for complete terms
Send questions, bug reports, any related requests to matt.e.garcia@gmail.com
See also 'README.md', 'DISCLAIMER.txt', 'CITATION.txt', 'ACKNOWLEDGEMENTS.txt'
Treat others as you would be treated. Pay it forward. Valar dohaeris.

PURPOSE: QA/QC of daily meteorological station data in GHCND datasets

DEPENDENCIES: h5py, numpy, pandas

USAGE: '$ python clean_GHCND.py GHCND_WIS_2015.csv ./.'

INPUT: Station meteorological data from GHCND in '.csv' format (one file)

NOTE: The labels in <metvals> (line 104) are the minimum information that you
      should request from the GHCND/GHCND data server. Check your '.csv' file
      header line to make sure these columns are present. Additional columns
      will be removed (inside this routine, not from your original data csv
      file) to reduce memory footprint.

OUTPUT: One '.csv' file with the 'cleaned' version of the input dataset
        One '.csv' file with an accounting of the 'errors' cleaned,
            by station and variable
        One '.csv' file with station metadata information
"""


import sys
import datetime
import numpy as np
import pandas as pd


def message(char_string):
    """
    prints a string to the terminal and flushes the buffer
    """
    print(char_string)
    sys.stdout.flush()
    return


def find_outliers(dset, stnlist, idxlist, nsigma):
    """
    Simplified (simplistic?) implementation of Chauvenet's criterion for
    outlier detection
    Reference: https://en.wikipedia.org/wiki/Chauvenet%27s_criterion
    Primary assumption: normally distributed data (may not be correct!)
    nsigma >= 1.96 --> p < 0.05  --> 95% chance of data being an outlier
              2.58 --> p < 0.01  --> 99%
              2.81 --> p < 0.005 --> 99.5%
              3.27 --> p < 0.001 --> 99.9%
    """
    outlier_stns = []
    outlier_idxs = []
    nobs = len(dset)
    dset_mean = np.mean(dset)
    dset_std = np.std(dset)
    dset_z = (dset - dset_mean) / dset_std
    for j in range(nobs):
        if abs(dset_z[j]) > nsigma:
            outlier_stns.append(stnlist[j])
            outlier_idxs.append(int(idxlist[j]))
    return outlier_stns, outlier_idxs

outlier_threshold = 3.27
#
message(' ')
message('clean_GHCND.py started at %s' % datetime.datetime.now().isoformat())
message(' ')
#
if len(sys.argv) < 3:
    message('input warning: no data directory path indicated, using ./data')
    path = './data'
else:
    path = sys.argv[2]
#
if len(sys.argv) < 2:
    message('input error: need CSV file containing GHCND weather data')
    sys.exit(1)
else:
    GHCNDfname = '%s/%s' % (path, sys.argv[1])
cleaneddatafile = '%s_cleaned.csv' % GHCNDfname[:-4]
errorsdatafile = '%s_errors.csv' % GHCNDfname[:-4]
stnmetadatafile = '%s_stnmeta.csv' % GHCNDfname[:-4]
#
message('reading input data file %s' % GHCNDfname)
stndata_df = pd.read_csv(GHCNDfname, low_memory=False)
ndatarows, ndatacols = np.shape(stndata_df)
message('- read %d total data rows with %d columns' % (ndatarows, ndatacols))
#
metvals = ['STATION', 'STATION_NAME', 'ELEVATION', 'LATITUDE', 'LONGITUDE',
           'DATE', 'PRCP', 'PRCP_M_FLAG', 'PRCP_Q_FLAG', 'TMAX', 'TMAX_M_FLAG',
           'TMAX_Q_FLAG', 'TMIN', 'TMIN_M_FLAG', 'TMIN_Q_FLAG']
idxs = list(stndata_df.columns.values)
for i, var in enumerate(idxs):
    if var == 'PRCP':
        if 'Measurement Flag' in idxs[i + 1]:
            idxs[i + 1] = 'PRCP_M_FLAG'
        if 'Quality Flag' in idxs[i + 2]:
            idxs[i + 2] = 'PRCP_Q_FLAG'
    if var == 'TMAX':
        if 'Measurement Flag' in idxs[i + 1]:
            idxs[i + 1] = 'TMAX_M_FLAG'
        if 'Quality Flag' in idxs[i + 2]:
            idxs[i + 2] = 'TMAX_Q_FLAG'
    if var == 'TMIN':
        if 'Measurement Flag' in idxs[i + 1]:
            idxs[i + 1] = 'TMIN_M_FLAG'
        if 'Quality Flag' in idxs[i + 2]:
            idxs[i + 2] = 'TMIN_Q_FLAG'
stndata_df.columns = idxs
idxs = list(stndata_df.columns.values)
idxs_extras = list(set(idxs) - set(metvals))
idxs_missing = list(set(metvals) - (set(idxs) - set(idxs_extras)))
if len(idxs_missing) > 0:
    message('NOTE: 1+ necessary data columns is absent from your dataset')
    message('      columns needed: %s' % str(metvals))
    message('      columns missing: %s' % str(idxs_missing))
    sys.exit(1)
#
for idx in idxs:
    if idx not in metvals:
        stndata_df = stndata_df.drop(idx, axis=1)
ndatarows, ndatacols = np.shape(stndata_df)
message('-- reduced to %d columns' % ndatacols)
#
# sort by station and date, assign index values, and add matching IDX column
stndata_df = stndata_df.sort_values(by=['STATION', 'DATE'])
stndata_df.index = pd.Index(np.arange(ndatarows))
stndata_df['IDX'] = np.arange(ndatarows)
message('-- assigned index and added index column')
message(' ')
#
message('original input dataset has:')
ndatarows, ndatacols = np.shape(stndata_df)
message('- %d total data rows' % ndatarows)
#
# dataset summary by stations
stns_all = stndata_df['STATION']
stn_id = list(sorted(stns_all.unique()))
message('- %d unique stations' % len(stn_id))
stn_counts = pd.value_counts(stns_all, sort=True)
min_dates = min(stn_counts)
n_min_dates = sum(stn_counts == min_dates)
message('-- minimum %d date entries for %d stations' %
        (min_dates, n_min_dates))
max_dates = max(stn_counts)
n_max_dates = sum(stn_counts == max_dates)
message('-- maximum %d date entries for %d stations' %
        (max_dates, n_max_dates))
avg_dates = np.mean(stn_counts)
message('-- average %.1f date entries per station' % avg_dates)
med_dates = np.median(stn_counts)
message('-- median %d date entries per station' % med_dates)
#
# dataset summary by dates
dates_all = stndata_df['DATE']
dates = list(sorted(dates_all.unique()))
message('- %d unique dates from %s to %s' %
        (len(dates), str(dates[0]), str(dates[-1])))
date_counts = pd.value_counts(dates_all, sort=True)
min_stns = min(date_counts)
n_min_stns = sum(date_counts == min_stns)
message('-- minimum %d station entries for %d dates' % (min_stns, n_min_stns))
max_stns = max(date_counts)
n_max_stns = sum(date_counts == max_stns)
message('-- maximum %d station entries for %d dates' % (max_stns, n_max_stns))
avg_stns = np.mean(date_counts)
message('-- average %.1f station entries per date' % avg_stns)
med_stns = np.median(date_counts)
message('-- median %d station entries per date' % med_stns)
message(' ')
#
message('cleaning input dataset')
message(' ')
noperations = 0
#
# drop entries with no location information
unk_lats_idxs = list(stndata_df.ix[stndata_df['LATITUDE'] == 'unknown', 'IDX'])
unk_lats_stns = list(stndata_df.ix[stndata_df['LATITUDE'] ==
                                   'unknown', 'STATION'].unique())
message('found %d records at %d stations with unknown latitude' %
        (len(unk_lats_idxs), len(unk_lats_stns)))
unk_lons_idxs = list(stndata_df.ix[stndata_df['LONGITUDE'] ==
                                   'unknown', 'IDX'])
unk_lons_stns = list(stndata_df.ix[stndata_df['LONGITUDE'] ==
                                   'unknown', 'STATION'].unique())
message('found %d records at %d stations with unknown longitude' %
        (len(unk_lons_idxs), len(unk_lons_stns)))
# set union, either lat or lon is missing
unk_locs_idxs = list(set(unk_lats_idxs) | set(unk_lons_idxs))
if len(unk_locs_idxs) > 0:
    stndata_df = stndata_df.drop(unk_locs_idxs)
message('dropped %d total records for lack of location information' %
        len(unk_locs_idxs))
noperations += len(unk_locs_idxs)
#
# re-index
stndata_df = stndata_df.drop('IDX', axis=1)
stndata_df = stndata_df.sort_values(by=['STATION', 'DATE'])
ndatarows, ndatacols = np.shape(stndata_df)
stns_all = stndata_df['STATION']
stn_id = list(sorted(stns_all.unique()))
dates_all = stndata_df['DATE']
dates = list(sorted(dates_all.unique()))
stndata_df.index = pd.Index(np.arange(ndatarows))
stndata_df['IDX'] = np.arange(ndatarows)
message('re-indexed remaining data records')
message(' ')
#
# *actual* null value is -9999, not 9999 as listed in the GHCND documentation
# find and adjust any 9999 observation values to -9999
prcp_null_idx = list(stndata_df.ix[stndata_df['PRCP'] == 9999, 'IDX'])
if len(prcp_null_idx) > 0:
    stndata_df.set_value((prcp_null_idx), 'PRCP', -9999)
    message('PRCP = 9999 set to PRCP = -9999 for %d observations' %
            len(prcp_null_idx))
    noperations += len(prcp_null_idx)
tmax_null_idx = list(stndata_df.ix[stndata_df['TMAX'] == 9999, 'IDX'])
if len(tmax_null_idx) > 0:
    stndata_df.set_value((tmax_null_idx), 'TMAX', -9999)
    message('TMAX = 9999 set to TMAX = -9999 for %d observations' %
            len(tmax_null_idx))
    noperations += len(tmax_null_idx)
tmin_null_idx = list(stndata_df.ix[stndata_df['TMIN'] == 9999, 'IDX'])
if len(tmin_null_idx) > 0:
    stndata_df.set_value((tmin_null_idx), 'TMIN', -9999)
    message('TMIN = 9999 set to TMIN = -9999 for %d observations' %
            len(tmin_null_idx))
    noperations += len(tmin_null_idx)
#
# set up DataFrame for error accounting and reporting
stnerr_df = pd.DataFrame(np.zeros((len(stn_id), 13)).astype(int),
                         index=np.arange(len(stn_id)),
                         columns=['STATION', 'IDX', 'PRCP_T_ADJ', 'PRCP_M_ERR',
                                  'PRCP_Q_ERR', 'PRCP_ZERO_ERR', 'TMAX_M_ERR',
                                  'TMAX_Q_ERR', 'TMIN_M_ERR', 'TMIN_Q_ERR',
                                  'T_REV_ERR', 'TMAX_OUTLIER', 'TMIN_OUTLIER'])
stnerr_df['STATION'] = stn_id
stnerr_df['IDX'] = np.arange(len(stn_id))
message('established error accounting table')
message(' ')
#
# process PRCP measurement flags
# - flags 'B' and 'D' are ok
# - flag 'T' means that a value of 0 *should* be 1 (= 0.1 mm ==> 0.01 cm)
# - flag 'P' means that value *should* be missing (-9999 instead of 0)
prcp_flags_all = stndata_df['PRCP_M_FLAG']
prcp_flags = list(sorted(prcp_flags_all.unique()))
message('found %d unique PRCP_M_FLAGS: %s' %
        (len(prcp_flags), str(prcp_flags)))
if len(prcp_flags) > 1:
    if 'B' in prcp_flags:
        message('- NOTE: PRCP flag B found, no data adjustments were made')
    if 'D' in prcp_flags:
        message('- NOTE: PRCP flag D found, no data adjustments were made')
    if 'T' in prcp_flags:
        prcp_t_idx = list(stndata_df.ix[stndata_df['PRCP_M_FLAG'] ==
                                        'T', 'IDX'])
        stndata_df.set_value((prcp_t_idx), 'PRCP', 1)
        message('- trace PRCP --> PRCP = 1 (= 0.01 cm) for %d observations' %
                len(prcp_t_idx))
        noperations += len(prcp_t_idx)
        prcp_t_stns_all = stndata_df.ix[stndata_df['PRCP_M_FLAG'] ==
                                        'T', 'STATION']
        prcp_t_stns = list(prcp_t_stns_all.unique())
        message('-- accounting PRCP adjustments for %d stations' %
                len(prcp_t_stns))
        for stn in prcp_t_stns:
            stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
            stnerr_count = sum(prcp_t_stns_all == stn)
            stnerr_df.set_value((stnerr_idx), 'PRCP_T_ADJ', stnerr_count)
    if 'P' in prcp_flags:
        prcp_m_idx = list(stndata_df.ix[stndata_df['PRCP_M_FLAG'] ==
                                        'P', 'IDX'])
        stndata_df.set_value((prcp_m_idx), 'PRCP', -9999)
        message('- presumed PRCP = 0 --> PRCP = -9999 for %d observations' %
                len(prcp_m_idx))
        noperations += len(prcp_m_idx)
        prcp_m_stns_all = stndata_df.ix[stndata_df['PRCP_M_FLAG'] ==
                                        'P', 'STATION']
        prcp_m_stns = list(prcp_m_stns_all.unique())
        message('-- accounting PRCP adjustments for %d stations' %
                len(prcp_m_stns))
        for stn in prcp_m_stns:
            stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
            stnerr_count = sum(prcp_m_stns_all == stn)
            stnerr_df.set_value((stnerr_idx), 'PRCP_M_ERR', stnerr_count)
    if len(prcp_flags) > 5:
        message('- NOTE: PRCP flag other than B/D/T/P was found, \
                but no data adjustments were made')
message(' ')
#
# process PRCP quality flags
# - anything but a blank here means that the observation failed one of
#   NOAA/GHCND's own QA checks
prcp_flags_all = stndata_df['PRCP_Q_FLAG']
prcp_flags = list(sorted(prcp_flags_all.unique()))
message('found %d unique PRCP_Q_FLAGS: %s' %
        (len(prcp_flags), str(prcp_flags)))
if len(prcp_flags) > 1:
    prcp_q_idx = list(stndata_df.ix[stndata_df['PRCP_Q_FLAG'] != ' ', 'IDX'])
    stndata_df.set_value((prcp_q_idx), 'PRCP', -9999)
    message('- QA-failed PRCP set to PRCP = -9999 for %d observations' %
            len(prcp_q_idx))
    noperations += len(prcp_q_idx)
    prcp_q_stns_all = stndata_df.ix[stndata_df['PRCP_Q_FLAG'] !=
                                    ' ', 'STATION']
    prcp_q_stns = list(prcp_q_stns_all.unique())
    message('-- accounting QA-based PRCP adjustments for %d stations' %
            len(prcp_q_stns))
    for stn in prcp_q_stns:
        stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
        stnerr_count = sum(prcp_q_stns_all == stn)
        stnerr_df.set_value((stnerr_idx), 'PRCP_Q_ERR', stnerr_count)
message(' ')
#
# process TMAX measurement flags
# - flag 'L' means that the observation may have been recorded some time later
#   than the actual TMAX occurrence
tmax_flags_all = stndata_df['TMAX_M_FLAG']
tmax_flags = list(sorted(tmax_flags_all.unique()))
message('found %d unique TMAX_M_FLAGS: %s' %
        (len(tmax_flags), str(tmax_flags)))
if len(tmax_flags) > 1:
    if 'L' in tmax_flags:
        tmax_m_idx = list(stndata_df.ix[stndata_df['TMAX_M_FLAG'] ==
                                        'L', 'IDX'])
        stndata_df.set_value((tmax_m_idx), 'TMAX', -9999)
        message('- lagged TMAX set to TMAX = -9999 for %d observations' %
                len(tmax_m_idx))
        noperations += len(tmax_m_idx)
        tmax_m_stns_all = stndata_df.ix[stndata_df['TMAX_M_FLAG'] ==
                                        'L', 'STATION']
        tmax_m_stns = list(tmax_m_stns_all.unique())
        message('-- accounting TMAX adjustments for %d stations' %
                len(tmax_m_stns))
        for stn in tmax_m_stns:
            stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
            stnerr_count = sum(tmax_m_stns_all == stn)
            stnerr_df.set_value((stnerr_idx), 'TMAX_M_ERR', stnerr_count)
    if (len(tmax_flags) > 2) or ('L' not in tmax_flags):
        message('- NOTE: a TMAX flag other than L was found, but no data \
                adjustments were made')
message(' ')
#
# process TMAX quality flags
# - anything but a blank here means that the observation failed one of
#   GHCND's own QA checks
tmax_flags_all = stndata_df['TMAX_Q_FLAG']
tmax_flags = list(sorted(tmax_flags_all.unique()))
n_tmax_flags = len(tmax_flags)
message('found %d unique TMAX_Q_FLAGS: %s' % (n_tmax_flags, str(tmax_flags)))
if n_tmax_flags > 1:
    tmax_q_idx = list(stndata_df.ix[stndata_df['TMAX_Q_FLAG'] != ' ', 'IDX'])
    stndata_df.set_value((tmax_q_idx), 'TMAX', -9999)
    message('- QA-failed TMAX set to TMAX = -9999 for %d observations' %
            len(tmax_q_idx))
    noperations += len(tmax_q_idx)
    tmax_q_stns_all = stndata_df.ix[stndata_df['TMAX_Q_FLAG'] !=
                                    ' ', 'STATION']
    tmax_q_stns = list(tmax_q_stns_all.unique())
    message('-- accounting QA-based TMAX adjustments for %d stations' %
            len(tmax_q_stns))
    for stn in tmax_q_stns:
        stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
        stnerr_count = sum(tmax_q_stns_all == stn)
        stnerr_df.set_value((stnerr_idx), 'TMAX_Q_ERR', stnerr_count)
message(' ')
#
# process TMIN measurement flags
# - flag 'L' means that the observation may have been recorded some time later
#   than the actual TMIN occurrence
tmin_flags_all = stndata_df['TMIN_M_FLAG']
tmin_flags = list(sorted(tmin_flags_all.unique()))
n_tmin_flags = len(tmin_flags)
message('found %d unique TMIN_M_FLAGS: %s' % (n_tmin_flags, str(tmin_flags)))
if n_tmin_flags > 1:
    if 'L' in tmin_flags:
        tmin_m_idx = list(stndata_df.ix[stndata_df['TMIN_M_FLAG'] ==
                                        'L', 'IDX'])
        stndata_df.set_value((tmin_m_idx), 'TMIN', -9999)
        message('- lagged TMIN set to TMIN = -9999 for %d observations' %
                len(tmin_m_idx))
        noperations += len(tmin_m_idx)
        tmin_m_stns_all = stndata_df.ix[stndata_df['TMIN_M_FLAG'] ==
                                        'L', 'STATION']
        tmin_m_stns = list(tmin_m_stns_all.unique())
        message('-- accounting TMIN adjustments for %d stations' %
                len(tmin_m_stns))
        for stn in tmin_m_stns:
            stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
            stnerr_count = sum(tmin_m_stns_all == stn)
            stnerr_df.set_value((stnerr_idx), 'TMIN_M_ERR', stnerr_count)
    if (n_tmin_flags > 2) or ('L' not in tmin_flags):
        message('- NOTE: a TMIN flag other than L was found, but no data \
                adjustments were made')
message(' ')
#
# process TMIN quality flags
# - anything but a blank here means that the observation failed one of
#   GHCND's own QA checks
tmin_flags_all = stndata_df['TMIN_Q_FLAG']
tmin_flags = list(sorted(tmin_flags_all.unique()))
n_tmin_flags = len(tmin_flags)
message('found %d unique TMIN_Q_FLAGS: %s' % (n_tmin_flags, str(tmin_flags)))
if n_tmin_flags > 1:
    tmin_q_idx = stndata_df.ix[stndata_df['TMIN_Q_FLAG'] != ' ', 'IDX']
    stndata_df.set_value((tmin_q_idx), 'TMIN', -9999)
    message('- QA-failed TMIN set to TMIN = -9999 for %d observations' %
            len(tmin_q_idx))
    noperations += len(tmin_q_idx)
    tmin_q_stns_all = stndata_df.ix[stndata_df['TMIN_Q_FLAG'] !=
                                    ' ', 'STATION']
    tmin_q_stns = list(tmin_q_stns_all.unique())
    message('-- accounting QA-based TMIN adjustments for %d stations' %
            len(tmin_q_stns))
    for stn in tmin_q_stns:
        stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
        stnerr_count = sum(tmin_q_stns_all == stn)
        stnerr_df.set_value((stnerr_idx), 'TMIN_Q_ERR', stnerr_count)
message(' ')
#
# find stations with potentially erroneous P values
# (evaluation criterion: all of that station's P values = 0)
# (outlier detection of extreme values not yet implemented)
message('checking PRCP values for %d individual stations' % len(stn_id))
for stn in stn_id:
    stn_df = stndata_df[stndata_df['STATION'] == stn]
    stn_idxs = list(stn_df['IDX'])
    stn_df = stn_df[stn_df['PRCP'] != -9999]
    stn_prcp = np.array(stn_df['PRCP'])
    if sum(stn_prcp) <= 0:
        message('- station %s may have erroneous PRCP data, setting all %d of \
                its PRCP values to -9999' % (stn, len(stn_idxs)))
        stndata_df.set_value((stn_idxs), 'PRCP', -9999)
        noperations += len(stn_idxs)
        stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
        stnerr_df.set_value((stnerr_idx), 'PRCP_ZERO_ERR', len(stn_idxs))
    else:
        message('- station %s PRCP data seem OK (no check for extreme values \
                yet)' % stn)
message(' ')
#
# find entries with potentially erroneous locations/values
# (evaluation criteria: reversed Tmax and Tmin values)
# (                     outlier detection of extreme T values)
message('checking PRCP/TMAX/TMIN values for %d individual dates' % len(dates))
for date in dates:
    date_df = stndata_df[stndata_df['DATE'] == date]
    stn_vals = list(date_df['STATION'])
    lats_vals = list(date_df['LATITUDE'])
    lons_vals = list(date_df['LONGITUDE'])
    prcp_vals = np.array(date_df['PRCP'])
    tmax_vals = np.array(date_df['TMAX'])
    tmin_vals = np.array(date_df['TMIN'])
    idx_vals = list(date_df['IDX'])
    #
    # find transposed Tmax/Tmin values
    # (only works if valid values are present for both)
    nvals = len(idx_vals)
    transposed = 0
    for i in range(nvals):
        if (tmax_vals[i] != -9999) and (tmin_vals[i] != -9999):
            if tmin_vals[i] > tmax_vals[i]:
                stndata_df.set_value((int(idx_vals[i])), 'TMAX', tmin_vals[i])
                stndata_df.set_value((int(idx_vals[i])), 'TMIN', tmax_vals[i])
                transposed += 1
                stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] ==
                                          stn_vals[i], 'IDX']
                stnerr_val = stnerr_df.ix[stnerr_df['STATION'] ==
                                          stn_vals[i], 'T_REV_ERR']
                stnerr_df.set_value((stnerr_idx), 'T_REV_ERR',
                                    (stnerr_val + 1))
    if transposed > 0:
        message('- date %d has %d entries where Tmax < Tmin, now reversed' %
                (date, transposed))
        noperations += transposed
    #
    # find outlier Tmax values (not spatial, purely arithmetic)
    valid_tmax_vals = date_df[date_df['TMAX'] != -9999]
    tmax_vals = np.array(valid_tmax_vals['TMAX'])
    stn_vals = list(valid_tmax_vals['STATION'])
    idx_vals = list(valid_tmax_vals['IDX'])
    tmax_outlier_stns, tmax_outlier_idxs = find_outliers(tmax_vals, stn_vals,
                                                         idx_vals,
                                                         outlier_threshold)
    if len(tmax_outlier_stns) > 0:
        stndata_df.set_value((tmax_outlier_idxs), 'TMAX', -9999)
        for stn in tmax_outlier_stns:
            stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
            stnerr_val = stnerr_df.ix[stnerr_df['STATION'] ==
                                      stn, 'TMAX_OUTLIER']
            stnerr_df.set_value((stnerr_idx), 'TMAX_OUTLIER', (stnerr_val + 1))
        message('- date %d has %d outlier Tmax values (p < 0.001), now set to \
                -9999' % (date, len(tmax_outlier_stns)))
        noperations += len(tmax_outlier_stns)
    #
    # find outlier Tmin values (not spatial, purely arithmetic)
    valid_tmin_vals = date_df[date_df['TMIN'] != -9999]
    tmin_vals = np.array(valid_tmin_vals['TMIN'])
    stn_vals = list(valid_tmin_vals['STATION'])
    idx_vals = list(valid_tmin_vals['IDX'])
    tmin_outlier_stns, tmin_outlier_idxs = find_outliers(tmin_vals, stn_vals,
                                                         idx_vals,
                                                         outlier_threshold)
    if len(tmin_outlier_stns) > 0:
        stndata_df.set_value((tmin_outlier_idxs), 'TMAX', -9999)
        for stn in tmin_outlier_stns:
            stnerr_idx = stnerr_df.ix[stnerr_df['STATION'] == stn, 'IDX']
            stnerr_val = stnerr_df.ix[stnerr_df['STATION'] ==
                                      stn, 'TMIN_OUTLIER']
            stnerr_df.set_value((stnerr_idx), 'TMIN_OUTLIER', (stnerr_val + 1))
        message('- date %d has %d outlier Tmin values (p < 0.001), now set to \
                -9999' % (date, len(tmin_outlier_stns)))
        noperations += len(tmin_outlier_stns)
    #
    # if no apparent errors/outliers found, report that
    if (transposed + len(tmax_outlier_stns) + len(tmin_outlier_stns)) == 0:
        message('- date %d data seem OK' % date)
message(' ')
#
# drop entries with no useful meteorological data (just for convenience)
no_prcp_data_idxs = list(stndata_df.ix[stndata_df['PRCP'] == -9999, 'IDX'])
no_tmax_data_idxs = list(stndata_df.ix[stndata_df['TMAX'] == -9999, 'IDX'])
no_tmin_data_idxs = list(stndata_df.ix[stndata_df['TMIN'] == -9999, 'IDX'])
# set intersection, see if all met data are missing
no_data_idxs = list(set(no_prcp_data_idxs) &
                    set(no_tmax_data_idxs) &
                    set(no_tmin_data_idxs))
if len(no_data_idxs) > 0:
    stndata_df = stndata_df.drop(no_data_idxs)
    message('dropped %d records for lack of useful meteorological data' %
            len(no_data_idxs))
    noperations += len(no_data_idxs)
#
# re-index
stndata_df = stndata_df.drop('IDX', axis=1)
stndata_df = stndata_df.sort_values(by=['STATION', 'DATE'])
ndatarows, ndatacols = np.shape(stndata_df)
stndata_df.index = pd.Index(np.arange(ndatarows))
stndata_df['IDX'] = np.arange(ndatarows)
message('re-indexed remaining data records')
message(' ')
#
message('a total of %d data value operations were performed' % noperations)
message(' ')
#
ndatarows, ndatacols = np.shape(stndata_df)
message('cleaned dataset has:')
message('- %d total data rows' % ndatarows)
#
# dataset summary by stations
stns_all = stndata_df['STATION']
stn_id = sorted(stns_all.unique())
message('- %d unique stations' % len(stn_id))
stn_counts = pd.value_counts(stns_all, sort=True)
min_dates = min(stn_counts)
n_min_dates = sum(stn_counts == min_dates)
message('-- minimum %d date entries for %d stations' %
        (min_dates, n_min_dates))
max_dates = max(stn_counts)
n_max_dates = sum(stn_counts == max_dates)
message('-- maximum %d date entries for %d stations' %
        (max_dates, n_max_dates))
avg_dates = np.mean(stn_counts)
message('-- average %.1f date entries per station' % avg_dates)
med_dates = np.median(stn_counts)
message('-- median %d date entries per station' % med_dates)
#
# dataset summary by dates
dates_all = stndata_df['DATE']
dates = sorted(dates_all.unique())
message('- %d unique dates from %s to %s' %
        (len(dates), str(dates[0]), str(dates[-1])))
date_counts = pd.value_counts(dates_all, sort=True)
min_stns = min(date_counts)
n_min_stns = sum(date_counts == min_stns)
message('-- minimum %d station entries for %d dates' %
        (min_stns, n_min_stns))
max_stns = max(date_counts)
n_max_stns = sum(date_counts == max_stns)
message('-- maximum %d station entries for %d dates' %
        (max_stns, n_max_stns))
avg_stns = np.mean(date_counts)
message('-- average %.1f station entries per date' % avg_stns)
med_stns = np.median(date_counts)
message('-- median %d station entries per date' % med_stns)
message(' ')
#
# dataset summary by variable
message('- %d total stations' % len(stn_id))
message('- %d unique dates from %s to %s' %
        (len(dates), str(dates[0]), str(dates[-1])))
lats = np.zeros((len(stn_id)))
lons = np.zeros((len(stn_id)))
ndates = np.zeros((len(stn_id))).astype(int)
nprcp = np.zeros((len(stn_id))).astype(int)
ntmax = np.zeros((len(stn_id))).astype(int)
ntmin = np.zeros((len(stn_id))).astype(int)
for i in range(len(stn_id)):
    message('-- processing station %s' % stn_id[i])
    stndata_subset = stndata_df[stndata_df['STATION'] == stn_id[i]]
    lats[i] = list(stndata_subset['LATITUDE'])[-1]
    lons[i] = list(stndata_subset['LONGITUDE'])[-1]
    ndates[i] = len(list(stndata_subset['DATE']))
    nprcp[i] = sum(list(stndata_subset['PRCP'] != -9999))
    ntmax[i] = sum(list(stndata_subset['TMAX'] != -9999))
    ntmin[i] = sum(list(stndata_subset['TMIN'] != -9999))
stnmeta_df = pd.DataFrame({'STATION': stn_id, 'LATITUDE': lats,
                           'LONGITUDE': lons, 'NDATES': ndates,
                           'NPRCP': nprcp, 'NTMAX': ntmax, 'NTMIN': ntmin})
message(' ')
#
# save cleaned csv dataset
stndata_df.to_csv(cleaneddatafile)
message('saved cleaned GHCND dataset to %s' % cleaneddatafile)
stnerr_df.to_csv(errorsdatafile)
message('saved station-by-station error counts to %s' % errorsdatafile)
stnmeta_df.to_csv(stnmetadatafile)
message('saved station metadataset to %s' % stnmetadatafile)
message(' ')
#
message('clean_GHCND.py completed at %s' % datetime.datetime.now().isoformat())
message(' ')
sys.exit(0)

# end clean_GHCND.py
