"""Utilities for parsing QUIP timing data (WIP)"""


import numpy as np


def store_timing_data(job, data, timing_key):
    """Store QUIP timing data and summary statistics in job"""
    for key, timer in data.items():
        datakey = '_'.join((timing_key, key))
        job.data[datakey] = timer
    natoms_list = data['natoms']
    del data['natoms']
    for key, timer in data.items():
        timer = np.sum(timer, axis=1)
        job.doc['{:s}_{:s}_time_peratom_mean'.format(timing_key, key)] = np.mean(
            np.array(timer) / natoms_list)
        job.doc['{:s}_{:s}_time_peratom_min'.format(timing_key, key)] = np.min(
            np.array(timer) / natoms_list)
        job.doc['{:s}_{:s}_time_peratom_max'.format(timing_key, key)] = np.max(
            np.array(timer) / natoms_list)
        job.doc['{:s}_{:s}_time_peratom_std'.format(timing_key, key)] = np.std(
            np.array(timer) / natoms_list)
        thesubset = slice(job.doc.test_subset[0],
                          job.doc.test_subset[0] + job.doc.test_subset[1])
        timer_sub = timer[thesubset]
        job.doc['{:s}_{:s}_timesub_peratom_mean'.format(timing_key, key)] = np.mean(
            timer_sub / natoms_list[thesubset])
        job.doc['{:s}_{:s}_timesub_peratom_min'.format(timing_key, key)] = np.min(
            timer_sub / natoms_list[thesubset])
        job.doc['{:s}_{:s}_timesub_peratom_max'.format(timing_key, key)] = np.max(
            timer_sub / natoms_list[thesubset])
        job.doc['{:s}_{:s}_timesub_peratom_std'.format(timing_key, key)] = np.std(
            timer_sub / natoms_list[thesubset])
