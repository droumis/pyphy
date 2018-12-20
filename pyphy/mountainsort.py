import h5py, socket, json, glob, os
import pandas as pd
import numpy as np
from collections import defaultdict
from mountainlab_pytools import mdaio

from . import tools

#### MDA utils for pyphy

def get_mda_list(anim, date, ntrode, data_location):
    date = str(date)
    mda_src_dict = defaultdict(dict)

    for epdirmda in os.listdir(os.path.join(data_location, date)):
        if '.mda' in epdirmda:
            # for each nt.mda file
            for eptetmda in os.listdir(os.path.join(data_location, date, epdirmda)):
                if '.nt' in eptetmda:
                    an = eptetmda.split('_')[1]
                    ep = eptetmda.split('_')[2].split('.')[0]
                    ntr = eptetmda.split('_')[-1].split('.')[1]
                    mda_src_dict[ntr][ep] = os.path.join(data_location, date, epdirmda, eptetmda)
    mda_list = list(mda_src_dict[f'nt{ntrode}'].values())
    mda_list.sort()
    return mda_list

def get_epoch_offsets(*,dataset_dir, opts={}):

    if 'mda_list' in opts:
        # initialize with 0 (first start time)
        lengths = [0]

        for idx in range(len(opts['mda_list'])):
            ep_path=opts['mda_list'][idx]
            ep_mda=mdaio.DiskReadMda(ep_path)
            #get length of the mda (N dimension)
            samplength = ep_mda.N2()
            #add to prior sum and append
            lengths.append(samplength + lengths[(idx)])    
    
    else:
    
        prv_list = dataset_dir + '/raw.mda.prv'

        with open(prv_list, 'r') as f:
            ep_files = json.load(f)

        # initialize with 0 (first start time)
        lengths = [0]

        for idx in range(len(ep_files['files'])):
            ep_path=ep_files['files'][idx]['prv']['original_path']
            ep_mda=mdaio.DiskReadMda(ep_path)
            #get length of the mda (N dimension)
            samplength = ep_mda.N2()
            #add to prior sum and append
            lengths.append(samplength + lengths[(idx)])

    #first entries (incl 0) are starttimes; last is total time
    total_samples =lengths[-1]
    sample_offsets=lengths[0:-1]

    return sample_offsets, total_samples


def get_samples_offset(animal,date, mda_preproc_dir, \
                       ntrode_4_offsets=1, epgap_sec=300, fs=30000):
    # create mapping from mda indices to offset indices
    mda_list = get_mda_list(animal, date, ntrode_4_offsets, mda_preproc_dir)
    epoch_sample_start, total_samples = get_epoch_offsets(dataset_dir=' ', \
                                                                opts={'mda_list':mda_list})
    # add 5 minute gaps between epochs- franklab convention
    samples_offset = np.arange(1,total_samples+1)
    epoch_sample_start_offset = np.array(epoch_sample_start)
    for ioff, soff in enumerate(epoch_sample_start[1:]):
        epoch_sample_start_offset[ioff+1:] += epgap_sec * fs
        samples_offset[soff+1:] += epgap_sec * fs
    return samples_offset, epoch_sample_start_offset

def load_firings(firings_path, samples_offset, ep_samp_start_offset, animal, date, ntrode, config_path, fs=30000):
    # load firings
    firings = mdaio.readmda(firings_path)
    #create spikes_df with firing times with inter-epoch gap
    spike_cols = ['animal', 'day', 'epoch', 'ntrode', 'cluster', \
                  'timedelta', 'sampleindex']
    spikes_df = pd.DataFrame(columns=spike_cols)
    spikes_sampleindex = samples_offset[firings[1,:].astype(int)]
    spikes_timedelta = pd.TimedeltaIndex(spikes_sampleindex/fs, unit='s', name='time')
    spikes_df['sampleindex'] = spikes_sampleindex
    spikes_df['timedelta'] = spikes_timedelta
    spikes_df['cluster'] = firings[2,:].astype(int)
    spikes_df['primary_chan'] = firings[0,:].astype(int)
    #add animal, day, epoch cols:
    spikes_df['animal'] = animal
    spikes_df['day'] = tools.convert_dates_to_days(animal, date, config_path)
    spikes_df['ntrode'] = ntrode
    spikes_df['epoch'] = 0
    ep_samp_start_timedelta = pd.TimedeltaIndex([ep_s/fs  for ep_s in \
                                                 ep_samp_start_offset], \
                                                unit='s', name='time')
    for epn, ep_samp_td in enumerate(ep_samp_start_timedelta):
        spikes_df.loc[spikes_df['timedelta'] >= \
                                    ep_samp_td, 'epoch'] = epn+1
    spikes_df.set_index(['animal', 'day', 'epoch', 'timedelta', \
                                       'ntrode', 'cluster'], inplace=True)
    return spikes_df

def load_marks(marks_path, spikes_df):
    marks = mdaio.readmda(marks_path)
    if marks.shape[-1] == spikes_df.shape[0]:
        ch_cols = [f'ch{c:>02d}' for c in np.arange(marks.shape[0])]
        marks_df = pd.DataFrame(marks[:,0,:].squeeze().T, columns=ch_cols, \
                                index=spikes_df.index)
        return marks_df
    return

def load_clips(clips_path, spikes_df):
    # load clips
    clips = mdaio.readmda(clips_path)
    if clips.shape[-1] == spikes_df.shape[0]:
        clips_2d_list = [clips[int(primary_chan-1),:,i] for i, primary_chan in enumerate(spikes_df.primary_chan.values)]
        clips_2d = np.array(clips_2d_list)
        cl_cols = [f'{c:>03d}' for c in np.arange(clips_2d.shape[1])]
        clips_df = pd.DataFrame(clips_2d, columns=cl_cols, \
                                index=spikes_df.index)
        return clips_df
    return

def load_metrics(metrics_path, animal, date, ntrode, config_path):
    with open(metrics_path) as f:
        metrics = json.load(f)
    di = {}
    for cluster in metrics['clusters']:
        cmetrics = cluster['metrics']
        cmetrics['cluster'] = str(cluster['label'])
        cmetrics['multiunit'] = int('mua' in cluster['tags'])
        di[cmetrics['cluster']] = cmetrics
    metrics_df = pd.DataFrame.from_dict(di,orient='index')
    metrics_df['animal'] = animal
    metrics_df['day'] = tools.convert_dates_to_days(animal, date, config_path)
    metrics_df['ntrode'] = ntrode
    return metrics_df.set_index(['animal', 'day', 'ntrode', 'cluster'])

def get_spikes_from_mda(animal, dates, ntrodes, mda_preproc_dir, ms_firings_dir,\
                        config_path, get_marks=True, get_clips=True, get_metrics=True, epgap_sec=300, \
                        fs=30000, ntrode_4_offsets=1):
    if dates == []:
        dates = config[animal]['day2date'].values()
    elif type(dates) is str or type(dates) is int:
        dates = [str(dates),]
    spikes_dfs = []; marks_dfs = []; clips_dfs = []; metrics_dfs = []
    
    for date in dates:
        # TODO: add epoch sample start and total samples to animal config,
        # then trash the epoch-wise mda files
        samples_offset, ep_samp_start_offset=get_samples_offset(animal, date, mda_preproc_dir, \
                                                                      ntrode_4_offsets=ntrode_4_offsets,\
                                                                      epgap_sec=epgap_sec, fs=fs)
        for ntrode in ntrodes:
            print(f'animal:{animal} date:{date} ntrode:{ntrode}')
            ms4_ntdir = os.path.join(ms_firings_dir, str(date),'ms4', f'nt{ntrode}')
            firings_path = os.path.join(ms4_ntdir, 'firings_burst_merged.mda')
            if os.path.exists(firings_path):
                spikes_dfs.append(load_firings(firings_path, samples_offset, \
                                                            ep_samp_start_offset, animal,  \
                                                            date, ntrode, config_path, fs=fs))
                marks_path = os.path.join(ms4_ntdir, 'marks.mda')
                if os.path.exists(marks_path) and get_marks:
                    marks_dfs.append(load_marks(marks_path, spikes_dfs[-1]))

                clips_path = os.path.join(ms4_ntdir, 'clips.mda')
                if os.path.exists(clips_path) and get_clips:
                    clips_dfs.append(load_clips(clips_path, spikes_dfs[-1]))

                metrics_path = os.path.join(ms4_ntdir, 'metrics_merged_tagged.json')
                if os.path.exists(metrics_path) and get_metrics:
                    metrics_dfs.append(load_metrics(metrics_path, animal, date, ntrode, config_path))
    out={}
    if spikes_dfs != {}:
        out['spikes'] = pd.concat(spikes_dfs, sort=True)
    if get_marks and marks_dfs != {}:
        out['marks'] = pd.concat(marks_dfs, sort=True)
    if get_clips and clips_dfs != {}:
        out['clips'] = pd.concat(clips_dfs, sort=True)
    if get_metrics and metrics_dfs != {}:
        out['metrics'] = pd.concat(metrics_dfs, sort=True)
    
    return out
