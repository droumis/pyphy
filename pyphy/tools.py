import h5py, socket, json, glob, os
import pandas as pd
import numpy as np
from collections import defaultdict
from mountainlab_pytools import mdaio
import loren_frank_data_processing as lfdp


import download

def get_computer_name(verbose=False):
    computer = socket.gethostname()
    if verbose:
        print(f'computer: {computer}')
    return computer

def make_animal_config(animal, animal_preprocessing_path):
    _dates = os.listdir(animal_preprocessing_path)
    _dates.sort()
    _days = np.arange(1,len(_dates)+1)

    _update = {}
    _update[animal] = {}
    _update[animal]['day2date'] = {}
    for dt, dy in zip(_dates, _days):
        _update[animal].update({str(dt): {'day':str(dy)}})
        _update[animal]['day2date'][str(dy)] = str(dt)
    
    return _update


def update_animal_config(animal, animal_preprocessing_path, config_path):
    config = read_config(config_path)
    if animal not in config['animals']:
        config['animals'].append(animal)
    update = make_animal_config(animal, animal_preprocessing_path)
    config.update(update)
    print(f'updating config for {animal}')
    write_config(config, config_path)
    return



def read_config(config_path):
    with open(config_path, 'r') as configf:
        config = json.load(configf)
    return config

def write_config(config, config_path, mode='w'):
    with open(config_path, mode) as f:
        json.dump(config, f, indent=4, sort_keys=True)
    return

def convert_dates_to_days(animal, dates, config_path):
    config = read_config(config_path)
    date2day =  {v:k for k,v in config[animal]['day2date'].items()}
    if type(dates) is list:
        days = [int(date2day[date]) for date in dates]
        days.sort()
    elif type(dates) is str:
        days = int(date2day[dates])
    return days

def download_data(animals_days, data_dir='', source='dropbox', fileformat='h5'):
    config = read_config('../config.json')
    if not data_dir:
        data_dir = config['data_dir']
    for animal, days in animals_days.items():
        for day in days:
            date = config[animal]['day2date'][str(day)]
            url = config[animal][date][source]
            if fileformat is 'h5':
                filename = f'{animal}_{date}.h5'
                print(f'animal:{animal} day:{day} date:{date} from:{url}/{filename} to:{data_dir}')
                download.download(url, filename, data_dir)
            else:
                print(f'fileformat {fileformat} not recognized')
                return

def get_local_data_catalogue(config_path = '../config.json'):
    config = read_config(config_path)
    subfs = [f.split('/')[-1] for f in glob.glob(config['data_dir'] + '/*')]
    available_data = {}
    for sub in subfs:
        if '.h5' in sub: #dropbox downloaded format
            animal = sub.split('_')[1]
            date = sub.split('_')[0]
            day = convert_dates_to_days(sub, date, config_path)
            available_data[animal] = day
        if sub not in ['lost+found']: # droumis local nested format
            dates = [f.split('/')[-1].split('_')[0] for f in glob.glob(config['data_dir'] + f'/{sub}/pyphy/*.h5')]
            if dates:
                days = convert_dates_to_days(sub, dates, config_path)
                available_data[sub] = days
    return available_data         

def get_remote_data_catalogue(config_path = '../config.json', source='dropbox'):
    config = read_config(config_path)
    available_data = {}
    for animal in config['animals']:
        available_data[animal] = []
        for date in config[animal].keys():
            try:
                url = config[animal][date][source]
                available_data[animal].append(config[animal][date]['day'])
            except:
                continue
    return available_data    

def get_animal_from_filepath(file):
    animal = file.split('/')[-1].split('_')[1].split('.')[0]
    return animal

def get_date_from_filepath(file):
    date = int(file.split('/')[-1].split('_')[0])
    return date

def get_day_from_filepath(file, config_path):
    animal = get_animal_from_filepath(file)
    date = get_date_from_filepath(file)
    config = read_config(config_path)
    day = config[animal][date]['day']
    return day

def get_h5_paths(animal_day_keys, config_path):
    h5_files = {}
    config = read_config(config_path)
    for animal, days in animal_day_keys.items():
        h5_files[animal] = []
        if days == []:
            # if empty grab all days
            days = config[animal]['day2date'].keys()
        for day in days:
            date = config[animal]['day2date'][str(day)]
            h5_files[animal].append(glob.glob(os.path.join(config['data_dir'], animal, 'pyphy') + f'/{date}*.h5')[0])
        h5_files[animal].sort()
    return h5_files

def read_h5_files(datafilter, datatype, config_path):
    animal_day_keys = datafilter['animals_days']
    anim_paths = get_h5_paths(animal_day_keys, config_path)      
#     print(anim_paths)
    out = []
    
    if datatype in ['ntrodeInfo', 'taskInfo']:
        # in info df, load all
        for animal, paths in anim_paths.items():
            an_df_list = []
            paths.sort()
            for path in paths:
#                 print(path)
                store = pd.HDFStore(path, 'r')
                a = store.select(datatype)
#                 a = pd.read_hdf(path, f'/{datatype}')
                an_df_list.append(a)
            out.append(pd.concat(an_df_list))
        out = pd.concat(out)

    else:
        taskInfo = read_h5_files(datafilter, datatype='taskInfo', config_path=config_path)
        ntrodeInfo = read_h5_files(datafilter, datatype='ntrodeInfo', config_path=config_path)
        
#         print(anim_paths)
        for animal, paths in anim_paths.items():
            an_df_list = []
            paths.sort()
            for path in paths:
#                 print(path)
#                 df = pd.read_hdf(path, f'/{datatype}')

                # get fields of data before loading it in
                try:
                    store = pd.HDFStore(path, 'r')
                    dtstore_cols = store.select(datatype, start=1, stop=2).reset_index().columns.tolist()
                except:
                    print(f'No object named {datatype} in {path}')
                    continue
                whereterms = []

                if 'epoch' in dtstore_cols:
                    # if any epoch filters can be applied to this df
                    epochs = []
                    for key,filters in datafilter.items():
                        if key in taskInfo.reset_index().columns.tolist():
                            #if the filter key can be applied to taskInfo
                            for filt in filters:
                            #for each filter condition
                                epochs.append(taskInfo.reset_index().query(filt).epoch.unique())
                    epochs = list(np.unique(np.concatenate(epochs)))
                    # epochs = list(set([i for x in epochs for i in x]))
                    whereterms.append('epoch == epochs') # '==' also functions as 'is in'. hdfstore query format is different than native pandas
#                     df = df.query('epoch==@epochs')
            
                if 'ntrode' in dtstore_cols:
                    # if any ntrode filters can be applied to this df
                    ntrodes = []
                    for key,filters in datafilter.items():
                        if key in ntrodeInfo.reset_index().columns.tolist():
                            #if the filter key can be applied to ntrodeInfo
                            for filt in filters:
                            #for each filter condition
                                ntrodes.append(ntrodeInfo.reset_index().query(filt).ntrode.unique())
#                     print(ntrodes)
                    ntrodes = list(np.unique(np.concatenate(ntrodes)))
#                     print(ntrodes)
                    whereterms.append('ntrode == ntrodes') # '==' also functions as 'is in'. hdfstore query format is different than native pandas
#                     df = df.query('ntrode==@ntrodes')
                an_df_list.append(store.select(datatype, where=whereterms))
            out.append(pd.concat(an_df_list)) #concat across days each animal
#             out.append(an_df_list)
        out = pd.concat(out) #concat across animals
    return out

def load_data(datafilter, datatypes=[], config_path='../config.json', source='h5'):
    if not datatypes:
        datatypes = datafilter['datatypes']
    df = {}
    for datatype in datatypes:
        if source == 'h5':
            df[datatype] = read_h5_files(datafilter, datatype, config_path)
    return df

### writing pandas df to pyphy format

def save_df_to_pyphy(df, data_key, config_path, mode='a'):
    
    config = read_config(config_path)
    if 'tetrode_number' in df.index.names:
        print('changing tetrode_number to ntrode')
        df.index.rename('ntrode', level='tetrode_number', inplace=True)

    
    if data_key in ['ntrodeInfo', 'taskInfo']:
        data_cols = df.iloc[[0]].reset_index().columns.tolist()
    elif data_key in ['position', 'dio', 'timeFilters', 'lfp', 'lfp_ripple', 'lfp_theta',]:
        data_cols = ['animal', 'day', 'epoch', 'timedelta'] 
    elif data_key in ['spikes', 'clips', 'marks']:
        data_cols = ['animal', 'day', 'epoch', 'ntrode', 'timedelta']
    print(f'data_key: {data_key}, data_cols: {data_cols}')
    
    for animal, an_df in df.groupby('animal'):
        data_dir = os.path.join(config['data_dir'], f'{animal}/pyphy/') 
        os.makedirs(data_dir, exist_ok=True)
        for day, df in an_df.groupby('day'):
            date = config[animal]['day2date'][str(day)]
            h5_file = os.path.join(data_dir, f'{date}_{animal}.h5')
            with pd.HDFStore(h5_file, mode) as f:
                print(f'saving {data_key} to {h5_file}')
#                 del f[data_key]
                f.put(data_key, df, format='table', data_columns=data_cols)

    return

#### Filter Framework utils
def _load_lfp_from_filterframework(ntrode_keys, animal_dict):
    print('loading lfp from ff')
    lfp = lfdp.get_LFPs(ntrode_keys, animal_dict)
    
    # broadcast indices from cols and reshape df to be: (an,day,ep,td):(ntrode)

    lfp_mi_list = [(str(a), int(d), int(e), int(n)) for col in lfp.columns.tolist() for (a, d, e, n) in [col.split('_')]]
    lfp_mi = pd.MultiIndex.from_tuples(lfp_mi_list, names=['animal', 'day', 'epoch', 'ntrode'])

    df = pd.DataFrame(lfp.T.values, index=lfp_mi).T #grab lfp values to avoid index merging
    df['timedelta'] = lfp.index
    df = df.set_index('timedelta').stack(level=['animal', 'day', 'epoch']).reorder_levels(['animal', 'day', 'epoch', 'timedelta'])
    return df

def load_from_filterframework(animal, datatype, filterframework_dir, epoch_keys=[], ntrode_keys=[]):
    animal_dict = {}
    animal_dict[animal] = lfdp.Animal(directory=filterframework_dir, short_name=animal)
    
    if datatype == 'ntrodeInfo':
        out = lfdp.make_tetrode_dataframe(animal_dict)
        out['subarea'] = out['subarea'].astype(str)
        
    elif datatype == 'taskInfo':
        out = lfdp.make_epochs_dataframe(animal_dict)
        
    elif datatype == 'position':
        if epoch_keys == []:
            print('ntrode_keys requred')
            return
        position_dict_df = {}
        for (animal, day, epoch) in epoch_keys:
            epoch_index = (animal, day, epoch)
            position_dict_df[(animal, day, epoch)] = lfdp.position._get_pos_dataframe(epoch_index, animal_dict)
        out = pd.concat(position_dict_df).reset_index().rename({'level_0':'animal', 'level_1':'day', 'level_2':'epoch'}, axis=1)
        out['timedelta'] = pd.TimedeltaIndex(out['time'], unit='ns')
        out.set_index(['animal', 'day', 'epoch', 'timedelta'], inplace=True)
        
    elif datatype == 'lfp':
        if ntrode_keys == []:
            print('ntrode_keys requred for LFP loading')
            return
        out = _load_lfp_from_filterframework(ntrode_keys, animal_dict)
        
    return out



#### MDA utils

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





def get_spikes_from_mda(animal, dates, ntrodes, mda_preproc_dir, ms_firings_dir, config_path, epgap_sec = 300, fs = 30000, cluster_tag_filters={'mua':False}):
    
    if dates == []:
        dates = config[animal]['day2date'].values()
        
    spikes_df = defaultdict(lambda: defaultdict(dict))
    
    for date in dates:
        # create mapping from mda indices to offset indices
        # TODO: add epoch sample start and total samples to animal config, then trash the epoch-wise mda files
        ntrode_4_offsets = 1
        mda_list = get_mda_list(animal, date, ntrode_4_offsets, mda_preproc_dir)
        epoch_sample_start, total_samples = get_epoch_offsets(dataset_dir=' ', opts={'mda_list':mda_list})
        # add 5 minute gaps between epochs- franklab convention
        samples_offset = np.arange(1,total_samples+1)
        epoch_sample_start_offset = np.array(epoch_sample_start)
        for ioff, soff in enumerate(epoch_sample_start[1:]):
            epoch_sample_start_offset[ioff+1:] += epgap_sec * fs
            samples_offset[soff+1:] += epgap_sec * fs

        for ntrode in ntrodes:
            print(f'animal:{animal} date:{date} ntrode:{ntrode}')
            
            # load firings
            firings_path = os.path.join(ms_firings_dir, str(date), 'ms4', f'nt{ntrode}', 'firings_burst_merged.mda')
            firings = mdaio.readmda(firings_path)
            clusters = set(firings[2,:].astype(int))
            
            # load cluster metrics 
            metrics_path = os.path.join(ms_firings_dir, str(date), 'ms4', f'nt{ntrode}', 'metrics_merged_tagged.json')
            with open(metrics_path) as f:
                metrics = json.load(f)
            cluster_tags = {metrics['clusters'][c]['label']: metrics['clusters'][c]['tags'] for c in np.arange(0,len(metrics['clusters']))}

            firings_bool = []
            for key,val in cluster_tag_filters.items():
                print(f'applying cluster filter: {key} : {val}')
                if val == False:
                    keep_clusters = [c for c,t in cluster_tags.items() if key not in t]
                    if keep_clusters != []:
                        firings_bool = [s in keep_clusters for s in firings[2,:]]
                if val == True:
                    keep_clusters = [c for c,t in cluster_tags.items() if key in t]
                    if keep_clusters != []:
                        firings_bool = [s in keep_clusters for s in firings[2,:]]
                
            print(f'keeping clusters: {keep_clusters}')
            print(f'{len(firings_bool)} firings across clusters')
            
            if firings_bool != []:
                spikes_df[date][ntrode] = pd.DataFrame(columns=['animal', 'day', 'epoch', 'ntrode', 'cluster', 'timedelta', 'sampleindex'])
                spikes_mda_inds = firings[1,firings_bool].astype(int)
                spikes_sampleindex = samples_offset[spikes_mda_inds]
                spikes_timedelta = pd.TimedeltaIndex(spikes_sampleindex/fs, unit='s', name='time')
                spikes_df[date][ntrode]['sampleindex'] = spikes_sampleindex
                spikes_df[date][ntrode]['timedelta'] = spikes_timedelta
                spikes_df[date][ntrode]['cluster'] = firings[2,firings_bool].astype(int)
                
                #add animal, day, epoch cols
                spikes_df[date][ntrode]['animal'] = animal
                spikes_df[date][ntrode]['day'] = convert_dates_to_days(animal, date, config_path)
                spikes_df[date][ntrode]['ntrode'] = ntrode
                spikes_df[date][ntrode]['epoch'] = 0
                ep_samp_start_timedelta = pd.TimedeltaIndex([ep_s/fs  for ep_s in epoch_sample_start_offset], unit='s', name='time')
                for epn, ep_samp_td in enumerate(ep_samp_start_timedelta):
                    spikes_df[date][ntrode].loc[spikes_df[date][ntrode]['timedelta'] >= ep_samp_td, 'epoch'] = epn+1
#                 spikes_df[day][ntrode].set_index(['animal', 'day', 'epoch', 'ntrode', 'cluster', 'spikes_timedelta'])
    # spikes = pd.concat(spikes_df[day], ignore_index=True)
    
    df_list = []
    for date in spikes_df.keys():
        for ntrode, df in spikes_df[date].items():
            df_list.append(df)
            
    df_out = pd.concat(df_list, sort=True).set_index(['animal', 'day', 'epoch', 'ntrode', 'cluster'])
        
    return df_out
