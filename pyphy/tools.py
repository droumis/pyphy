import h5py, socket, json, glob, os
import pandas as pd
import numpy as np
from collections import defaultdict
import loren_frank_data_processing as lfdp
from scipy.io import loadmat


from pyphy import download

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
                    with pd.HDFStore(path, 'r') as store:
                        dtstore_cols = store.select(datatype, start=1, stop=2).reset_index().columns.tolist()
#                     print(dtstore_cols)
                except:
                    print(f'No object named {datatype} in {path}')
                    continue
                whereterms = []
                epochs = []
                ntrodes = []
                if 'epoch' in dtstore_cols:
                    # if any epoch filters can be applied to this df
                    for key,filters in datafilter.items():
                        if key in taskInfo.reset_index().columns.tolist():
                            #if the filter key can be applied to taskInfo
                            for filt in filters:
                            #for each filter condition
                                epochs.append(taskInfo.reset_index().query(filt).epoch.unique())
                    if epochs != []:
                        epochs = list(np.unique(np.concatenate(epochs)))
                    # epochs = list(set([i for x in epochs for i in x]))
                        whereterms.append('epoch == epochs') # '==' also functions as 'is in'. hdfstore query format is different than native pandas
#                     df = df.query('epoch==@epochs')
            
                if 'ntrode' in dtstore_cols:
                    
                    # if any ntrode filters can be applied to this df
                    for key,filters in datafilter.items():
                        if key in ntrodeInfo.reset_index().columns.tolist():
                            #if the filter key can be applied to ntrodeInfo
                            for filt in filters:
                            #for each filter condition
                                ntrodes.append(ntrodeInfo.reset_index().query(filt).ntrode.unique())
                    if ntrodes != []:
                        ntrodes = list(np.unique(np.concatenate(ntrodes)))
                        whereterms.append('ntrode == ntrodes') # '==' also functions as 'is in'. hdfstore query format is different than native pandas
#                     df = df.query('ntrode==@ntrodes')
                if whereterms != []:
                    print(f'loading {datatype} where {whereterms}')
    #                 import pdb; pdb.set_trace()
                    with pd.HDFStore(path, 'r') as store:
                        an_df_list.append(store.select(datatype, where=whereterms))
                else:
                    print(f'loading {datatype}')
    #                 import pdb; pdb.set_trace()
                    with pd.HDFStore(path, 'r') as store:
                        an_df_list.append(store.select(datatype))
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
    elif data_key in ['cluster_metrics']:
        data_cols = ['animal', 'day', 'ntrode', 'cluster']
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
    print(ntrode_keys)
    lfp = lfdp.get_LFPs(ntrode_keys, animal_dict)
    
    # broadcast indices from cols and reshape df to be: (an,day,ep,td):(ntrode)

    lfp_mi_list = [(str(a), int(d), int(e), int(n)) for col in lfp.columns.tolist() for (a, d, e, n) in [col.split('_')]]
    lfp_mi = pd.MultiIndex.from_tuples(lfp_mi_list, names=['animal', 'day', 'epoch', 'ntrode'])

    df = pd.DataFrame(lfp.T.values, index=lfp_mi).T #grab lfp values to avoid index merging
    df['timedelta'] = lfp.index
    df = df.set_index('timedelta').stack(level=['animal', 'day', 'epoch']).reorder_levels(['animal', 'day', 'epoch', 'timedelta'])
    return df

def _get_linearcoord_tasksegments(df, animal_dict):
    track_segments = {}
    linearcoord_arms = {}
    center_well_position = {}
    animals_days_iter = df.reset_index().groupby(['animal', 'day'])
    for (animal, day), _ in animals_days_iter:
#         print(animal, day)
        file_name = os.path.join(animal_dict[animal].directory, f'{animal}task{day:0>2d}.mat')
#         print(file_name)
        data = loadmat(file_name, variable_names=('task'))['task']
        day = data.shape[-1]
        epochs = data[0, -1][0]
        for i, epoch_data in enumerate(epochs):
            epoch = i+1
            if 'linearcoord' in epoch_data.dtype.names:
                linearcoord = epoch_data[0]['linearcoord'].item().squeeze()
                linearcoord_arms[(animal, day, epoch)] = [np.round(arm[:, :, 0], decimals=2) for arm in linearcoord]
                track_segments[(animal, day, epoch)] = [np.round(np.stack(((arm[:-1, :, 0], arm[1:, :, 0])), axis=1), decimals=2)
                                  for arm in linearcoord]
                center_well_position[(animal, day, epoch)] = track_segments[(animal, day, epoch)][0][0][0]
                track_segments[(animal, day, epoch)] = np.concatenate(track_segments[(animal, day, epoch)])
                _, unique_ind = np.unique(track_segments[(animal, day, epoch)], return_index=True, axis=0)
                track_segments[(animal, day, epoch)] = [track_segments[(animal, day, epoch)][ind] for ind in unique_ind]
    df['linearcoord'] = df.index.map(linearcoord_arms)
    df['track_segments'] = df.index.map(track_segments)
    df['center_well_position'] = df.index.map(center_well_position)
    return df

def load_from_filterframework(animal, datatype, filterframework_dir, index_keys=[]):
    if type(index_keys) != list:
        index_keys = [index_keys,]
    animal_dict = {}
    animal_dict[animal] = lfdp.Animal(directory=filterframework_dir, short_name=animal)
    
    if datatype == 'ntrodeInfo':
        out = lfdp.make_tetrode_dataframe(animal_dict)
        out['subarea'] = out['subarea'].astype(str)
        
    elif datatype == 'taskInfo':
        out = lfdp.make_epochs_dataframe(animal_dict)
    
    elif datatype == 'linearcoord' or datatype == 'task_segments':
        out = lfdp.make_epochs_dataframe(animal_dict)
        out = _get_linearcoord_tasksegments(out, animal_dict)
        
    elif datatype == 'position':
        if len(index_keys[0]) != 3:
            print('epoch_keys requred as list of (animal, day, epoch)')
            return
        position_dict_df = {}
        for (animal, day, epoch) in index_keys:
            epoch_index = (animal, day, epoch)
            position_dict_df[(animal, day, epoch)] = lfdp.get_position_dataframe(epoch_index, animal_dict)
        out = pd.concat(position_dict_df).reset_index().rename({'level_0':'animal', 'level_1':'day', 'level_2':'epoch'}, axis=1)
        out['timedelta'] = pd.TimedeltaIndex(out['time'], unit='ns')
        out.set_index(['animal', 'day', 'epoch', 'timedelta'], inplace=True)
        out['is_correct'] = out.is_correct.astype('float')
        
    elif datatype == 'lfp':
        if len(index_keys[0]) != 4:
            print('ntrode_keys requred as list of (animal, day, epoch, ntrode)')
            return
        out = _load_lfp_from_filterframework(index_keys, animal_dict)
        
    return out

