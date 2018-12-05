import h5py, socket, json, glob, os
import pandas as pd
import numpy as np
from collections import defaultdict

import download

def get_computer_name(verbose=False):
    computer = socket.gethostname()
    if verbose:
        print(f'computer: {computer}')
    return computer

def read_config(config_path):
    with open(config_path, 'r') as configf:
        config = json.load(configf)
    return config

def write_config(config, config_path, mode='w'):
    with open(config_path, mode) as f:
        json.dump(config, f, indent=4, sort_keys=True)
    return

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

def get_data_catalogue(config_path = '../config.json', source='dropbox'):
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
        if not days:
            # if empty grab all days
            days = config[animal]['day2date'].keys()
        for day in days:
            date = config[animal]['day2date'][str(day)]
            h5_files[animal] = glob.glob(os.path.join(config['data_dir'], animal, 'pyphy') + f'/{date}*.h5')
        h5_files[animal].sort()
            
    return h5_files

def read_h5_files(datafilter, datatype, config_path):
    animal_day_keys = datafilter['animals_days']
    out = []
    if datatype in ['ntrodeInfo', 'taskInfo']:
        anim_paths = get_h5_paths(animal_day_keys, config_path)        
        for animal, paths in anim_paths.items():
            an_df_list = []
            paths.sort()
            for path in paths:
                an_df_list.append(pd.read_hdf(path, f'/{datatype}'))
            out.append(pd.concat(an_df_list))
        out = pd.concat(out)

    else:
        taskInfo = read_h5_files(datafilter, datatype='taskInfo', config_path=config_path)
        ntrodeInfo = read_h5_files(datafilter, datatype='ntrodeInfo', config_path=config_path)
        
        anim_paths = get_h5_paths(animal_day_keys, config_path)        
        for animal, paths in anim_paths.items():
            an_df_list = []
            paths.sort()
            for path in paths:
                df = pd.read_hdf(path, f'/{datatype}')
                
                if 'epoch' in df.reset_index().columns.tolist():
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
                    df = df.query('epoch==@epochs')
                    
                if 'ntrode' in df.reset_index().columns.tolist():
                    # if any ntrode filters can be applied to this df
                    ntrodes = []
                    for key,filters in datafilter.items():
                        if key in ntrodeInfo.reset_index().columns.tolist():
                            #if the filter key can be applied to ntrodeInfo
                            for filt in filters:
                            #for each filter condition
                                ntrodes.append(ntrodeInfo.reset_index().query(filt).tetrode_number.unique())
                    print(ntrodes)
                    ntrodes = list(np.unique(np.concatenate(ntrodes)))
                    print(ntrodes)
                    df = df.query('ntrode==@ntrodes')
                
                an_df_list.append(df)
            out.append(pd.concat(an_df_list))
        out = pd.concat(out)
    return out

def load_data(datafilter, datatypes=[], config_path='../config.json', source='h5'):
    if not datatypes:
        datatypes = datafilter['datatypes']
    df = {}
    for datatype in datatypes:
        if source == 'h5':
            df[datatype] = read_h5_files(datafilter, datatype, config_path)
    return df