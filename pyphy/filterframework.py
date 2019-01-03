
#standard lib imports
import h5py, socket, json, glob, os, sys
import pandas as pd
import numpy as np
from scipy.io import loadmat

### lfdp imports
import loren_frank_data_processing as lfdp

### pyphy imports
import pyphy


""" Filter Framework utils

"""

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
        
    if 'tetrode_number' in out.index.names:
        print('renaming "tetrode_number" to "ntrode" in multiindex')
        out.index.rename('ntrode', level='tetrode_number', inplace=True)

    return out

