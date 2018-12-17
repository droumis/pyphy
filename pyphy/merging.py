import pandas as pd

def get_position_of_spikes(position, spikes_df, tolerance_str='100ms'):
    spikes_w_pos = get_position_at_times(position, spikes_df, tolerance_str=tolerance_str)
    return spikes_w_pos.set_index(['animal', 'day', 'epoch', 'ntrode', 'cluster', 'timedelta']).sort_index()

def get_position_at_times(position, df, tolerance_str='100ms'):
    # basically a hacky merge_asof to deal with multiindex, timedelta, tolerance
    
    has_cols = [x for x in df.columns.tolist() if x in ['x_position', 'y_position', 'speed', 'head_direction']]
    if has_cols:
        print('pos cols found in df, recalculating pos')
        df = df.drop(has_cols, axis=1)
    
    tolerance = pd.Timedelta(tolerance_str)
    
    if 'timedelta' not in position.reset_index().columns.tolist():
        position.reset_index(inplace=True)
        position['timedelta'] = pd.TimedeltaIndex(position['time'], unit='ns')

    day_df_wpos_list = []
    for animal, an_df in df.groupby('animal'):
        for day, day_df in an_df.groupby('day'):
            
            pos_an_day = position.query('animal == @animal and day == @day')
            pos_an_day = pos_an_day.reset_index().set_index('timedelta').sort_index()
            
            day_df = day_df.reset_index().sort_values('timedelta')
            
            pos_indexer = pos_an_day.index.get_indexer(day_df['timedelta'], method='nearest', tolerance=tolerance)
            times_pos = pos_an_day[['x_position', 'y_position', 'head_direction', 'speed']].iloc[pos_indexer]
            day_df_wpos = pd.merge_asof(day_df, times_pos.sort_index(), on='timedelta')
            
            day_df_wpos_list.append(day_df_wpos)
    out_df = pd.concat(day_df_wpos_list, sort=True)
    return out_df