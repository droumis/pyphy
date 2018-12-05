import cv2
import pandas as pd
import numpy as np
### pyphy imports
import merging

def calc_spatial_firing_rate_maps(spikes, position, blur=9, pos_fs=30):
    
    spikes_w_position = merging.get_position_of_spikes(position, spikes)
    
    occnorm_FR_maps =  {}
    ep_grps = spikes_w_position.groupby(['animal', 'day', 'epoch'])
    for (animal, day, epoch), ep_spikes in ep_grps:
        # get occupancy for this epoch
        pos_ep = position.xs((animal, day, epoch), drop_level=False)
        x_edges = np.arange(pos_ep.x_position.min(),pos_ep.x_position.max(), 1) 
        y_edges = np.arange(pos_ep.y_position.min(),pos_ep.y_position.max(), 1)
        occ, _, _ = np.histogram2d(pos_ep.x_position.values, pos_ep.y_position.values, [x_edges, y_edges])
        smooth_occupancy = cv2.GaussianBlur(occ+1e-6, (blur,blur), 0) / pos_fs

        for (ntrode, cluster), cl_spikes in ep_spikes.groupby(['ntrode', 'cluster']):
            spk_map, _, _ = np.histogram2d(cl_spikes.x_position.values, cl_spikes.y_position.values, [x_edges, y_edges])
            smooth_spk_map = cv2.GaussianBlur(spk_map+1e-6, (blur,blur), 0)

            occnorm_FR_maps[(animal,day,epoch,ntrode,cluster)] = smooth_spk_map / smooth_occupancy
    return occnorm_FR_maps