from imagej_tiff_meta import TiffFile
import glob
from scipy.spatial import distance

def get_pixels_per_micron(file):
    t = TiffFile(file)
    info = t.pages[0].info().split('\n')
    res_entry = info[[idx for idx, s in enumerate(info) if 'x_resolution' in s][0]]
    res_clean = res_entry.split('(')[-1].strip(')').split(',')
    num,den = [int(s) for s in res_clean]
    t.close()
    return num/den

def get_annotations(file):
    pix_per_mic = get_pixels_per_micron(file)
    t = TiffFile(file)
    overlays = {}
    for roi in t.pages[0].imagej_tags['parsed_overlays']:
        overlays[roi['name']] = {}
        overlays[roi['name']]['x'] = roi['left']/pix_per_mic
        overlays[roi['name']]['y'] = roi['top']/pix_per_mic
    t.close()
    return overlays

def get_ntrode_distances(image_dir):
    distances = {}
    files = glob.glob(image_dir)
    for file in files:
        overlays = get_annotations(file)
        d_MEC = [overlays['d_MEC']['x'], overlays['d_MEC']['y']]
        for name, xy  in overlays.items():
            if 'nt_' in name:
                xyl = [xy['x'], xy['y']]
                distances[int(name.split('_')[-1])] = distance.euclidean(xyl, d_MEC)
    return distances