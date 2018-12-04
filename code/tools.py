def get_computer_name(verbose=False):
    import socket
    computer = socket.gethostname()
    if verbose:
        print(f'computer: {computer}')
    return computer

def read_config(config_path):
    import json
    with open(config_path, 'r') as configf:
        config = json.load(configf)
    return config

def write_config(config, config_path, mode='w'):
    import json
    with open(config_path, mode) as f:
        json.dump(config, f, indent=4, sort_keys=True)
    return

def download_data(animals_days, download_dir='', source='dropbox', fileformat='h5'):
    import download
    config = read_config('../config.json')
    if not download_dir:
        download_dir = config['download_dir']
    for animal, days in animals_days.items():
        for day in days:
            date = config[animal]['day2date'][str(day)]
            url = config[animal][date][source]
            if fileformat is 'h5':
                filename = f'{animal}_{date}.h5'
                print(f'animal:{animal} day:{day} date:{date} from:{url}/{filename} to:{download_dir}')
                download.download(url, filename, download_dir)
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
                pass
    return available_data    
