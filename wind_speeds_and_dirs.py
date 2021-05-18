import os
import pandas as pd


def read_sogrape_data(station_name):

    dpath = '/home/h03/hadmi/observations/MED-GOLD/Portugal/'
    filename = '{}dados.xlsx'.format(st1)

# select columns containing dates, hours of day, radiation, wind speed and wind direction
    usecols = [1, 2, 5, 8, 9]
    column_names = ['date', 'hour', 'radiation', 'wind speed', 'direction']

# Read in the observations
    df = pd.read_excel(os.path.join(dpath, filename), names=column_names, skiprows=1, header=None, usecols=usecols)

# Convert the radiation and wind speed data to strings. Some entries do not have a decimal comma and are
# read as floats.
    for col_name in ['radiation', 'wind speed']:
        df[col_name] = df[col_name].astype(str)
# The radiation and wind speed data use a decimal comma. Replace by a dot and convert to a float.
        df[col_name] = [x.replace(',', '.') for x in df[col_name]]
        df[col_name] = df[col_name].astype(float)

    df['compass'] = convert_direction_to_compass(df['direction'].to_numpy())

    return df.set_index('date', inplace=True)


def convert_direction_to_compass(d):
    '''
    Converts a wind direction (in degrees) to one of the cardinal and ordinal points

    :param d: A numpy array containing wind direction(s) in degrees
    :returns: A string containing the abbreviated compass direction, 'N', 'SE', etc.
    '''

    compass_directions = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
    n_compass_dirs = len(compass_directions)
    delta = 360.0 / n_compass_dirs
    edges = [(i*delta - delta/2) for i in list(range(len(compass_directions)+1))]
    edges[0] += 360

    compass_dirs = np.full_like(d, 'N', dtype=str)
    for i in list(range(1, len(compass_directions))):
        w = (d >= edges[i]) & (d < edges[i+1])
        if sum(w) > 0:
            compass_dirs[w] = compass_directions[i]

    return compass_dirs


def read_cap_events(filename):



def ??(df_cap, df_s1):

    cap_dates = df_cap.index.values

    df_s1_cap = df_s1.loc[df_s1['dates'].isin(cap_dates)]
    df_s1_nocap = df_s1.loc[~df_s1['dates'].isin(cap_dates)]

# Analyses wind speeds on CAP days, non-CAP days, all days
    plot_num = 1
    ws = df_s1_cap['wind speed'].to_numpy()

    ax = fig.add_subplot(3, 3, 1)
    bin_max = np.ceil(np.max(ws)/2) * 2
    bins = np.arange(0, bin_max+2, 2)
    n, bin_edges, ? = ax.hist(ws, bins, width=0.8)
