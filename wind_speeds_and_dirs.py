import os
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Sun import Sun


def read_sogrape_data(station_name):
    '''
    Reads in the raw data from a Sogrape weather station

    :param station_name: The name of the station
    '''

    dpath = '/home/h03/hadmi/observations/MED-GOLD/Portugal/'
    filename = '{}dados.xlsx'.format(station_name)

# Select columns containing dates, hours of day, wind speed and wind direction
    usecols = [1, 2, 8, 9]
    column_names = ['date', 'hour', 'wind speed', 'direction']

# Read in the observations
    df = pd.read_excel(os.path.join(dpath, filename), names=column_names, skiprows=1, header=None, usecols=usecols)

    for col_name in ['wind speed']:
        df[col_name] = df[col_name].astype(str)
# The wind speed data use a decimal comma. Replace by a dot and convert to a float.
        df[col_name] = [x.replace(',', '.') for x in df[col_name]]
        df[col_name] = df[col_name].astype(float)

# Shift the dates back by half day, so each date covers the period noon to noon,
# and the night will not be divided between two days. Most (all?) CAP events should
# happen at night.
# https://stackoverflow.com/questions/62472689/pandas-shift-date-time-columns-back-one-hour

    df['dandt'] = df['date'] + ' ' + df['hour']
    date_list = [datetime.datetime.strptime(dt, "%Y-%m-%d %H:%M") for dt in df['dandt']]
    new_dandt = [d - datetime.timedelta(hours=12) for d in date_list]
    new_dates = [datetime.datetime.strftime(d, "%Y-%m-%d") for d in new_dandt]
    df['new_date'] = new_dates
    df.drop('dandt', axis=1, inplace=True)

# Add a column containing the hours as decimals.
    factors = [1, 1/60]  #  To convert hours and minutes to decimal hours.
    df['hour_dec'] = [sum(i*j for i, j in zip(map(int, hr.split(':')), factors)) for hr in df['hour']]

    return df


def calc_nighttime_means(df, coords):
    '''
    Calculates night-time mean values of variables

    :param df: Dataframe containing the variables at hourly or sub-hourly intervals
    :param coords: Dictionary containing the longitude and latitude of the weather station
    '''

    sun = Sun(coords)

    sunset_hours = []
    sunrise_hours = []
    ws_mean = []
    wd_mean = []

    df['Date'] = pd.to_datetime(df['date'])
    unique_dates = np.unique(df["Date"].dt.strftime('%Y-%m-%d'))
    for the_date in unique_dates:
        prev_date = (datetime.datetime.strptime(the_date, '%Y-%m-%d') -
            datetime.timedelta(days=1)).strftime('%Y-%m-%d')

# Calculate sunset time on the previous day and sunrise time on this day
        sunset_time = sun.getSunsetTime(prev_date)
        sunrise_time = sun.getSunriseTime(the_date)
        hr_of_sunset = sunset_time['hr'] + sunset_time['min'] / 60
        hr_of_sunrise = sunrise_time['hr'] + sunrise_time['min'] / 60
        sunset_hours.append(hr_of_sunset)
        sunrise_hours.append(hr_of_sunrise)

# Get the wind speeds and directions on this night
        df_sd = df.loc[df['new_date'] == the_date, ['hour_dec', 'wind speed', 'direction']]

# Calculate the night-time mean values
        df_tmp = df_sd.loc[(df_sd['hour_dec'] >= hr_of_sunset) | (df_sd['hour_dec'] <= hr_of_sunrise)]
        ws_mean.append(df_tmp['wind speed'].mean())
        avg_dir = calc_average_wind_direction(df_tmp['wind speed'].to_numpy(), df_tmp['direction'].to_numpy())
        wd_mean.append(avg_dir)

    df_out = pd.DataFrame({'date': unique_dates, 'sunset': sunset_hours, 'sunrise': sunrise_hours,
        'wind speed': ws_mean, 'direction': wd_mean})

    return df_out


def calc_average_wind_direction(ws, wd):
    '''
    Calculates the average wind direction using a vector approach

    :param ws: list or array containing wind speeds
    :param wd: list or array containing wind directions

    :returns: Average wind direction in degrees (range 0 - 359)
    '''
# https://math.stackexchange.com/questions/44621/calculate-average-wind-direction

# Check for NaNs
    n = len(ws)
    if np.sum(np.isnan(np.array(ws))) > n/2:
        return -999
    else:
        j = np.argwhere(~np.isnan(np.array(ws)))

# Calculate the mean u- and v- wind vectors
        u = np.mean(np.array([ws[i] * np.sin(wd[i] * np.pi/180) for i in j]))
        v = np.mean(np.array([ws[i] * np.cos(wd[i] * np.pi/180) for i in j]))
        mean_wd = np.arctan2(u, v) * 180/np.pi

        return (360 + mean_wd) % 360


def convert_direction_to_compass(wd_in, compass_directions):
    '''
    Converts a wind direction (in degrees) to one of the cardinal and ordinal points

    :param wd_in: A list containing wind direction(s) in degrees
    :returns: A string containing the abbreviated compass direction, 'N', 'SE', etc.
    '''

    d = np.array(wd_in)
    n_compass_dirs = len(compass_directions)
    delta = 360.0 / n_compass_dirs
    edges = [delta * (i-0.5) for i in list(range(len(compass_directions)+1))]
    edges[0] += 360

# Set up the compass directions array, fill with 'N' (north)
    compass_dirs = np.full_like(d, 'N', dtype=str)

    dirn_summary = {}
# Replace 'N' with other directions.
    for i in list(range(1, len(compass_directions))):
        w = (d >= edges[i]) & (d < edges[i+1])
        sw = sum(w)
        dirn_summary[compass_directions[i]] = sw
        print(i, compass_directions[i], sw)
        if sw > 0:
            compass_dirs[w] = compass_directions[i]

    w = (d == -999)
    if sum(w) > 0:
        compass_dirs[w] = 'MISSING'

    n = sum(dirn_summary.values())
    dirn_summary['N'] = len(wd_in) - sum(w) - n

    return compass_dirs, dirn_summary


def read_cap_events(datadir, filename, c_limit=1.0, mdi=99.99):
    '''
    Read in a file containing all temperature differences into a dataframe
    Return a subset containg dates when cold air pooling events occurred.
    '''

    df = pd.read_csv(os.path.join(datadir, filename))
    df = df.loc[(df['diff'] >= c_limit) & (df['hour'] != mdi)]
    df.set_index('date', inplace=True)

    return df


def count_compass_directions(wd, compass_directions):
    '''
    :param wd: List of wind directions on an 8-point compass
    '''

    counts = [sum(1 for r in wd if r == l) for l in compass_directions]
#   for l in compass_directions:
#       c = sum(1 for r in wd if r == l)
#       counts.append(c)

    return counts


def plot_wind_speeds_dirs(station, df_s1, df_cap, compass_directions):
    '''
    Analyses wind speeds and directions on CAP days, non-CAP days, all days
    Calculates the mean wind speed and direction in each night.

    :param df_s1: Dataframe containing night-time average wind speeds and directions at Sairrao 1.
    :param df_cap: Dataframe holding dates of CAP events
    '''

    cap_dates = df_cap.index.values

# Drop rows (i.e. dates) with missing data
    df_s1.dropna(axis=0, inplace=True)
    df_s1_cap = df_s1.loc[df_s1['date'].isin(cap_dates)]
    df_s1_nocap = df_s1.loc[~df_s1['date'].isin(cap_dates)]

    fig = plt.figure(figsize=(8,8))

# Histograms of wind speeds on CAP days and non-CAP days
    column_name = 'wind speed'
    ax = fig.add_subplot(2, 2, 1)
    ws = df_s1_cap[column_name].to_numpy()
    bin_max = np.ceil(np.max(ws))
    bins = np.arange(0, bin_max+1)
    n, bin_edges, _ = ax.hist(ws, bins, width=0.8, align='mid', color='grey')
    ax.set_title('CAP days')
    ax.set_xlabel('Wind Speed / m s-1')
    ax.set_ylabel('Count')
    ax.set_ylim(0, 800)

    bx = fig.add_subplot(2, 2, 2)
    ws = df_s1_nocap[column_name].to_numpy()
    bin_max = np.ceil(np.max(ws))
    bins = np.arange(0, bin_max+1)
    n, bin_edges, _ = bx.hist(ws, bins, width=0.8, align='mid', color='grey')
    bx.set_title('Non-CAP days')
    bx.set_xlabel('Wind Speed / m s-1')
    bx.set_ylabel('Count')
    bx.set_ylim(0, 800)

#   cx = fig.add_subplot(3, 3, 3)
#   ws = df_s1[column_name].to_numpy()
#   bin_max = np.ceil(np.max(ws))
#   bins = np.arange(0, bin_max+1, 1)
#   n, bin_edges, _ = cx.hist(ws, bins, width=0.8)
#   cx.set_title('All days')

# Histograms of wind directions on CAP days and non-CAP days
    column_name = 'direction'
    xticks = list(range(1, len(compass_directions)+1))
    dx = fig.add_subplot(2, 2, 3)
# Get the cardinal or ordinal points which corresponds to the mean wind directions
    wd = df_s1_cap[column_name].tolist()
    _, dirn_summary = convert_direction_to_compass(wd, compass_directions)
    sc = sum(dirn_summary.values())
    counts = [100 * dirn_summary[c] / sc for c in compass_directions]
    dx.bar(xticks, counts, width=0.8)
    dx.set_xticks(xticks)
    dx.set_xticklabels(compass_directions)
    dx.set_ylabel('Percent')
    dx.set_title('CAP days')

    ex = fig.add_subplot(2, 2, 4)
    wd = df_s1_nocap[column_name].tolist()
    _, dirn_summary = convert_direction_to_compass(wd, compass_directions)
    sc = sum(dirn_summary.values())
    counts = [100 * dirn_summary[c] / sc for c in compass_directions]
    ex.bar(xticks, counts, width=0.8)
    ex.set_xticks(xticks)
    ex.set_xticklabels(compass_directions)
    ex.set_ylabel('Percent')
    ex.set_title('Non-CAP days')

#   fx = fig.add_subplot(3, 3, 6)
#   wd = df_s1[column_name].tolist()
#   _, dirn_summary = convert_direction_to_compass(wd, compass_directions)
#   counts = [dirn_summary[c] for c in compass_directions]
#   fx.bar(xticks, counts, width=0.8)
#   fx.set_xticks(xticks)
#   fx.set_xticklabels(compass_directions)
#   fx.set_title('CAP days')

    plt.subplots_adjust(hspace=0.5)

    fpath = '/home/h03/hadmi/Python/MedGOLD/cold_air_pooling/figures/'
    filename = '{}.png'.format(station)
    plt.savefig(os.path.join(fpath, filename), dpi=150)
    plt.close()


def get_station_coords(station):

    filename = '/home/h03/hadmi/observations/MED-GOLD/Portugal/sogrape_station_metadata.txt'
    widths = [14, 8, 8, 5]

    df_metadata = pd.read_fwf(filename, skiprows=2, widths=widths)

    for i in df_metadata.index.values:
        sta = df_metadata.at[i, 'Station'].lower()
        sta = ''.join(sta.split())
#       print (i, sta)
        if station.lower() == sta:
            lon = df_metadata.at[i, 'Long']
            lat = df_metadata.at[i, 'Lat']

    return {'longitude': lon, 'latitude': lat}


def main():

    datadir = '/home/h03/hadmi/Python/MedGOLD/cold_air_pooling/data/'
    compass_directions = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
    station_name = 'SAIRRAO1'

    df_sairrao1 = read_sogrape_data(station_name)
    sairrao1_coords = get_station_coords(station_name)

    df_night = calc_nighttime_means(df_sairrao1, sairrao1_coords)

    sta_extra = ['leda', 'sairrao']
    for i, cap_file in enumerate(['LEDA3_LEDA2.csv', 'SEIXO_SAIRRAO3.csv']):
        filename = '_'.join(['nighttime_tdiffs', cap_file])
        df_cap = read_cap_events(datadir, filename)
        print(df_cap[:10])
        plot_wind_speeds_dirs('_'.join([station_name, sta_extra[i]]), df_night, df_cap, compass_directions)


if __name__ == '__main__':
    main()
