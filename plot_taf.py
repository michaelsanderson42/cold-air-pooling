import os
import pandas as pd
import matplotlib.pyplot as plt


def read_taf(path, valley_head):

#   df = pd.DataFrame({'width': width, 'height': height, 'area': valley_area, 'taf': taf})
    filename = 'tributary_{:d}_{:d}_taf.csv'.format(valley_head[0], valley_head[1])

    return pd.read_csv(os.path.join(path, filename), header=0)


def plot_taf(df):

    df['WdivA'] = df['width'] / df['area']

    data_names = list(df.columns)
    xvalues = df.index.values

    for col in data_names:
        yvalues = df[col].to_numpy()
        fig = plt.figure()
        plt.plot(xvalues, yvalues, marker='o', linestyle='None', color='green')
        plt.title(col)

        plt.show()


def main():

    datadir = '/home/h03/hadmi/Python/MedGOLD/cold_air_pooling/data/'
    valley_head = [314, 323]

    df = read_taf(datadir, valley_head)

    plot_taf(df)


if __name__ == '__main__':
    main()
