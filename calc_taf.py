import os
import numpy as np
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
from scipy.stats import linregress


def read_srtm_data(datadir, filename):
    '''
    Reads the SRTM data.

    Returns an Iris cube
    '''

# Filename containing the elevation data
    cube = iris.load_cube(os.path.join(datadir, filename))

    return cube


def read_valley_coordinates(datadir, valley_head):
    '''
    Reads the valley coordinates (array indices)
    :param datadir: Directory where the coordinate files are located
    :param valley_head: List of length 2 containing the x,y coordinates of the head of the valley

    Returns a n x 2 list, where n is the number of points
    '''

    filename = 'tributary_{:d}_{:d}_coords.dat'.format(valley_head[0], valley_head[1])

    coords = np.genfromtxt(os.path.join(datadir, filename), delimiter=',').astype(int)

    return coords.tolist()


def calc_perpendicular_line(coord, m, c):
    '''
    Calculates the perpendicular to the line defined by slope m and intercept c
    at point coord
    :param coord: Point through which the perpendicular must pass
    :param m, c: slope and intercept of a straight line

    :returns: slope and intercept of the perpendicular to the straight line
    '''

    perp_slope = -1 / m
    perp_int = coord[1] - perp_slope * coord[0]

    return (perp_slope, perp_int)


def points_along_line(slope, intercept, centre_point, radius):
    '''
    Returns the coordinates of pixels along the straight line through centre_point
    within the given distance

    :param slope, intercept: slope and intercept of the straight line
    :param centre_point: Point at the centre of the line
    :param radius: Extend the line 'distance' pixels either side of the centre point
 
    :returns: A list of pixel coordinates along the line
    '''

    coords = []
    alpha = np.arctan(slope)
    for d in list(range(-radius, radius+1)):
        x = np.int(0.5 + d * np.cos(alpha) + centre_point[0])
        y = np.int(0.5 + d * np.sin(alpha) + centre_point[1])
        if [x,y] not in coords:
            coords.append([x, y])

    return coords


def calculate_taf_pt(cube, pts, centre_point):
    '''
    Calculates the TAF at the given location

    cube -- Iris cube containing the elevation data
    pts -- 2D list of x,y coords across the valley, as [2,n]
    centre_point -- Coordinates of the centre of the valley bottom

    Returns: The TAF
    '''

    npts = len(pts[0])
    elev_valley_bottom = cube.data[centre_point[0], centre_point[1]]


def diagnostic_plots(cube, xcoords, ycoords, line0, pts):

    centre_point = [xcoords[1], ycoords[1]]

# Create plots to check has worked
    fig = plt.figure()
# Show all the elevation data
    plt.subplot(2,1,1)
    im = plt.imshow(cube.data)

# Plot the straight line fitted to 3 points in the valley
    y_sl0 = line0[1] + line0[0] * (centre_point[0]-10)
    y_sl1 = line0[1] + line0[0] * (centre_point[0]+10)
    plt.plot([centre_point[0]-10, centre_point[0], centre_point[0]+10], [y_sl0, centre_point[1], y_sl1],
        marker='o', markersize='2', linestyle='-', color='black')

# Plot the perpendicular line
    perp_line_x = [pts[i][0] for i in list(range(len(pts)))]
    perp_line_y = [pts[i][1] for i in list(range(len(pts)))]
    plt.plot([perp_line_x[0], perp_line_x[-1]], [perp_line_y[0], perp_line_y[-1]],
        marker='^', markersize='2', linestyle='-', color='red')

    plt.subplot(2,1,2)
    xvalues = list(range(len(perp_line_x)))
    yvalues = [cube.data[px, py] for px, py in zip(perp_line_x, perp_line_y)]

    plt.plot(xvalues, yvalues, 'b-')

    plt.show()
        

def main():

    datadir = '/home/h03/hadmi/Python/MedGOLD/cold_air_pooling/data/'

# distance is in pixels, roughly 90 m per pixel
    distance = 30

    elev_filename = 'srtm_douro.nc'
    cube = read_srtm_data(datadir, elev_filename)

# Coordinates of top of one of the tributaries of the Douro
    valley_head = [314, 323]

# Get the coordinates of the entire valley, as list of size (n, 2)
    valley_coords = read_valley_coordinates(datadir, valley_head)

# Add a copy of the first and last coords to the ends
    valley_coords.insert(0, valley_coords[0])
    valley_coords.extend(valley_coords[-1])

# Loop over the valley coordinates
    for n in list(range(1, len(valley_coords)-1)):
        xcoords = [valley_coords[i][0] for i in [n-1, n, n+1]]
        ycoords = [valley_coords[i][1] for i in [n-1, n, n+1]]
        centre_point = [xcoords[1], ycoords[1]]

# Fit a straight line through the three points
        result = linregress(xcoords, ycoords)

# Find the slope and intercept of the perpendicular line through the centre point
        m_p, c_p = calc_perpendicular_line(centre_point, result[0], result[1])

# Get the x,y coordinates along the perpendicular line
        pts = points_along_line(m_p, c_p, centre_point, distance)
        print('Points along perpendicular line')
        print(pts)

        diagnostic_plots(cube, xcoords, ycoords, result, pts)

# Calculate the TAF at this location
#       taf = calculate_taf_pt(cube, pts, centre_point)


if __name__ == '__main__':
    main()
