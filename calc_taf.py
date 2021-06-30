import os
import numpy as np
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema


def read_srtm_data(datadir, filename):
    '''
    Reads the SRTM data.

    Returns an Iris cube
    '''

# Filename containing the elevation data
    cube = iris.load_cube(os.path.join(datadir, filename))

    print(cube)

    return cube


def read_valley_coordinates(datadir, valley_head):
    '''
    Reads the valley coordinates (array indices)
    :param datadir: Directory where the coordinate files are located
    :param valley_head: List of length 2 containing the x,y coordinates of the head of the valley

    Returns a n x 2 list, where n is the number of points
    '''

    path = '/home/h03/hadmi/Python/MedGOLD/cold_air_pooling/data/'
    filename = 'tributary_{:d}_{:d}_coords.dat'.format(valley_head[0], valley_head[1])

    coords = np.genfromtxt(os.path.join(path, filename), delimiter=',').astype(int)

    return coords.tolist()


def calc_perpendicular_line(coord, slope):
    '''
    Calculates the perpendicular to the line defined by slope m and intercept c
    at point coord
    :param coord: Point through which the perpendicular must pass
    :param slope: slope of the straight line; a value of 10 means the line is vertical
     (i.e. the slope would be infinite)

    :returns: slope of the perpendicular to the straight line
    '''

    if slope == 0.0:
        perp_slope = 10.0
    elif slope == 10.0:
        perp_slope = 0.0
    else:
        perp_slope = -1.0 / slope

    return perp_slope


def points_along_line(slope, centre_point, radius):
    '''
    Returns the coordinates of pixels along the straight line through centre_point
    within the given radius.

    :param slope: slope of the straight line
    :param centre_point: Point at the centre of the line
    :param radius: Extend the line 'radius' pixels either side of the centre point
 
    :returns: A list of pixel coordinates along the line
    '''

    npts = 2*radius + 1

    if slope == 10.0:
# The line is vertical
        x = [centre_point[0]] * npts
        y = [centre_point[1] + d for d in list(range(-radius, radius+1))]
    elif slope == 0.0:
# The line is horizontal
        x = [centre_point[0] + d for d in list(range(-radius, radius+1))]
        y = [centre_point[1]] * npts
    else:
# The line is diagonal
        x = [centre_point[0] + d for d in list(range(-radius, radius+1))]
        y = [centre_point[1] + (d*int(slope)) for d in list(range(-radius, radius+1))]

    return (x, y)


def find_local_maxima(cube, xpts, ypts, centre_point, order=3):
    '''
    Finds the coordinates of the maxima in the DEM data along the given line,
    each side of the centre point.
    These maxima are taken to be the heights of the valley sides.
    Returns the elevations of the valley cross-section, excluding points outside
    of the valley.

    cube -- Iris cube containing the DEM data for the entire domain.
    xpts, ypts -- Array indices of a line across the valley (the valley cross-section).
    centre_point -- Array indices of the centre of the line, corresponding to the valley floor.
    order -- Number of points on each side of a given point to use for the comparison.
    '''

    elev = np.array([cube.data[y,x] for x,y in zip(xpts, ypts)])
    ileft = 99
    iright = 99

# Find array indices of the local maxima
    idx_maxima, = argrelextrema(elev, np.greater_equal, order=order)
    if len(idx_maxima) < 2:
        idx_maxima = np.array([ileft, iright])

# Check the number of maxima - if more than 2, choose points closest to the centre
    if len(idx_maxima) > 2:
        idx_centre = len(elev) // 2
        for i in list(range(len(idx_maxima)-1)):
            if idx_maxima[i] < idx_centre and idx_maxima[i+1] > idx_centre:
                ileft = idx_maxima[i]
                iright = idx_maxima[i+1]

        idx_maxima = np.array([ileft, iright])

# Choose the lowest of the two maxima as the valley height, and return the elevations
# equal to or below this maximum.
    if idx_maxima[0] < 99:
        if elev[idx_maxima[0]] < elev[idx_maxima[1]]:
            v = elev[idx_maxima[0]]
        else:
            v = elev[idx_maxima[1]]

        elev_valley_xsection = elev[(elev <= v)]
    else:
        elev_valley_xsection = []

    return elev_valley_xsection


def calculate_taf(cube, xpts, ypts, centre_point):
    '''
    Calculates the topographical amplification factor (TAF) at the given location.

    elev -- Elevations at each point across the valley
    pts -- 2D list of x,y coords of points across the valley, as [2,n]
    centre_point -- Coordinates of the centre of the valley bottom

    Returns: The TAF
    '''

# Find the indices of the local maxima along this line, which will be the tops of the valley sides
    elev_xsection = find_local_maxima(cube, xpts, ypts, centre_point)
    if len(elev_xsection) == 0:
        print('Local maxima not found for {:d},{:d}, cannot calculate TAF'.format(centre_point[0], centre_point[1]))
        taf = -99.0
    else:

# Calculate the valley width (W) and height (H); height sometimes labelled D.
        W = calc_valley_width(elev_xsection, xpts, ypts)
        H, elev_top = calc_valley_height(elev_xsection)

# Calculate the cross-sectional area of the valley (Ayz) using the Trapezium rule
        A = calc_valley_xsect_area(elev_xsection, xpts, ypts, elev_top)

# TAF = (W / Ayz) / (1/H)
        taf = (W / A) * H

    return taf


def euclidian_distance(xpts, ypts):
    '''
    Returns the Euclidian distance between two adjacent points

    xpts -- The x-coordinates
    ypts -- The y-coordinates
    '''

    dx = xpts[1] - xpts[0]
    dy = ypts[1] - ypts[0]

    return np.sqrt((dx*dx) + (dy*dy))


def calc_valley_width(elev, xpts, ypts):
    '''
    Calculate the width of the valley

    elev -- 1D array containing the elevation data for the valley cross-section
    xpts -- The x-coordinates of the transect across the valley
    ypts -- The y-coordinates of the transect across the valley
    '''

    delta_x = euclidian_distance(xpts[:2], ypts[:2])
    npts = len(elev)

    return delta_x * (npts-1)


def calc_valley_height(elev):
    '''
    Calculate the height of the valley

    elev -- 1D array containing the elevation data for the valley cross-section
    '''

# Find the lowest and highest points within the valley
    elev_lo = np.min(elev)
    elev_top = np.max(elev)

# Valley height, often referred to as H or D
    valley_height = elev_top - elev_lo + 1

    return valley_height, elev_top


def calc_valley_xsect_area(elev, xpts, ypts, elev_top):
    '''
    Calculate the cross-sectional area of the valley (Ayz) using the Trapezium rule

    elev -- 1D array containing the elevation data for the valley cross-section
    xpts -- The x-coordinates of the transect across the valley
    ypts -- The y-coordinates of the transect across the valley
    elev_top -- Elevation of the top of the valley sides
    '''

    npts = len(elev)
    delta_x = euclidian_distance(xpts[:2], ypts[:2])

# Heights of the transect across the valley relative to the valley floor
    h = elev_top - np.array(elev)
    term0 = np.sum(h[1:npts])
    term1 = (h[0] + h[-1]) / 2.0
    Ayz = delta_x * (term0 + term1)

    return Ayz


def diagnostic_plots(cube, xpts, ypts, centre_point):
    '''
    Creates diagnostic plots to check valley cross-sections

    cube -- Iris cube containing DEM for entire Douro valley domain
    xpts -- x-coordinates (in array indices) of line across valley
    ypts -- y-coordinates (in array indices) of line across valley
    centre_point -- Coordinates of point in centre of valley (should be the lowest point)
    '''

# Create plots to check has worked
    fig = plt.figure()

# Find the edges of a block of pixels centred on the valley mid-point.
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    radius = len(xpts) // 2
    left_edge = lons[centre_point[0] - radius]
    right_edge = lons[centre_point[0] + radius]
    top_edge = lats[centre_point[1] + radius]
    bottom_edge = lats[centre_point[1] - radius]

    print('centre_point = ', centre_point)
    print('bottom-edge: centre_point[1] - radius = ', centre_point[1] - radius)
    print('top edge:    centre_point[1] + radius = ', centre_point[1] + radius)
    print('lats = ', lats[centre_point[1]-radius: centre_point[1]+radius+1])

    xcoords = [lons[x] for x in xpts]
    ycoords = [lats[y] for y in ypts]

# Set up constraints and extract the block of pixels
#   lat_con = iris.Constraint(latitude = lambda y: bottom_edge <= y <= top_edge)
#   lon_con = iris.Constraint(longitude = lambda x: left_edge <= x <= right_edge)
    cube_region = cube.intersection(longitude=(left_edge, right_edge), latitude=(top_edge, bottom_edge))

#   xcoords = [lons[x] for x in xpts]
#   ycoords = [lats[y] for y in ypts]

# Coordinates of perpendicular line (latitudes and longitudes)
    xcentre = lons[centre_point[0]]
    ycentre = lats[centre_point[1]]

# First subplot showing elevations around centre point and the perpendicular line
    plt.subplot(2,1,1)
    iplt.pcolormesh(cube_region)
    plt.plot(xcoords, ycoords, marker='None', markersize=2, linestyle='-', color='white')
    plt.plot(xcentre, ycentre, marker='o', markersize=2, linestyle='None', color='red')

# Second subplot showing the elevations along the perpendicular line and at the centre point
    plt.subplot(2,1,2)
    xvalues = list(range(len(xpts)))
    xmid = xvalues[len(xvalues) // 2]
    yvalues = [cube.data[py, px] for px, py in zip(xpts, ypts)]
    print(yvalues)
    plt.plot(xvalues, yvalues, 'bo')
    plt.plot(xmid, cube.data[centre_point[1], centre_point[0]], marker='s', linestyle='None', color='green')
    plt.title('Elevation along perpendicular')

    plt.show()
        

def main():

    datadir = '/data/users/hadmi/SRTM/'

# Distance over which to find the valley top in pixels.
# distance is in pixels, roughly 90 m per pixel
    radius = 10

    elev_filename = 'srtm_douro_90m.nc'
    cube = read_srtm_data(datadir, elev_filename)

# Coordinates of top of one of the tributaries of the Douro
    valley_head = [314, 323]

# Get the coordinates of the entire valley, as list of size (n, 2)
    valley_coords = read_valley_coordinates(datadir, valley_head)

# Add a copy of the first and last coords to the ends
#   valley_coords.insert(0, valley_coords[0])
#   valley_coords.extend(valley_coords[-1])

# Loop over the valley coordinates
    for n in list(range(1, len(valley_coords))):
        x0, y0 = valley_coords[n-1]
        x1, y1 = valley_coords[n]

# Following Lunquist et al. (2009) find the slope between adjacent points
        if x0 == x1:
# slope between points is infinite as line is vertical; set to 10.
            slope = 10.0
        elif y0 == y1:
# slope between points is zero as line is horizontal
            slope = 0.0
        else:
            slope = (y1 - y0) / (x1 - x0)

        print(x0, y0, x1, y1, slope)

# Find the slope of the perpendicular line through the centre point
        centre_point = [x1, y1]
        m_p = calc_perpendicular_line(centre_point, slope)

# Get the x,y coordinates along the perpendicular line
        xpts, ypts = points_along_line(m_p, centre_point, radius)
        print('Points along perpendicular line')
        print([(x,y) for x, y in zip(xpts, ypts)])
        print('Centre point:', centre_point)

# Plot the elevation data around the valley bottom and the perpendicular line
# for a visual check.
        diagnostic_plots(cube, xpts, ypts, centre_point)
        taf = calculate_taf(cube, xpts, ypts, centre_point)


if __name__ == '__main__':
    main()
