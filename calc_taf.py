import os
import numpy as np
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
import pandas as pd


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
    idx_centre = len(elev) // 2

# Find array indices of the local maxima
    idx_maxima, = argrelextrema(elev, np.greater_equal, order=order)
    if len(idx_maxima) < 2:
        idx_maxima = np.array([ileft, iright])

# Check the number of maxima - if more than 2, choose points closest to the centre
    if len(idx_maxima) > 2:
        for i in list(range(len(idx_maxima)-1)):
            if idx_maxima[i] < idx_centre and idx_maxima[i+1] > idx_centre:
                ileft = idx_maxima[i]
                iright = idx_maxima[i+1]

        idx_maxima = np.array([ileft, iright])

# Choose the lowest of the two maxima as the elevation of the top of the side of the valley.
# Return the elevations # equal to or below this maximum either side of the centre point.
# These latter data are the elevations of the valley cross-section.
    atol = 1.0  #  1 m tolerance
    if idx_maxima[0] < 99:
        if elev[idx_maxima[0]] < elev[idx_maxima[1]]:
            v = elev[idx_maxima[0]]
        else:
            v = elev[idx_maxima[1]]

# Find all points either side of the centre point whose elevations are smaller than the valley edge elevation (v)
# Include a point whose elevation is higher if two points straddle the valley edge elevation.
        elev_valley_xsection = [elev[idx_centre]]
        for i in list(range(1, idx_centre)):
            ileft = idx_centre - i
            if elev[ileft] <= v:
                elev_valley_xsection.insert(0, elev[ileft])
            if elev[ileft] == v:
                break
            if elev[ileft] < v and elev[ileft-1] >= v:
                elev_valley_xsection.insert(0, elev[ileft-1])
                ileft -= 1
                break
        for i in list(range(1, idx_centre)):
            iright = idx_centre + i
            if elev[iright] <= v:
                elev_valley_xsection.append(elev[iright])
            if elev[iright] == v:
                break
            if elev[iright] < v and elev[iright+1] >= v:
                elev_valley_xsection.append(elev[iright+1])
                iright += 1
                break

    else:
        elev_valley_xsection = []
        ileft = 0
        iright = 0

#   print('Elevation along perpendicular:')
#   print(list(elev))
#   print('v =', v)
#   print('Elevation along valley cross-section:')
#   print(list(elev_valley_xsection))

# Diagnostic plot
#   plot_valley_edges(cube, xpts, ypts, centre_point, v, [ileft, iright])

    return v, elev_valley_xsection


def plot_valley_edges(cube, xpts, ypts, centre_point, elev_valley_top, idx_maxima):
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

#   print('centre_point = ', centre_point)
#   print('bottom-edge: centre_point[1] - radius = ', centre_point[1] - radius)
#   print('top edge:    centre_point[1] + radius = ', centre_point[1] + radius)
#   print('lats = ', lats[centre_point[1]-radius: centre_point[1]+radius+1])

    xcoords = [lons[x] for x in xpts]
    ycoords = [lats[y] for y in ypts]

# Set up constraints and extract the block of pixels
    cube_region = cube.intersection(longitude=(left_edge, right_edge), latitude=(top_edge, bottom_edge))

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
    plt.plot([xvalues[0], xvalues[-1]], [elev_valley_top, elev_valley_top], marker='None', linestyle='--', color='grey')
    plt.plot(xmid, cube.data[centre_point[1], centre_point[0]], marker='s', linestyle='None', color='green')
    plt.plot(xvalues[idx_maxima[0]], cube.data[ypts[idx_maxima[0]], xpts[idx_maxima[0]]], marker='<', linestyle='None', color='orange')
    plt.plot(xvalues[idx_maxima[1]], cube.data[ypts[idx_maxima[1]], xpts[idx_maxima[1]]], marker='>', linestyle='None', color='orange')
    plt.title('Elevation along perpendicular')

    plt.show()


def calculate_taf(cube, xpts, ypts, centre_point):
    '''
    Calculates the topographical amplification factor (TAF) at the given location.

    elev -- Elevations at each point across the valley
    pts -- 2D list of x,y coords of points across the valley, as [2,n]
    centre_point -- Coordinates of the centre of the valley bottom

    Returns: The TAF, valley cross-sectional area, width and height
    '''

# Find the elevation of the valley top and the subset of the elevations of the transect across the valley.
    elev_top, elev_xsection = find_local_maxima(cube, xpts, ypts, centre_point)
    if len(elev_xsection) == 0:
        print('Local maxima not found for {:d},{:d}, cannot calculate TAF'.format(centre_point[0], centre_point[1]))
        W = -99.0
        H = -99.0
        A = -99.0
        taf = -99.0
    else:

# Calculate the valley width (W) and height (H); height sometimes labelled D.
        W, H = calc_valley_WH(elev_top, elev_xsection, xpts, ypts)

# Calculate the cross-sectional area of the valley (Ayz) using the Trapezium rule
        A = calc_valley_xsect_area(elev_xsection, xpts, ypts, elev_top)

# TAF = (W / Ayz) / (1/H)
        taf = (W / A) * H
        print('Width=',W,' Height=',H,' Area=',A, 'TAF=',taf)

    return (W, H, A, taf)


def euclidian_distance(xpts, ypts):
    '''
    Returns the Euclidian distance between two adjacent points

    xpts -- The x-coordinates
    ypts -- The y-coordinates
    '''

    dx = xpts[1] - xpts[0]
    dy = ypts[1] - ypts[0]

    return np.sqrt((dx*dx) + (dy*dy))


def calc_valley_WH(elev_top, elev, xpts, ypts):
    '''
    Calculate the width and height of the valley

    elev_top -- elevation of top of valley
    elev -- 1D array containing the elevation data for the valley cross-section
    xpts -- The x-coordinates of the transect across the valley
    ypts -- The y-coordinates of the transect across the valley
    '''

# Find the lowest point within the valley
    elev_lo = np.min(elev)

    valley_height = elev_top - elev_lo

# Calculate the valley width
    delta_x = euclidian_distance(xpts[:2], ypts[:2])
    npts = len(elev)

    valley_width = delta_x * (npts-1)

    return (valley_width, valley_height)


def calc_valley_xsect_area(elev, xpts, ypts, elev_top):
    '''
    Calculate the cross-sectional area of the valley (Ayz) using the Trapezium rule

    elev -- 1D array containing the elevation data for the valley cross-section
    xpts -- The x-coordinates of the transect across the valley
    ypts -- The y-coordinates of the transect across the valley
    elev_top -- Elevation of the top of the valley sides
    '''

    delta_x = euclidian_distance(xpts[:2], ypts[:2])

# Heights of the transect across the valley relative to the valley floor
    if elev[0] > elev_top:
        h = elev_top - np.array(elev[1:])
    elif elev[-1] > elev_top:
        h = elev_top - np.array(elev[:-1])
    else:
        h = elev_top - np.array(elev)

# Calculate the cross-sectional area for all points at or below elev_top
# using the Trapezium rule.
    npts = len(h)
    term0 = np.sum(h[1:npts-1])
    term1 = (h[0] + h[-1]) / 2.0
    Ayz = delta_x * (term0 + term1)

# Add on the final part of the area where the valley heights bracket the valley top
    if elev[0] > elev_top:
        x_int = delta_x * (elev_top - elev[1]) / (elev[0] - elev[1])
        Ayz += ((elev_top + elev[1])/2) * x_int
    elif elev[-1] > elev_top:
        x_int = delta_x * (elev_top - elev[-2]) / (elev[-1] - elev[-2])
        Ayz += ((elev_top + elev[-2])/2) * x_int

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

#   print('centre_point = ', centre_point)
#   print('bottom-edge: centre_point[1] - radius = ', centre_point[1] - radius)
#   print('top edge:    centre_point[1] + radius = ', centre_point[1] + radius)
#   print('lats = ', lats[centre_point[1]-radius: centre_point[1]+radius+1])

    xcoords = [lons[x] for x in xpts]
    ycoords = [lats[y] for y in ypts]

# Set up constraints and extract the block of pixels
#   lat_con = iris.Constraint(latitude = lambda y: bottom_edge <= y <= top_edge)
#   lon_con = iris.Constraint(longitude = lambda x: left_edge <= x <= right_edge)
    cube_region = cube.intersection(longitude=(left_edge, right_edge), latitude=(top_edge, bottom_edge))

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
    width = []
    height = []
    valley_area = []
    taf = []

# Distance over which to find the valley top in pixels.
# distance is in pixels, roughly 90 m per pixel
    radius = 10

    elev_filename = 'srtm_douro_90m.nc'
    cube = read_srtm_data(datadir, elev_filename)

# Coordinates of top of one of the tributaries of the Douro
    valley_head = [314, 323]

# Get the coordinates of the entire valley, as list of size (n, 2)
    valley_coords = read_valley_coordinates(datadir, valley_head)

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

#       diagnostic_plots(cube, xpts, ypts, centre_point)
# Get the valley width (W), height (H), cross-sectional area (A) and TAF (T)
        W, H, A, T = calculate_taf(cube, xpts, ypts, centre_point)

        width.append(W)
        height.append(H)
        valley_area.append(A)
        taf.append(T)

    df = pd.DataFrame({'width': width, 'height': height, 'area': valley_area, 'taf': taf})
    opath = '/home/h03/hadmi/Python/MedGOLD/cold_air_pooling/data/'
    ofilename = 'tributary_{:d}_{:d}_taf.csv'.format(valley_head[0], valley_head[1])
    df.to_csv(os.path.join(opath, ofilename), index=False)


if __name__ == '__main__':
    main()
