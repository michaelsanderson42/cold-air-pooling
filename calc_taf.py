import os
import numpy as np
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy.stats import linregress
from scipy import interpolate
from Haversine import haversine
from utils import rotate_points
import pandas as pd


def read_srtm_data(datadir, resolution='90'):
    '''
    Reads the SRTM data.

    Returns an Iris cube
    '''

# Filename containing the elevation data
    elev_filename = f'srtm_douro_{resolution}m.nc'
    return iris.load_cube(os.path.join(datadir, elev_filename))


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


def calc_perpendicular_line(slope):
    '''
    Calculates the slope of a line perpendicular to that defined by 'slope'.

    :param slope: slope of the straight line; a value of 999 means the line is vertical
     (i.e. the slope would be infinite)

    :returns: slope of the perpendicular to the straight line
    '''

    if slope == 0.0:
        perp_slope = 999.0
    elif slope == 999.0:
        perp_slope = 0.0
    else:
        perp_slope = -1.0 / slope

    return perp_slope


def points_along_line(centre_point, radius, slope):
    '''
    Returns the coordinates of pixels along the straight line through centre_point
    within the given radius.

    :param centre_point: Point at the centre of the line
    :param radius: Extend the line a 'radius' number of pixels either side of the centre point
    :param slope: slope of the straight line fitted to points along the valley floor
 
    :returns: A list of pixel coordinates along the line, shape(n, 2)
    '''

# Find the slope of a perpendicular line (m_p), which will define a transect across the valley.
    m_p = calc_perpendicular_line(slope)

# Calculate the coordinates of points along a transect of length 2*radius
    if m_p == 0:
# Transect is horizontal
        pline_coords = [[centre_point[0]+r, centre_point[1]] for r in list(range(-radius, radius+1))]
    elif m_p == 999.0:
# Transect is vertical
        pline_coords = [[centre_point[0], centre_point[1]+r] for r in list(range(-radius, radius+1))]
    else:
# Transect is diagonal
        pline_coords = []
        alpha = np.arctan(m_p)
        for r in list(range(-radius, radius+1)):
            h = r * np.sin(alpha)
            b = r * np.cos(alpha)
            pline_coords.append([centre_point[0] + b, centre_point[1] + h])

    return pline_coords


def find_local_maxima(cube, z_transect, pts, centre_point, order=3):
    '''
    Finds the coordinates of the maxima in the DEM data along the given line,
    each side of the centre point.
    These maxima are taken to be the heights of the valley sides.
    Returns the elevations of the valley cross-section, excluding points outside
    of the valley.

    cube -- Iris cube containing DEM data for the Douro region (only needed for plotting).
    z_transect -- Elevations of all points along the transect across the valley.
    pts -- Array indices of the transect across the valley.
    centre_point -- Array indices of the centre of the line, corresponding to the valley floor.
    order -- Number of points on each side of a given point to use for the comparison.
    '''

    ileft = 99
    iright = 99
    idx_centre = len(z_transect) // 2

# Find array indices of the local maxima
    idx_maxima, = argrelextrema(z_transect, np.greater_equal, order=order)
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
# Return the elevations equal to or below this maximum either side of the centre point.
# These latter data are the elevations of the valley cross-section.
    atol = 1.0  #  1 m tolerance
    if idx_maxima[0] < 99:
# v is the elevation of the valley top
        if z_transect[idx_maxima[0]] < z_transect[idx_maxima[1]]:
            v = z_transect[idx_maxima[0]]
        else:
            v = z_transect[idx_maxima[1]]

# Find all points either side of the centre point whose elevations are smaller than the valley edge elevation (v)
# Include a point whose elevation is higher if two points straddle the valley edge elevation.
        elev_valley_xsection = [z_transect[idx_centre]]
        for i in list(range(1, idx_centre)):
            ileft = idx_centre - i
            if z_transect[ileft] <= v:
                elev_valley_xsection.insert(0, z_transect[ileft])
            if z_transect[ileft] == v:
                break
            if z_transect[ileft] < v and z_transect[ileft-1] >= v:
                elev_valley_xsection.insert(0, z_transect[ileft-1])
                ileft -= 1
                break
        for i in list(range(1, idx_centre)):
            iright = idx_centre + i
            if z_transect[iright] <= v:
                elev_valley_xsection.append(z_transect[iright])
            if z_transect[iright] == v:
                break
            if z_transect[iright] < v and z_transect[iright+1] >= v:
                elev_valley_xsection.append(z_transect[iright+1])
                iright += 1
                break
    else:
        v = 0.0
        elev_valley_xsection = []

#   print('Elevation along perpendicular:')
#   print(list(elev))
#   print('v =', v)
#   print('Elevation along valley cross-section:')
#   print(list(elev_valley_xsection))

# Diagnostic plot
    if abs(iright-ileft+1) >= 10:
        plot_valley_edges(cube, pts, centre_point, v, z_transect, [ileft, iright])

    return v, [ileft, iright], elev_valley_xsection


def get_coords_of_indices(cube, xpts, ypts):
    '''
    Finds the longitudes and latitudes corresponding to the array indices.
    The array indices may be floats.

    cube -- Iris cube containing DEM for entire Douro valley domain
    xpts -- x-coordinates (in array indices) of line across valley
    ypts -- y-coordinates (in array indices) of line across valley
    '''

    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points

# Coordinates of line (latitudes and longitudes)
    dlon = lons[1] - lons[0]
    dlat = lats[1] - lats[0]
    xcoords = [lons[int(x)] + dlon*(x - int(x)) for x in xpts]
    ycoords = [lats[int(y)] + dlat*(y - int(y)) for y in ypts]

    return (xcoords, ycoords)


def plot_valley_edges(cube, pts, centre_point, elev_valley_top, z_transect, idx_maxima):
    '''
    Creates diagnostic plots to check valley cross-sections

    cube -- Iris cube containing DEM for entire Douro valley domain
    pts -- x- and y-coordinates (in array indices) of line across valley
    centre_point -- Coordinates of point in centre of valley (should be the lowest point)
    '''

# Create plots to check has worked
    fig = plt.figure()

# Find the edges of a block of pixels centred on the valley mid-point.
    radius = len(z_transect) // 2
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    left_edge = lons[centre_point[0] - radius]
    right_edge = lons[centre_point[0] + radius]
    top_edge = lats[centre_point[1] + radius]
    bottom_edge = lats[centre_point[1] - radius]

#   print('centre_point = ', centre_point)
#   print('bottom-edge: centre_point[1] - radius = ', centre_point[1] - radius)
#   print('top edge:    centre_point[1] + radius = ', centre_point[1] + radius)
#   print('lats = ', lats[centre_point[1]-radius: centre_point[1]+radius+1])

# Coordinates of perpendicular line (latitudes and longitudes)
    xpts, ypts = np.array(pts).T
    xcoords, ycoords = get_coords_of_indices(cube, xpts, ypts)

# Set up constraints and extract the block of pixels
    cube_region = cube.intersection(longitude=(left_edge, right_edge), latitude=(top_edge, bottom_edge))

    xcentre = lons[centre_point[0]]
    ycentre = lats[centre_point[1]]

# First subplot showing elevations around centre point and the perpendicular line
    plt.subplot(2,1,1)
    iplt.pcolormesh(cube_region)
    plt.plot(xcoords, ycoords, marker='None', linestyle='-', color='white')
    plt.plot(xcentre, ycentre, marker='o', markersize=2, linestyle='None', color='red')

# Second subplot showing the elevations along the perpendicular line and at the centre point
    plt.subplot(2,1,2)
    xvalues = list(range(len(z_transect)))
    xmid = xvalues[radius]

    plt.plot(xvalues, z_transect, 'bo')
    plt.plot([xvalues[0], xvalues[-1]], [elev_valley_top, elev_valley_top], marker='None', linestyle='--', color='grey')
    plt.plot(xmid, cube.data[centre_point[1], centre_point[0]], marker='s', linestyle='None', color='green')
    plt.plot(xvalues[idx_maxima[0]], z_transect[idx_maxima[0]], marker='<', linestyle='None', color='orange')
    plt.plot(xvalues[idx_maxima[1]], z_transect[idx_maxima[1]], marker='<', linestyle='None', color='orange')
    plt.title('Elevation along perpendicular')

    plt.show()


def calculate_taf(cube, centre_point, slope, radius, interpol_dem_2d):
    '''
    Calculates the topographical amplification factor (TAF) at the given location.

    z_transect -- Elevations at each point on the transect of the valley
    pts -- x, y indices of points on the transect of the valley, a 2D list of size (n,2)
    centre_point -- Indices of the centre of the transect, the valley bottom

    Returns: The TAF, valley cross-sectional area, width and height
    '''

# order is the number of points used to identify local maxima.
# Order = 2 means 2 points either side of a given point will be tested.
    order = 2
    z = cube.data.data

# Get the x,y coordinates of points along a transect across the valley, which will be
# perpendicular to the line following the valley floor.
    pts = points_along_line(centre_point, radius, slope)
    print('Coords of points along perpendicular line')
    print(pts)
    print('Centre point:', centre_point)

# Get the elevations along the transect, using bilinear interpolation if needed.
    if slope == 0.0 or slope == 999.0:
# Transect is vertical or horizontal, no interpolation needed.
        z_transect = np.array([z[p[1],p[0]] for p in pts])
    else:
# Set a list of points ordered (y,x), as the DEM data are ordered (lat, lon)
        pts_r = [(p[1],p[0]) for p in pts]
        z_transect = interpol_dem_2d(pts_r)

    print('Elevations along transect across valley:')
    print(z_transect)

#   diagnostic_plots(cube, pts, centre_point, radius, z_transect)

# Find the elevation of the valley top and the subset of the elevations that describe the valley only.
    elev_top, idx_edges, elev_xsection = find_local_maxima(cube, z_transect, pts, centre_point, order=order)
    print(idx_edges, elev_top)

# If valley width is too wide, the transect may follow the valley instead of crossing it,
# or follow the line of a tributary.
    n_attempts = 1
    while (abs(idx_edges[1] - idx_edges[0]) > 10) and (n_attempts < 7):
# Rotate the transect ...
        angle = 22.5 * n_attempts
        new_pts = rotate_points(pts, centre_point, angle)

# ... and recalculate the valley width
        pts_r = [(p[1],p[0]) for p in new_pts]
        z_transect = interpol_dem_2d(pts_r)
        elev_top, idx_edges, elev_xsection = find_local_maxima(cube, z_transect, new_pts, centre_point, order=order)
        n_attempts += 1
        print('After attempt ', n_attempts)
        print(z_transect)
        print(idx_edges, elev_top)

    if n_attempts > 1:
        pts = new_pts[:]

    if len(elev_xsection) == 0:
        print('Local maxima not found for {:d},{:d}, cannot calculate TAF'.format(centre_point[0], centre_point[1]))
        W = -99.0
        H = -99.0
        A = -99.0
        taf = -99.0
    else:

# Calculate the valley width (W), height (H) and cross-sectional area(A) [Ayz]
        W, H, A = calc_valley_WHA(cube, elev_xsection, pts, centre_point, elev_top)

# TAF = (W / Ayz) / (1/H)
        taf = (W / A) * H
        print('Width=',W,' Height=',H,' Area=',A, 'TAF=',taf)

    return (W, H, A, taf)


def calc_valley_WHA(cube, elev, pts, centre_point, elev_top):
    '''
    Calculate the width, height and area of the valley

    elev -- 1D array containing the elevation data for the valley cross-section
            only (not necessarily the whole transect)
    pts -- The coordinates of the transect across the valley
    centre_point -- The x- and y-coordinates of the centre of the valley
    elev_top -- elevation of top of valley
    '''

# Find the lowest point within the valley
    elev_lo = np.min(elev)
    valley_height = elev_top - elev_lo

# Get the longitudes and latitudes of the array indices in xpts, ypts
    xpts, ypts = np.array(pts).T
    tlons, tlats = get_coords_of_indices(cube, xpts, ypts)

# Calculate the distance between two points located at the centre of the transect in metres
# using the Haversine formula.
    j = len(tlons) // 2
    t_distance = haversine([tlats[j], tlons[j]], [tlats[j+1], tlons[j+1]]) * 1000.0

# Calculate the valley width and cross-sectional area (Ayz)
    valley_width, Ayz = calc_valley_xsect_area(elev, t_distance, elev_top)

    return (valley_width, valley_height, Ayz)


def calc_valley_xsect_area(elev, t_distance, elev_top):
    '''
    Calculate the width (W) and cross-sectional area of the valley (Ayz) using the Trapezium rule.
    The distance between pairs of points of the transect is assumed to be constant, although
    in reality it will vary slightly.

    elev -- 1D array containing the elevation data for the valley cross-section
    t_distance -- Distance between two points located at the centre of the transect.
    elev_top -- Elevation of the top of the valley sides
    '''

# Heights of the transect across the valley relative to the valley floor
    if elev[0] > elev_top:
        h = elev_top - np.array(elev[1:])
    elif elev[-1] > elev_top:
        h = elev_top - np.array(elev[:-1])
    else:
        h = elev_top - np.array(elev)

    valley_width = t_distance * (len(h)-1)

# Calculate the cross-sectional area of the valley for all points at or below elev_top
# using the Trapezium rule.
    npts = len(h)
    term0 = np.sum(h[1:npts-1])
    term1 = (h[0] + h[-1]) / 2.0
    Ayz = t_distance * (term0 + term1)

# Add on the final part of the area where the valley heights bracket the valley top,
# using the simple formula for area of a triangle. x_int is the base of the triangle and
# the height is elev_top - elev[1]  or  elev_top - elev[-2]
    if elev[0] > elev_top:
        x_int = t_distance * (elev_top - elev[1]) / (elev[0] - elev[1])
        Ayz += (elev_top - elev[1]) * x_int / 2
    elif elev[-1] > elev_top:
        x_int = t_distance * (elev_top - elev[-2]) / (elev[-1] - elev[-2])
        Ayz += (elev_top - elev[-2]) * x_int / 2
    else:
        x_int = 0.0

# Add the base of the triangle above to the valley width
    return valley_width+x_int, Ayz


def diagnostic_plots(cube, pts, centre_point, radius, z_transect):
    '''
    Creates diagnostic plots to check valley cross-sections

    cube -- Iris cube containing DEM for entire Douro valley domain
    pts -- x- and y-coordinates (in array indices) of line across valley
    centre_point -- Coordinates of point in centre of valley (should be the lowest point)
    radius --
    z_transect -- Elevations of points along the transect across the valley.
    '''

# Create plots to check has worked
    fig = plt.figure()

# Find the edges of a block of pixels centred on the valley mid-point.
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points

    left_edge = lons[centre_point[0] - radius]
    right_edge = lons[centre_point[0] + radius]
    top_edge = lats[centre_point[1] + radius]
    bottom_edge = lats[centre_point[1] - radius]

#   print('centre_point = ', centre_point)
#   print('bottom-edge: centre_point[1] - radius = ', centre_point[1] - radius)
#   print('top edge:    centre_point[1] + radius = ', centre_point[1] + radius)
#   print('lats = ', lats[centre_point[1]-radius: centre_point[1]+radius+1])

#   xcoords = [lons[x] for x in xpts]
#   ycoords = [lats[y] for y in ypts]

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
#   plt.plot(xcoords, ycoords, marker='None', markersize=2, linestyle='-', color='white')
    plt.plot(xcentre, ycentre, marker='o', markersize=2, linestyle='None', color='red')

# Second subplot showing the elevations along the perpendicular line and at the centre point
    plt.subplot(2,1,2)
    xvalues = list(range(len(z_transect)))
    xmid = xvalues[len(xvalues) // 2]
#   print(z_transect)
    plt.plot(xvalues, z_transect, 'bo')
    print(xmid, z_transect[xmid])
    plt.plot(xmid, z_transect[xmid], marker='s', linestyle='None', color='green')
    plt.title('Elevation along transect')

    plt.show()


def calc_slope(valley_coords, n, k=5):
    '''
    Calculates the slope of a straight line through the given points

    valley_coords -- Array indices of the valley floor, a list of size(n, 2)
    n -- Index of a point. The straight line will be fitted to k points
         centred on point n
    k -- Number of points to fit the line to (must be 2 or larger)

    Returns the slope of the fitted line
    '''

    if k == 2:
        k2 = 0
    else:
        k2 = k // 2

    lo = n-k2
    if lo < 0:
        lo = 0
    hi = n+k2
    if hi > len(valley_coords):
        hi = len(valley_coords)

    x = [valley_coords[j][0] for j in list(range(lo, hi))]
    y = [valley_coords[j][1] for j in list(range(lo, hi))]

# Check for horizontal / vertical points
    if all(i == x[0] for i in x):
# All points have same x-value, so line is vertical
        slope = 999.0
    elif all(j == y[0] for j in y):
# All points have same y-value, so line is horizontal
        slope = 0.0
    else:
        result = linregress(x, y)
        slope = result[0]

    return slope
        

def main():

    datadir = '/data/users/hadmi/SRTM/'
    width = []
    height = []
    valley_area = []
    taf = []

# Distance over which to find the valley top in pixels.
# distance is in pixels, roughly 90 m per pixel
    radius = 10

# Read in the DEM (elevation) data
    cube = read_srtm_data(datadir, resolution='90')

# Set up a bilinear interpolation function for the DEM data
    z = cube.data.data
    ny, nx = cube.shape
    xi = list(range(nx))
    yi = list(range(ny))
    interpol_dem_2d = interpolate.RegularGridInterpolator((yi, xi), z, method='linear')

# Coordinates of top of one of the tributaries of the Douro
    valley_head = [314, 323]

# Get the coordinates of the entire valley, as list of size (n, 2)
    valley_coords = read_valley_coordinates(datadir, valley_head)

# A straight line will be fitted through 'k' points centred at 'centre-point'
    k = 5

# Loop over the valley coordinates
    n0 = len(valley_coords)-10
    for n in list(range(1, len(valley_coords))):
# Find the slope of a line through 'k' points centred at 'centre-point'
        print('n=', n)
        centre_point = [valley_coords[n][0], valley_coords[n][1]]
        slope = calc_slope(valley_coords, n, k=k)
        print('Slope through ',k,' points is ', slope)

# Get the valley width (W), height (H), cross-sectional area (A) and TAF (T)
        W, H, A, T = calculate_taf(cube, centre_point, slope, radius, interpol_dem_2d)

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
