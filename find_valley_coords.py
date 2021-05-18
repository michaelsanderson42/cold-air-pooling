# https://stackoverflow.com/questions/53492508/find-plateau-in-numpy-array
# Could find plateau of each vertical strip of shuttle SRTM data.
# Then link closest plateaux, draw a line through their centres.
# Then calculate the TAF.

# Might be useful for determining coords of valleys?
# https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1002/2017WR020418

'''
    # Make a copy in coordinate format
    B = A.tocoo()
    rows = B.row
    cols = B.col
    data = B.data
    
    # rows, cols, data are numpy arrays with size 437766
    
    B.data[:6]
    # array([0.50117265, 0.49882735, 0.64171535, 0.35828465, 0.3878429 , 0.6121571 ])
    cols[:6]
    # array([ 8,  8,  9,  9, 10, 10], dtype=int32)
    rows[:6]
    # array([730,   9, 731,  10, 732,  11], dtype=int32)
    
    # data[4] + data[5] = 0.3878429 + 0.6121571 = 1.0
    # Cell 10 drains to cells 732, 11. Cell 10 has coords (10, 0)
    # Noting that elev has shape (361, 721)
    # i,j coords are:
    # i = 732 % 721 = 11    and    11 % 721 = 11
    # j = 732 // 721 = 1    and    11 // 721 = 0
    # So point (10,0) drains to (11,1) and (11,0)
'''

import os
import numpy as np
import iris
from scipy import sparse
from PyDEM.pydem.dem_processing import DEMProcessor
import matplotlib.pyplot as plt


def read_srtm(datadir, resolution):

# Filename containing the elevation data
    filename = f'srtm_douro_{resolution}.nc'
    cube = iris.load_cube(os.path.join(datadir, filename))

    return cube


def make_connectivity_matrix(cube):
    '''
    Gets the connectivity matrix from elevation data

    :param cube: An Iris cube containing the elevation data
    :returns: The connectivity matrix and the elevation data
    '''

    elev = cube.data
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points

# Create the DEM processor
    dem_proc = DEMProcessor((elev, lats, lons), dx_dy_from_file=False)

# or:
# dem_proc = DEMProcessor(filename, dx_dy_from_file=True)

# Calculate the magnitudes and directions of the slopes (needed to calculate
# the connectivity matrix)
    mag, direction = dem_proc.calc_slopes_directions()

# Calculate the connectivity matrix
    A = dem_proc.calc_connectivity_matrix()

    return A, elev


def find_coords_drainage_channel(A, x, y, nx):
    '''
    From the starting coordinates (x, y), find the coordinates of the downstream
    drainage channel
    
    :param A: Connectivity matrix, in compressed sparse column format.
    :param x: x-coordinate of starting point (as an array index)
    :param y: y-coordinate of starting point (as an array index)
    :param nx: Size of elevation array in the x-direction

    :returns: A list of coordinates (array indices) of the downstream channel
    '''

# Make a copy of the connectivity matrix in coordinate format
    B = A.tocoo()
    rows = B.row
    cols = B.col
    data = B.data

#   print('cols = ', cols)
#   print('rows = ', rows)
#   print('data = ', data)

# Cell number of starting point.
    cell_no = nx*y + x
    print('Starting coords: x=', x, ' y = ', y,' cell number = ', cell_no)
    idx = cell_no
    cell_numbers = [cell_no]
    
    # Limit the loop to 200 to avoid any chance of getting stuck in an infinite loop
    for k in list(range(200)):
        print('Loop count = ', k, 'idx = ', idx)
# Find the column where outflow data for this cell are stored. "ravel" the array to remove any dimensions equal to 1.
        a = np.argwhere(cols == idx)
        a = a.ravel()
        print('a = ', a)
# If a column is not found, have reached the end of the channel, so exit.
        if len(a) == 0:
            break
# get the outflow fractions
        fracs = data[a]
        print('Fractions = ', fracs)
        if len(a) == 1:
            idx_prev = idx
            idx, = rows[a]
            print(idx_prev, ' drains to single cell, next cell index is ', idx)
        else:
            im = fracs.argmax()
            idx = rows[a[im]]
            print('Next cell index is ', idx)
        cell_numbers.append(idx)
    
    # Convert the cell numbers to (i,j) coordinates
    i = []
    j = []
    for c in cell_numbers:
        i.append(c % nx)
        j.append(c // nx)
    
    return (i, j)


def fill_in_missing_points(x, y, elev_in):
    '''
    Fill in missing points in the valley bottom coordinates. The drainage calculations in
    pyDEM do not output the coordinates of pixels in flat areas, so need to add them in.

    :param x: list of x-coordinates of the valley pixels
    :param y: list of y-coordinates of the valley pixels
    :param elev_in: A 2D numpy array containing the elevation data

    Returns: x, y expanded with the missing points
    '''

    moore_offsets = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
    coords_out = [[x[0], y[0]]]
    npts = len(x)

# Ensure the elevation data are a numpy array of integers
    elev = np.asarray(elev_in).astype(int)

    for i in list(range(1, npts)):
        x0, x1 = [x[i-1], x[i]]
        y0, y1 = [y[i-1], y[i]]
        print('Pixels are ',x0,y0,' and ',x1,y1)
# Max distance between 2 adjacent points is sqrt(2) = 1.41
        if euclidian_distance(x0, y0, x1, y1) < 1.5:
# Points are adjacent, just copy them over
            print('No infilling needed.')
            coords_out.append([x1, y1])
        else:
# Need to fill in the gap by choosing surrounding pixel with the lowest elevation
# Check the elevations of the pixels surrounding point x0,y0 == x[i-1], y[i-1]
            print('Infilling from ', x0, y0, ' to ', x1, y1)
            print('Elevations are ', elev[y0,x0],' and ', elev[y1,x1])
            coords_old = [[x0, y0]]
            flag = True
            while flag:
                elev_surround = []
                elev0 = elev[y0, x0]
                print('Central pixel = ', x0, y0, 'elev = ', elev0)
                print('Distance to end coordinate = ', euclidian_distance(x0, y0, x1, y1))
# Get the elevations of the 8 surrounding pixels.
                for n, neighbour in enumerate(moore_offsets):
                    dx, dy = neighbour
                    xnew = x0 + dx
                    ynew = y0 + dy
                    if ynew >= 0 and xnew >= 0:
# If this pixel's coordinates equal the pixel just infilled, set the difference to a large
# value so we don't choose it.
                        if [xnew, ynew] in coords_old:
                            e = 9999.0
                        else:
                            e = elev[ynew, xnew]
                    else:
                        e = 9999.0

                    elev_surround.append(e)

# Find the lowest elevation and how many times it occurs.
                lowest_elev = min(elev_surround)
                print('elev surrounding ', elev0, ' = ',elev_surround)
                print('lowest_elev = ',lowest_elev)
                n_lowest = elev_surround.count(lowest_elev)
# NEED TO ADD A CHECK HERE FOR THE DISTANCE TO x1,y1 - IF BIGGER CHOOSE NEXT LOWEST ELEVATION
                if n_lowest == 1:
                    k = elev_surround.index(lowest_elev)
                    print('Only 1 pixel with lowest elev, k=',k, 'moore_offsets=',moore_offsets[k] )
                else:
# If 2 or more of the surrounding pixels have the same elevation choose the one nearest to
# pixel at x1, y1.
                    print('n_lowest = ',n_lowest)
                    elevs_around = np.array(elev_surround)
                    l_where_min = (elevs_around == lowest_elev)
                    print(l_where_min)
                    print(sum(l_where_min), ' pixels with lowest elev')
                    i_moore = [idx for idx in list(range(8)) if l_where_min[idx]]
                    distance = 9999.0
                    for idx in i_moore:
                        dx, dy = moore_offsets[idx]
                        x_pixel = x0 + dx
                        y_pixel = y0 + dy
                        d = euclidian_distance(x1, y1, x_pixel, y_pixel)
                        if d < distance:
                            distance = d
                            k = idx

# Calculate the coordinates of the pixel adjacent to x0,y0 with the lowest elevation
                print('k=',k)
                dx, dy = moore_offsets[k]
                xnew = x0 + dx
                ynew = y0 + dy
                print('New pixel coordinates: ',xnew,ynew, 'Elev=', elev[ynew,xnew])
                coords_out.append([xnew,ynew])
# If this pixel is adjacent to x1,y1, then finish.
                if euclidian_distance(x1, y1, xnew, ynew) < 1.5:
                    coords_out.append([x1,y1])
                    flag = False
                    print('Finished infilling, final coords are ', x1, y1,' elev=', elev[y1,x1])
                else:
                    x0 = xnew
                    y0 = ynew

                s = input('Press any key to continue:')

    npts = len(coords_out)
    xout = [coords_out[k][0] for k in list(range(1, npts))]
    yout = [coords_out[k][1] for k in list(range(1, npts))]

    return (xout, yout)


def euclidian_distance(x0, y0, x1, y1):
    '''
    Calculates the Euclidian distance between two points

    :param x0, y0: x, y values of first point
    :param x1, y1: x, y values of second point
    '''

    return np.sqrt((x1 - x0)**2 + (y1 - y0)**2)


def main():

# Select resolution of SRTM data, 30m or 90m.
    srtm_resolution = '90m'
    datadir_in = '/data/users/hadmi/SRTM/'
    datadir_out = '/home/h03/hadmi/Python/MedGOLD/cold_air_pooling/data/'

    test_run = False

    if test_run:
# For testing, use sample matrix from pyDEM paper
# x=0, y=0 are the coordinates of the pixel with the highest elevation
        x = 0
        y = 0
        data = [0.3, 0.7, 1.0, 1.0, 1.0, 0.4, 0.6, 1.0, 1.0, 1.0]
        row_ind = [3, 4, 4, 5, 4, 5, 8, 8, 7, 8]
        col_ptr = [0, 2, 3, 4, 5, 7, 8, 9, 10, 10]
        A = sparse.csc_matrix((data, row_ind, col_ptr), shape=(9, 9))
        nx = 3

    else:
        cube = read_srtm(datadir_in, srtm_resolution)

# Coordinates of top of one of the tributaries of the Douro
        x = 314
        y = 323

# Calculate the connectivity matrix
        A, elev = make_connectivity_matrix(cube)
# elev has shape (361, 721)
        nx = elev.shape[1]

    print('Read in elev data, calculated connectivity matrix')

    i, j = find_coords_drainage_channel(A, x, y, nx)
    print('Found drainage channel coords')

# Plot the elevation as a function of the position along the tributary
    npts = len(i)
    xc = list(range(npts))
    ydata = [cube.data[j[k], i[k]] for k in xc]
# reverse the elevation data so the left-hand side will be the lowest end of the valley
    ydata = ydata[::-1]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xc, ydata, marker='s', markersize=2, mfc='black', mec='black', linestyle='-', color='green')
    ax.set_title('Elevation along tributary')
    plt.show()

# Owing to 'flat' areas, and the particular method used within pyDEM, the
# valley coordinates returned from the connectivity matrix are not necessarily
# contiguous.
# Fill in the gaps by tracing a route between pixels following a line of lowest elevations
    ifill, jfill = fill_in_missing_points(i, j, cube.data)

# Write the valley coordinates as comma-separated pairs, which are easily
# read in directly to a numpy array using np.genfromtxt
#   ofilename = f'tributary_{x}_{y}_coords.dat'
#   with open(os.path.join(datadir_out, ofilename), 'w') as ofp:
#       for k in list(range(len(i))):
#           s = '{:d},{:d}\n'.format(i[k], j[k])
#           ofp.write(s)

# Plot the coordinates of the valley on the elevation data, to check all has gone well.
    n = 11
    fig = plt.figure()
    plt.subplot(2,1,1)
    plt.imshow(elev)
    plt.plot(ifill, jfill, marker='o', markersize=2, linestyle='None', color='white')
    plt.plot(i, j, marker='o', markersize=2, linestyle='None', color='red')
    plt.subplot(2,1,2)
    left_edge = min(i[:n])
    right_edge = max(i[:n])
    top_edge = max(j[:n])
    bottom_edge = min(j[:n])
    elev_region = elev[bottom_edge:top_edge+1, left_edge:right_edge+1]
    plt.imshow(elev_region)
    i_region = i[:n] - left_edge
    j_region = j[:n] - bottom_edge
    plt.plot(i_region, j_region, marker='o', markersize=2, linestyle='None', color='red')
    print([cube.data[j[k], i[k]] for k in list(range(n))])
    print(elev_region.astype(int))
    plt.show()


if __name__ == '__main__':
    main()
