import  numpy as np
import matplotlib.pyplot as plt
from numpy.random import random
import pandas as pd


def rect_area(ra_min, ra_max, dec_min, dec_max):
    """
    Takes the boundaries of the rectangle and caluclates the area of the rectangle
    Inputs
    ------
    ra_min: minimum RA in degrees
    ra_max: maximum RA in degrees
    dec_min: minimum declination in degrees
    dec_max: maximum declination in degrees

    Returns: area of the rectangle bounded by these values in square degrees
    """

    area  = (ra_max - ra_min) * np.rad2deg(np.sin(np.deg2rad(dec_max)) -
                                           np.sin(np.deg2rad(dec_min)))
    return(area)


def plot_rect(rect1, rect2, rect3, rect4, saveloc):
    """
    Takes a list  of ra and dec boundaries for a rectangle
    and plots the boundaries and calculates are of the rectangle.
    
    Inputs
    ------
    rect1: list 
    - a list of [ra_min, ra_max, dec_min, dec_max]
    rect2, rect3, rect4: same as rect1
    saveloc: str
    - full path of location to save plot

    Output
    ------
    A plot of the four rectangles on a sphere, saved to saveloc
    """

    # KFH Plot all sides  of rectangle
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='aitoff')

    # KFH Plot all 4 boundaries, converted to radians for aitoff
    for rect in [rect1, rect2, rect3, rect4]:

        # KFH Determine area of rectangle
        area = rect_area(rect[0], rect[1], rect[2], rect[3])

        # KFH Plot boundaries of rectangle
        rect_rad = np.deg2rad(rect)
        ax.plot([rect_rad[0], rect_rad[0]], [rect_rad[2], rect_rad[3]], color='red')
        ax.plot([rect_rad[1], rect_rad[1]], [rect_rad[2], rect_rad[3]], color='red')
        ax.plot([rect_rad[0], rect_rad[1]], [rect_rad[2], rect_rad[2]], color='red')
        ax.plot([rect_rad[0], rect_rad[1]], [rect_rad[3], rect_rad[3]], color='red')
        ax.fill_between(np.arange(rect_rad[0], rect_rad[1], .01), rect_rad[2], rect_rad[3],
                        label=str('%.2f' % area), alpha=0.5)
        
    ax.grid(True)
    ax.legend(title="Area (arcsec^2)")

    fig.savefig(saveloc)
    print("Saved 4 rectangle figure to: " + saveloc)

    
def populate_rect(ra_min, ra_max, dec_min, dec_max, points=10000):
    """
    A function to find the subset of a spherically distributed set of random points
    that fall within the perimeter of a rectangle with the limits given.

    Inputs
    ------
    ra_min, ra_max, dec_min, dec_max : float
    - the boundaries of the box in degrees
    points: int
    - total number of points to populate the sphere with
    Returns
    -------
    points_in_box : df
    - All points inside the rectangle in a dataframe
    """

    # KFH First populate the entire sphere with random even points
    ra = 360.*(random(points)) - 180
    dec = (180/np.pi)*np.arcsin(1.-random(points)*2.)
    coords = pd.DataFrame({'ra': ra, 'dec' : dec})

    # KFH Subset to only points in a given rectangle
    points_in_box = coords[(ra_min < coords['ra']) &
                           (coords['ra'] < ra_max) &
                           (dec_min < coords['dec']) &
                           (coords['dec'] < dec_max)]

    return(points_in_box)


def compare_area(rect, pts):
    """
    Compare the fraction of points in a box with the fraction of arcseconds
    in area of a sphere.
    
    Inputs
    ------
    rect: list
    - a list of [ra_min, ra_max, dec_min, dec_max]

    Returns
    -------
    fraction_points : float
    - fraction of the randomly generated points inside given rect
    fraction_area : float
    - fraction of area of sphere inside given rect
    """

    # KFH Confirm that the rectangle contains the correct area relative to sphere                                                               
    points_in_box = populate_rect(rect[0], rect[1], rect[2], rect[3], points=pts)
    fraction_points = points_in_box.shape[0] / pts

    area_of_box = rect_area(rect[0], rect[1], rect[2], rect[3])
    fraction_area = area_of_box / rect_area(0, 360, -90, 90)

    return(fraction_points, fraction_area)


if __name__ == '__main__':
    from argparse import ArgumentParser

    ap = ArgumentParser(description="Plot rectangles on a sphere")
    ap.add_argument('plot_loc', help='The full filepath for output plot')
    args = ap.parse_args()
    
    # KFH Show that the hemisphere has the correct number of points
    print('Checking area of hemisphere:' + str(rect_area(0,360,0,90)))

    # Plot four boxes on an aitoff plot
    plot_rect([0, 20, 0, 15],
              [0, 20, 20, 35],
              [0, 20, 40, 55],
              [0, 20, 60, 75],
              saveloc=args.plot_loc)

    # KFH Confirm that the rectangle contains the correct area relative to sphere
    fraction_points, fraction_area  = compare_area([0,15,0,20], 100000)
    print("Fraction of points: " + str(fraction_points))
    print("Fraction of sphere: " + str(fraction_area))
