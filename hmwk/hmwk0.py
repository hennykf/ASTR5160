import random
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt


def mb(m, b):
    """
    Given an m and b, generates 10 random numbers between 0-10, finds their y values according to the equation y=mx+b, and adds some random gaussian-distributed error centered at 0 with a standard dev of 0.5 to the caluclated y values.
    Inputs:
    m: a value that will be the slope of the equation y=mx+b
    b: a value that will be the y-intercept of the linear equation y=mx+b
    Returns:
    x values, y values, and the error added to the y values 
    """
    # KFH First generate 10 floating point numbers and find y vals for them
    xvals = np.random.uniform(0, 10, size=10)
    yvals = linfunc(xvals, m, b)

    # KFH Add a gaussian-distributed error to y values
    scatter = np.random.normal(0,.5,10)
    scattered_yvals = yvals + scatter
    return(xvals, scattered_yvals, scatter)

def linfunc(x, m, b):
    # KFH Creates a linear function
    return(m*x + b)

def linfunc_plot(m, b):
    # KFH Accepts values of m and b, randomly generates x data and y data with error accordingly, finds a best fit for the generated x and y data, and plots the original line, the data, and the best fit line.

    # KFH Create random data with error, and find a best fit for it following the form mx+b
    x, y, yerr  = mb(m,b)
    popt, pcov = scipy.optimize.curve_fit(linfunc, x, y)

    # KFH Plot original data, original line with m and b, and line with fitted m and b.
    plt.clf()
    plt.plot(x, linfunc(x, m, b), '--', label = "Original line", color = "blue")
    plt.errorbar(x, y, yerr=yerr, label =  "Original data", fmt="o", ms=3, color = "green")
    plt.plot(x, linfunc(x, popt[0], popt[1]), ':', label =  "Linear fit of data", color="red")
    plt.legend()
    plt.xlabel("x value")
    plt.ylabel("y value")
    plt.title("Linear fitting  of randomly generated data")
    plt.savefig('/d/www/kianah/public_html/hmwk0_line_fitting.png')
    print("You will find the plot at: http://faraday.uwyo.edu/~kianah/hmwk0_line_fitting.png")


m = float(input("Please enter a value for slope m: "))
b = float(input("Please enter a value for y intercept b: "))
linfunc_plot(m, b)
# Note: I'm not sure if you meant that we should have it accept values at the command line, but I've interpretted it as meaning that if you load this python file, you can run the function linfunc_plot() and have everything happen.
