"""
Creates a meshgrid of points in [-10,10]x[-10,10] and solves the
minimum-perimeter problem at each point.
"""

import argparse
import sys

import numpy as np
import pandas as pd
import sympy as sp

from numba import njit

sys.path.append('/home/amit/WorkSpace/JHU/minimum-perimeter')
from isoperimetric.polar.mobile.areajacobian import totalarea, centroid
from isoperimetric.polar.mobile.lengthjacobian import totallength
from isoperimetric.polar.mobile.generatecode import (generatefunction,
        generateEFGfunc)


def calc_obj_func(summaryfile, S, k=1e9, v=2.0):
    """
    Calculate the total objective function value
    from parts.

    Parameter:
    ----------
    summaryfile: a CSV file containing Perimeter, AreaConstraint and ViscousTerm
    S: the system size used in the simulation
    k: the penalty term used in the simulation
    v: the viscosity used in the simulation

    Output:
    -------
    objfunc: objective function values for each row of summaryfile.
    """
    df = pd.read_csv(summaryfile)
    obj = df.Perimeter/S + 0.5*(k/S**4)*df.AreaConstraint**2 + \
            df.ViscousTerm/S**2
    return obj.to_numpy()


zexpr = '(111/20)*sin(pi*x/10)*cos(pi*y/10) + ' + \
        '(2/5)*(cos(pi*x/2)*sin(pi*y/2))'
symmap={'sin': sp.sin, 'cos':sp.cos, 'pi':sp.pi}
ffft = generatefunction(zexpr, symmap)
efg = generateEFGfunc(zexpr, symmap)
S = 2.0
k = 1e9

inputdf = pd.read_csv('Input.txt')

viscosity = 2.0

for i, row in inputdf.iterrows():
    pointid = int(row.PointId)
    area = row.Area
    a = row.X
    b = row.Y

    resultfile = 'path-{0}.csv'.format(pointid)
    resultsarr = np.loadtxt(resultfile, delimiter=',')

    outputfile = 'summary-{0}.csv'.format(pointid)
    with open(outputfile, 'w') as out:
        out.write('Xc,Yc,Perimeter,AreaConstraint,ViscousTerm,ObjFunc\n')
        numpts = resultsarr.shape[1] - 2

        # Starting position
        R0 = 0.8*np.sqrt(area/np.pi)*np.ones(numpts)
        x, y = centroid(R0, a, b, efg)
        perimeter = totallength(R0, a, b, efg)
        areaconst = totalarea(R0, a, b, efg) - area
        viscousterm = 0.0
        obj_fun = perimeter/S + 0.5*(k/S**4)*areaconst**2

        out.write('{0:.4f},{1:.4f},{2:.4f},'.format(x, y, perimeter) + \
                '{0:.4g},{1:.4g},{2:.4g}\n'.format(areaconst, viscousterm,
                    obj_fun))

        # Later positions
        for j, result in enumerate(resultsarr):
            R = resultsarr[j, :-2]
            a, b = resultsarr[j, -2:]
            x, y = centroid(R, a, b, efg)
            perimeter = totallength(R, a, b, efg)
            areaconst = totalarea(R, a, b, efg) - area
            viscousterm = 0.5*viscosity*((R - R0)**2).mean()
            obj_fun = perimeter/S + 0.5*(k/S**4)*areaconst**2 + \
                    viscousterm/S**2
            out.write('{0:.4f},{1:.4f},{2:.4f},'.format(x, y, perimeter) + \
                    '{0:.4g},{1:.4g},{2:.4g}\n'.format(areaconst, viscousterm,
                        obj_fun))
            R0 = R.copy()
