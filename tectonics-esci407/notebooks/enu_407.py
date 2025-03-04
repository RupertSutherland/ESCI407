#!/usr/bin/env python
# python3
"""
    Tools to manage orientations, stress, etc. in an ENU coordinate system.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=10,family='Arial')

__author__ = "Rupert Sutherland"

def dirVecLine(plunge,trend):
    '''
    Returns a unit direction vector from (plunge,trend) of a line.
    Coordinates are ENU (East,North,Up).

    Parameters
    ----------
    plunge : number or array
        The angle (degrees) downwards from horizontal plane.
    trend : number or array
        The angle (degrees) clockwise from North in the horizontal plane.

    Returns
    -------
    array
        Vector of unit length in direction specified.

    '''
    p = np.radians(plunge)
    t = np.radians(trend)
    # z is +ve up, but plunge is +ve down
    vz = -np.sin(p)
    vHorizontal = np.cos(p)
    vx = np.sin(t) * vHorizontal
    vy = np.cos(t) * vHorizontal
    return np.array([vx,vy,vz])

def dirVecPlane(strike,dip):
    '''
    Returns a unit direction vector perpendicular to the plane.
    The vector is upwards-pointing.
    Coordinates are ENU (East,North,Up).

    Parameters
    ----------
    strike : number or array
        Angle (degrees) clockwise from North of a horizontal line in the plane.
    dip : number or array
        The dip angle (degrees) downwards from horizontal plane.

    Returns
    -------
    array
        Unit vector perpendicular to plane(s).

    '''
    trend = strike + 90.0
    plunge = dip - 90.0
    return dirVecLine(plunge,trend)

def rotationMatrixAroundZ(angle):
    '''
    Returns a matrix for clockwise rotation by angle (degrees) around Z axis.

    Parameters
    ----------
    angle : float
        Angle (degrees) of clockwise rotation.

    Returns
    -------
    array
        3x3 rotation matrix.

    '''
    a = np.radians(angle)
    return np.array([[ np.cos(a), np.sin(a), 0],
                     [-np.sin(a), np.cos(a), 0],
                     [0, 0, 1]])
    
    
    
    
    