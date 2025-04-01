#!/usr/bin/env python
# python3
"""
    Tools to manage orientations, stress, etc. in ENU (XYZ) coordinates.
"""
import numpy as np
import matplotlib.pyplot as plt

__author__ = "Rupert Sutherland"

def vecDotProduct(arrayOfColumnVectors1,arrayOfColumnVectors2):
    '''
    Given 2 arrays of column vectors of the same dimensions.
    each with structure something like:
    [[Ux0,Ux1,...Uxn],  and [[Vx0,Vx1,...Vxn],
     [Uy0,Uy1,...Uyn],       [Vy0,Vy1,...Vyn],
     [Uz0,Uz1,...Uzn]]       [Vz0,Vz1,...Vzn]]
    
    Returns the dot products
    -------
    [Ux0Vx0+Uy0Vy0+Uz0Vz0, ... UxnVxn+UynVyn+UznVzn]

    '''
    return np.sum(arrayOfColumnVectors1*arrayOfColumnVectors2, axis=0)

def vecMagnitude(arrayOfColumnVectors):
    '''
    arrayOfColumnVectors has structure something like:
    [[x0,x1,x2,...xn],
     [y0,y1,y2,...yn],
     [z0,z1,z2,...zn]]
    '''
    return np.sqrt(vecDotProduct(arrayOfColumnVectors,arrayOfColumnVectors))

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
    vHorizontal = np.cos(p)
    vx = np.sin(t) * vHorizontal
    vy = np.cos(t) * vHorizontal
    # z is +ve up, but plunge is +ve down
    vz = -np.sin(p)
    v = np.array([vx,vy,vz])
    # force the array to have 3 rows, even if it is a single vector
    if len(v.shape) == 1: 
        v.shape = (3,1)
    return v

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
    trend = np.array(strike) + 90.0
    plunge = np.array(dip) - 90.0
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
                     [0,          0,         1]])
      
def stressTensorDiagonal(sigma_E, sigma_N, sigma_U):
    '''
    Returns a stress tensor (3x3 matrix) with principal stresses (in order ENU)
    as consecutive elements in a diagonal matrix. 
    '''
    return np.array([[sigma_E, 0,       0      ],
                     [0,       sigma_N, 0      ],
                     [0,       0,       sigma_U]])
 
def makeStressTensor(SHmin,SHmax,SV,azimuth_SHmax):
    '''
    Constructs a stress tensor, given the principal stress magnitudes
    and azimuth of the maximum horizontal principal stress.

    Returns
    -------
    array
        3x3 rotation matrix.

    '''
    stressTensor_principal = stressTensorDiagonal(SHmin,SHmax,SV)
    rotationMatrix = rotationMatrixAroundZ(azimuth_SHmax)
    return rotationMatrix @ stressTensor_principal @ rotationMatrix.T

def stressOnPlane(stressTensor,vecNormal):
    # stress 3x3 matrix multiplied by 3xn matrix of normal vectors
    vecTraction = stressTensor @ vecNormal
    # dot product - works for an array of vectors
    normalStress = vecDotProduct(vecNormal,vecTraction)
    shearTraction = vecTraction - np.multiply(normalStress,vecNormal)
    shearStress = vecMagnitude(shearTraction) 
    return normalStress,shearStress

def strengthCoulomb(cohesion,coefficientFriction,normalStress):
    '''
    Coulomb failure criterion
    '''
    return cohesion + coefficientFriction * normalStress

    
    
