# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:22:02 2025

@author: agarz
"""
import numpy as np
import gauss_legendre
from scipy.integrate import quad
import matplotlib.pyplot as plt
import sys

pi = np.pi

def gl_quad(fun, a, b):
    """
    performs integration of a function by 6-point Gauss-Legendre quadrature

    Parameters:
    fun: function that takes as argument a one-dimensional ndarray with N
    elements and outputs another one-dimensional ndarray with N elements
    a: lower integration limit
    b: upper integration limit

    Returns:
    Scalar with value of integral
    """
    midpoint = 0.5*(a+b)
    half_length = 0.5*(b-a)
    x = midpoint + half_length*gauss_legendre.z
    return half_length*np.sum(fun(x)*gauss_legendre.w)

def fun(x):
    return np.tan(x)

# testing gl_quad
a, b = 0.0, 0.5
I1 = gl_quad(fun, a, b)
I2, err = quad(fun, a, b)
print(I1)
print(I2)


def greens_function(x, y, x0, y0):
    """
    Computes free-space Green's function G
    and its partial derivatives dG/dx, dG/dy
    to be used in the boundary element method

    Parameters:
    x: one-dimensional ndarray of N elements
    y: one-dimensional ndarray of N elements
    x0: scalar
    y0: scalar

    Returns:
    G: one-dimensional ndarray of N elements
    Gx: one-dimensional ndarray of N elements with dG/dx
    Gy: one-dimensional ndarray of N elements with dG/dy
    """
    
    dx = x - x0
    dy = y - y0
    rs = dx**2 + dy**2
    r = np.sqrt(rs)
    G = -np.log(r)/(2*pi)

    denom = 2*pi*rs
    Gx = - dx/denom
    Gy = - dy/denom

    return G, Gx, Gy

def greens_parametrized(z, point_a, point_b, point_0):
    """
    Computes Green's function and its normal derivative,
    parametrized along a straight element by the local
    coordinate z in [-1, 1]

    Parameters:
    z: one-dimensional ndarray of N elements
    xa: scalar
    xb: scalar
    ya: scalar
    yb: scalar

    Returns:
    G: one-dimensional ndarray of N elements
    Gx: one-dimensional ndarray of N elements
    Gy: one-dimensional ndarray of N elements
    """
    xa, ya = point_a
    xb, yb = point_b
    x0, y0 = point_0
    
    xm = 0.5*(xa+xb) # midpoint
    dx = xb-xa
    dxh = 0.5*dx
    
    ym = 0.5*(ya+yb) # midpoint
    dy = yb-ya
    dyh = 0.5*dy

    x = xm + dxh*z
    y = ym + dyh*z

    return greens_function(x, y, x0, y0)

def gl_quad_greens(point_a, point_b, point_0):
    """
    performs integration by 6-point Gauss-Legendre quadrature of Green's function G and
    its partial derivatives dG/dx, dG/dy

    Parameters:
    fun: function that takes a one dimensional ndarray with N elements and
    returns 3 ndarrays with N elements with the values of G, dG/dx,
    and dG/dy at each of the N input values.
    a: lower integration limit
    b: upper integration limit
    
    Returns:
    intG: scalar
    intGx: scalar
    intGy: scalar
    """
    xa, ya = point_a
    xb, yb = point_b

    h = 0.5*np.sqrt((xb-xa)**2 + (yb-ya)**2) # half element length
    G, Gx, Gy = greens_parametrized(gauss_legendre.z,
                                    point_a, point_b, point_0)
    intG = h*np.sum(G*gauss_legendre.w)
    intGx = h*np.sum(Gx*gauss_legendre.w)
    intGy = h*np.sum(Gy*gauss_legendre.w)

    return intG, intGx, intGy
    
def influence_coefficients(point_a, point_b, point_0, is_singular):
    """
    influence coefficients computed according to Eqs. (3.2.11) and (3.2.12)
    in C. Pozrikidis book
    """
    xa, ya = point_a
    xb, yb = point_b
    
    # components of unit vector normal to the element (nx, ny)
    dx = xb-xa
    dy = yb-ya
    dr = np.sqrt(dx**2 + dy**2)
    nx = dy/dr
    ny = -dx/dr

    if not is_singular:
        alpha, intGx, intGy = gl_quad_greens(point_a, point_b, point_0)
        beta = nx*intGx + ny*intGy
    else: # is singular
        alpha = -dr*(np.log(0.5*dr)-1)/(2*pi)
        beta = 0
        
    return alpha, beta

    
def is_even(n: int) -> bool:
    """Return True if n is even, False otherwise."""
    return n % 2 == 0

def elements_line(point_a, point_b, N, ratio):
    """
    Divides a straight-line segment from (xa, xb) to (ya, yb) into N
    elements of varying length such that
    ratio = (length of element N/2)/(length of element 1)
          = (length of element N/2+1)/(length of elemnt N)
    """
    xa, ya = point_a
    xb, yb = point_b

    if not is_even(N):
        print(f"{N} is not even. Stopping execution.")
        sys.exit(1)

    Nh = int(N/2) # N must be even
    alpha = ratio**(1/(Nh-1))
    factor = (1 - alpha)/(1 - alpha**Nh)
    delta_x = 0.5*factor*(xb - xa)
    delta_y = 0.5*factor*(yb - ya)
    
    X = np.zeros(N + 1)
    Y = np.zeros(N + 1)

    X[0] = xa
    Y[0] = ya
    
    for i in range(1, Nh+1):
        X[i] = X[i-1] + delta_x
        Y[i] = Y[i-1] + delta_y
        delta_x = alpha*delta_x
        delta_y = alpha*delta_y

    delta_x /= alpha
    delta_y /= alpha
        
    for i in range(Nh+1, N + 1):
        X[i] = X[i-1] + delta_x
        Y[i] = Y[i-1] + delta_y
        delta_x = delta_x/alpha
        delta_y = delta_y/alpha

    return X, Y

def potential(point_0, bound_param):
    """
    Computes potential for point outside (not on) obstacle boundary
    using the boundary integral form

    Parameters:
    point_0: point where the potential must be evaluated
    should not sit on a boundary element, otherwise the result is not
    reliable.
    P_long: element nodes.
    dfdn: normal derivative of potential on boundary
    f: potential on boundary

    Returns:
    potential at point_0
    """
    P_long, dfdn, f = bound_param
    Ne = P_long.shape[1] - 1 # number of elements

    sum = 0
    for i in range(Ne):
        alpha, beta = influence_coefficients(P_long[:,i], P_long[:,i+1],
                                             point_0, is_singular=False)
        sum += -alpha*dfdn[i] + beta*f[i]

    return sum

def velocity(point_0, parameters):
    Ux, Uy = parameters[:2]
    bound_param = parameters[2:]
    eps = 1e-2
    
    point_xa = point_0 + eps*np.array([-1,0])
    point_xb = point_0 + eps*np.array([1,0])
    p_xa = potential(point_xa, bound_param) 
    p_xb = potential(point_xb, bound_param) 
    vx = (p_xb - p_xa)/(2*eps)

    point_ya = point_0 + eps*np.array([0,-1])
    point_yb = point_0 + eps*np.array([0,1])
    p_ya = potential(point_ya, bound_param) 
    p_yb = potential(point_yb, bound_param) 
    vy = (p_yb - p_ya)/(2*eps)

    return np.array([vx + Ux, vy + Uy])

def streamline(point_0, x_range, y_range, parameters):
    """
    Compute streamline, integrating differential equation
    by Heun's method
    """
    d = 0.05

    xa, xb = x_range
    ya, yb = y_range
    
    x, y = point_0
    
    if not ((xa <= x <= xb) and (ya <= y <= yb)):
        print("Initial condition out of range")
        sys.exit(1)
    
    point_0 = np.array(point_0)
    point_list = [point_0]

    
    while (xa <= x <= xb) and (ya <= y <= yb):
        vel_0 = velocity(point_0, parameters)
        dt = d/np.linalg.norm(vel_0) # adaptive time step
        point_1 = point_0 + dt*vel_0
        vel_1 = velocity(point_1, parameters)
        point_0 = point_0 + dt/2*(vel_0 + vel_1) # new point
        point_list.append(point_0)
        x, y = point_0

    return np.array(point_list)
    

if __name__ == "__main__":
    # Defining rectangular obstacle boundary
    half_size_x = 0.8
    half_size_y = 0.5
    x_center = 0
    y_center = 0.3
    theta = 0.125*pi # pi/8 radians


    v1 = [x_center + half_size_x, y_center - half_size_y]
    v2 = [x_center + half_size_x, y_center + half_size_y]
    v3 = [x_center - half_size_x, y_center + half_size_y]
    v4 = [x_center - half_size_x, y_center - half_size_y]

    vertices = np.array([v1, v2, v3, v4]).T # vertices in columns

    # defining rotation matrix
    cs = np.cos(theta)
    sn = np.sin(theta)
    rot_mat = np.array([[cs, sn],[-sn, cs]])

    rot_vert = rot_mat.dot(vertices) # vertices in columns

    plt.figure(1)
    plt.clf()
    plt.plot(vertices[0,:], vertices[1,:], 'o-')
    plt.plot(rot_vert[0,:], rot_vert[1,:], 'o-')
    plt.grid(True)

    N = 16 # number of elements of each side, must be even
    ratio = 2
    # each rectangle side will be divided into N elements
    # such that
    # ratio = (length of element N/2)/(length of element 1)
    #       = (length of element N/2+1)/(length of elemnt N)
    
    X = np.zeros((N+1,4))
    Y = np.zeros((N+1,4))
    for i in range(4):
        X[:,i], Y[:,i] = elements_line(rot_vert[:,i], rot_vert[:, (i+1) % 4], N, ratio)

    plt.figure(1)
    for i in range(4):
        plt.plot(X[:,i], Y[:,i], 's')
        plt.grid(True)

    # turn matrix with elements of rectangle sides into a single
    # two-column matrix
    X_long = np.concatenate((X[:,0], X[1:,1], X[1:,2], X[1:,3]))
    Y_long = np.concatenate((Y[:,0], Y[1:,1], Y[1:,2], Y[1:,3]))

    plt.figure(2)
    plt.clf()
    plt.plot(X_long, Y_long, 's-')
    plt.grid(True)

    Ne = X_long.shape[0] - 1

    A = np.zeros((Ne, Ne))
    B = np.zeros((Ne, Ne))

    P_long = np.row_stack((X_long, Y_long))
    
    for i in range(Ne):
        point_0 = 0.5*(P_long[:,i] + P_long[:,i+1])
        for j in range(Ne):
            point_a = P_long[:,j]
            point_b = P_long[:,j+1]
            A[i,j], B[i,j] = influence_coefficients(point_a, point_b,
                                                    point_0, i == j)


    M = B - 0.5*np.identity(Ne)
    
    # velocity of uniform flow
    Ux = 1.0
    Uy = 0.0
    
    dfdn = np.zeros(Ne)
    
    for i in range(Ne):
        xa, ya = P_long[:,i]
        xb, yb = P_long[:,i+1]
        
        dx = xb - xa
        dy = yb - ya
        dr = np.sqrt(dx**2 + dy**2)
        
        # unit vector normal to the element that points outwards
        nx = dy/dr
        ny = -dx/dr
        
        dfdn[i] = - (Ux*nx + Uy*ny)
    
    rhs = A.dot(dfdn)
    f = np.linalg.solve(M, rhs)


    parameters = (Ux, Uy, P_long, dfdn, f)
    x_range = [-2,2]
    y_range = [-2,2]

    # line = streamline([-2,0], x_range, y_range, parameters)
    
    plt.figure(2)
    # plt.clf()
    # plt.plot(line[:,0], line[:,1])
    # plt.grid(True)

    y_list = [-1.3, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2,
              0.4, 0.6, 0.8, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    icond_list = [(-2,y) for y in y_list]

    lines_list = []
    for icond in icond_list:
        print(icond)
        line = streamline(icond, x_range, y_range, parameters)
        lines_list.append(line)
        plt.plot(line[:,0], line[:,1])
