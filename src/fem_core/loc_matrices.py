import numpy as np

def locA(x1, y1, x2, y2, x3, y3):
    """
    Evaluate the element diffusion matrix for P1 element.
    
    Parameters:
    x1, y1: Coordinates of the first node
    x2, y2: Coordinates of the second node
    x3, y3: Coordinates of the third node
    
    Returns:
    np.ndarray: 3x3 element diffusion matrix
    """
    # Calculate differences between coordinates
    dx23 = x2 - x3
    dy23 = y2 - y3
    dx31 = x3 - x1
    dy31 = y3 - y1
    dx12 = x1 - x2
    dy12 = y1 - y2
    
    # Calculate triangle area
    A = 0.5 * (dx31 * dy12 - dy31 * dx12)
    
    # Initialize the element matrix
    elmA = np.zeros((3, 3))
    
    # Fill the element matrix
    elmA[0, 0] = 0.25 * (dx23 * dx23 + dy23 * dy23) / A
    elmA[0, 1] = 0.25 * (dx23 * dx31 + dy23 * dy31) / A
    elmA[0, 2] = 0.25 * (dx23 * dx12 + dy23 * dy12) / A
    
    elmA[1, 0] = 0.25 * (dx31 * dx23 + dy31 * dy23) / A
    elmA[1, 1] = 0.25 * (dx31 * dx31 + dy31 * dy31) / A
    elmA[1, 2] = 0.25 * (dx31 * dx12 + dy31 * dy12) / A
    
    elmA[2, 0] = 0.25 * (dx12 * dx23 + dy12 * dy23) / A
    elmA[2, 1] = 0.25 * (dx12 * dx31 + dy12 * dy31) / A
    elmA[2, 2] = 0.25 * (dx12 * dx12 + dy12 * dy12) / A
    
    return elmA

import numpy as np

def locB(vx, vy, x1, y1, x2, y2, x3, y3):
    """
    Calculate the element matrix from convection term for P1 element.
    
    Parameters:
    vx, vy: Convection velocity components
    x1, y1: Coordinates of the first node
    x2, y2: Coordinates of the second node
    x3, y3: Coordinates of the third node
    
    Returns:
    np.ndarray: 3x3 element convection matrix
    """
    # Calculate differences between coordinates
    dx32 = x3 - x2
    dy32 = y3 - y2
    dx13 = x1 - x3
    dy13 = y1 - y3
    dx21 = x2 - x1
    dy21 = y2 - y1
    
    # Initialize the element matrix
    elmB = np.zeros((3, 3))
    
    # Calculate first row entries
    elmB[0, 0] = (-vx * dy32 + vy * dx32) / 6.0
    elmB[0, 1] = (-vx * dy13 + vy * dx13) / 6.0
    elmB[0, 2] = (-vx * dy21 + vy * dx21) / 6.0
    
    # Fill remaining entries (symmetric)
    elmB[1, 0] = elmB[0, 0]  # Second row, first column
    elmB[2, 0] = elmB[0, 0]  # Third row, first column
    elmB[1, 1] = elmB[0, 1]  # Second row, second column
    elmB[2, 1] = elmB[0, 1]  # Third row, second column
    elmB[1, 2] = elmB[0, 2]  # Second row, third column
    elmB[2, 2] = elmB[0, 2]  # Third row, third column
    
    return elmB

import numpy as np

def locM(x1, y1, x2, y2, x3, y3):
    """
    Calculate the element mass matrix for P1 element.
    
    Parameters:
    x1, y1: Coordinates of the first node
    x2, y2: Coordinates of the second node
    x3, y3: Coordinates of the third node
    
    Returns:
    np.ndarray: 3x3 element mass matrix
    """
    # Calculate differences between coordinates
    dx23 = x2 - x3
    dy23 = y2 - y3
    dx31 = x3 - x1
    dy31 = y3 - y1
    dx12 = x1 - x2
    dy12 = y1 - y2
    
    # Calculate triangle area
    A = 0.5 * (dx31 * dy12 - dy31 * dx12)
    
    # Calculate constants for matrix elements
    c_diag = A / 6      # Diagonal constant
    c_off = A / 12      # Off-diagonal constant
    
    # Initialize and populate the element matrix
    elmM = np.zeros((3, 3))
    for j in range(3):
        for i in range(3):
            if i == j:
                elmM[i, j] = c_diag
            else:
                elmM[i, j] = c_off
                
    return elmM
