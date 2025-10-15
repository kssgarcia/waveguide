"""
Mesh Generation with Pygmsh

This script creates a 2D mesh for a rectangular domain, center Y, with a perturbe circular hole using
Pygmsh and OpenCASCADE geometry kernel. It demonstrates how to:
- Create complex geometries with boolean operations (rectangle minus circle)
- Define physical groups for boundary conditions and domain identification
- Generate triangular meshes with controlled element size
- Visualize the generated mesh using matplotlib
- Export meshes in GMSH format for use with finite element solvers

The geometry consists of a rectangle centered at origin with a perturbe circular hole.

@author: kssgarcia
"""
# %%
import numpy as np
import pygmsh
import matplotlib.pyplot as plt
import os

# Create folder if it doesn't exist
folder = "../../meshes/"
file = "mesh_with_perturbed_hole_refine"
os.makedirs(folder, exist_ok=True)

# Parameters
L = 100.0    # rectangle width
b = 5.0      # half-height
xc, yc = 10.0, 0.0  # center of the perturbed circle
lc = 0.5     # mesh size (increased for better stability)
beta = 0.2   # perturbation parameter (must be small)
n_points = 800  # reduced number of points for better stability

# Parametric boundary of perturbed circle
t = np.linspace(0, 2*np.pi, n_points, endpoint=False)
X = np.sin(t) - (beta/2) * np.sin(2*t)  
Y = -np.cos(t) + (beta/2) * np.cos(2*t)  

# Scale and shift to desired center (xc, yc) and radius r
r = 2.0
X = xc + r * X
Y = yc + r * Y

def remove_duplicates(X, Y, tol=1e-6):
    """Remove duplicate points that are too close together"""
    pts = [(X[0], Y[0])]
    for x, y in zip(X[1:], Y[1:]):
        if (x - pts[-1][0])**2 + (y - pts[-1][1])**2 > tol**2:
            pts.append((x, y))
    return np.array(pts)[:,0], np.array(pts)[:,1]

X, Y = remove_duplicates(X, Y, tol=1e-6)

# Create mesh using pygmsh
with pygmsh.occ.Geometry() as geom:
    geom.characteristic_length_min = lc
    geom.characteristic_length_max = lc
    
    # Outer rectangle
    rectangle = geom.add_rectangle([0.0, -b, 0.0], L, 2*b)
    
    # Create perturbed hole using polygon approximation
    # Convert to 3D points (z=0) for pygmsh
    hole_points = [[x, y, 0.0] for x, y in zip(X, Y)]
    
    # Create polygon from points
    polygon_points = []
    for point in hole_points:
        polygon_points.append(geom.add_point(point, lc))
    
    # Create lines connecting the points
    lines = []
    for i in range(len(polygon_points)):
        next_i = (i + 1) % len(polygon_points)
        lines.append(geom.add_line(polygon_points[i], polygon_points[next_i]))
    
    # Create curve loop and surface for the hole
    hole_loop = geom.add_curve_loop(lines)
    hole_surface = geom.add_plane_surface(hole_loop)
    
    # Create domain by boolean difference
    domain = geom.boolean_difference([rectangle], [hole_surface])
    
    # Add physical groups
    geom.add_physical(domain, "domain")
    
    # Generate mesh
    mesh = geom.generate_mesh()
    pygmsh.write(f"{folder}{file}.msh")

# Extraer puntos y triángulos
points = mesh.points[:, :2]
cells = mesh.get_cells_type("triangle")

# Mostrar malla
plt.figure(figsize=(8, 4))
plt.triplot(points[:, 0], points[:, 1], cells, color="black", linewidth=0.5)
plt.gca().set_aspect("equal")
plt.title("Triangular Mesh")
plt.xlabel("x")
plt.ylabel("y")
plt.show()
