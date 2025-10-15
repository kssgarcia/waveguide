"""
Mesh Generation with Pygmsh

This script creates a 2D mesh for a rectangular domain, center in XY, with a circular hole using
Pygmsh and OpenCASCADE geometry kernel. It demonstrates how to:
- Create complex geometries with boolean operations (rectangle minus circle)
- Define physical groups for boundary conditions and domain identification
- Generate triangular meshes with controlled element size
- Visualize the generated mesh using matplotlib
- Export meshes in GMSH format for use with finite element solvers

The geometry consists of a rectangle centered at origin with a circular hole.

@author: kssgarcia
"""
# %%
import pygmsh
import matplotlib.pyplot as plt

folder = "../../meshes/"

# Parámetros geométricos
L = 50.0    # Mitad de ancho del rectángulo
b = 5.0    # Mitad de alto del rectángulo
a = 0.0    # Altura del centro del agujero
r = 2.0    # Radio del agujero
lc = 0.05  # Tamaño de malla

with pygmsh.occ.Geometry() as geom:
    geom.characteristic_length_min = lc
    geom.characteristic_length_max = lc

    # Rectángulo centrado en (0, 0)
    rectangle = geom.add_rectangle([-L, -b, 0.0], 2 * L, 2 * b)

    # Agujero: disco centrado en (0, a)
    hole = geom.add_disk([0.0, a, 0.0], r)

    # Región final = rectángulo menos agujero
    domain = geom.boolean_difference([rectangle], [hole])
    # Etiquetas físicas (para condiciones de borde o dominio)
    geom.add_physical(domain, label="domain")
    geom.add_physical(hole, label="hole")
    geom.add_physical(rectangle, label="outer")

    # Generar y guardar malla
    mesh = geom.generate_mesh()
    pygmsh.write(f"{folder}out_r_bigger_center.msh")
# meshio.write("out.vtk", mesh)

points = mesh.points[:, :2]  # drop z-coordinates
cells = mesh.get_cells_type("triangle")  # extract triangle elements
print(cells)

plt.figure(figsize=(6, 6))
plt.triplot(points[:, 0], points[:, 1], cells, color="black", linewidth=0.5)
plt.gca().set_aspect("equal")
plt.title("Triangular Mesh")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# @author: kssgarcia
