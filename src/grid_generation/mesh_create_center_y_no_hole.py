"""
Mesh Generation with Pygmsh

This script creates a 2D mesh for a rectangular domain, center Y
Pygmsh and OpenCASCADE geometry kernel. It demonstrates how to:
- Create complex geometries with boolean operations (rectangle minus circle)
- Define physical groups for boundary conditions and domain identification
- Generate triangular meshes with controlled element size
- Visualize the generated mesh using matplotlib
- Export meshes in GMSH format for use with finite element solvers

The geometry consists of a rectangle centered at origin.

@author: kssgarcia
"""
# %%
import pygmsh
import matplotlib.pyplot as plt

folder = "../../meshes/"

# Parámetros geométricos
L = 100.0    # Ancho total del rectángulo
b = 5.0     # Mitad de alto → eje y de -b a b
a = 0.0     # Altura del centro del agujero
xc, yc = 10.0, 0.0 # Center of circle
r = 2.0     # Radio del agujero
lc = 0.05   # Tamaño de malla

with pygmsh.occ.Geometry() as geom:
    geom.characteristic_length_min = lc
    geom.characteristic_length_max = lc

    # Rectángulo: x de 0 a L, y de -b a b
    rectangle = geom.add_rectangle([0.0, -b, 0.0], L, 2*b)

    # Región final = rectángulo menos agujero
    domain = rectangle

    # Etiquetas físicas
    geom.add_physical(domain, label="domain")
    geom.add_physical(rectangle, label="outer")

    # Generar y guardar malla
    mesh = geom.generate_mesh()
    pygmsh.write(f"{folder}mesh_no_hole.msh")

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
