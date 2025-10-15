# %%
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.sparse import lil_matrix
from matplotlib.tri import Triangulation
from scipy.sparse import csr_matrix, bmat
from scipy.sparse.linalg import spsolve, eigs, eigsh
import src.helmholtz.utils.utils as utils
import src.grid_generation.gen_mesh_center_y as gen
from scipy.sparse.linalg import svds

import pygmsh
from meshio import CellBlock
import meshio
# import matplotlib.pyplot as plt
# import numpy as np

def mesh_with_obstacle_center(L=5.0, b=1.0, a=0.0, xc=0.0, yc=0.0, r=0.1, lc=0.05, plot=False):
    #=========== Creating Mesh ===========
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = lc
        geom.characteristic_length_max = lc

        # Rectángulo: x de 0 a L, y de -b a b
        rectangle = geom.add_rectangle([-L, -b, 0.0], 2 * L, 2 * b)

        # Agujero: rectángulo centrado en x medio del rectángulo
        hole = geom.add_rectangle([-r, -r, 0.0], 2*r, 2*r)
        # Agujero: disco centrado en x medio del rectángulo
        # hole = geom.add_disk([xc, yc, 0.0], r)

        # Región final = rectángulo menos agujero
        domain = geom.boolean_difference([rectangle], [hole])

        # Etiquetas físicas
        geom.add_physical(domain, label="domain")
        geom.add_physical(hole, label="hole")
        geom.add_physical(rectangle, label="outer")

        # Generar y guardar malla
        mesh = geom.generate_mesh()

        #=========== Reading Mesh ===========
        points = mesh.points

        # Get triangle data
        triangle_data = mesh.get_cells_type("triangle")

        # Create CellBlock properly - this is the correct format
        triangle_cells = [CellBlock("triangle", triangle_data)]

        # Create the new mesh with proper structure
        mesh_tri_only = meshio.Mesh(
            points=points,  # Keep 3D points
            cells=triangle_cells,
            # Copy over any relevant data from original mesh
            point_data=mesh.point_data if hasattr(mesh, 'point_data') else {},
            cell_data={} if not hasattr(mesh, 'cell_data') else {
                key: [val for i, val in enumerate(values) if mesh.cells[i].type == "triangle"]
                for key, values in mesh.cell_data.items()
            }
        )

        #=========== loading data ===========
        nodes = mesh.points  # Shape: (num_nodes, 3)
        elements = mesh.get_cells_type("triangle")  # Shape: (num_elements, 3)
        node_result = nodes[:, :2]

        print("="*60)
        print("Generating Mesh")
        print("="*60)
        print(mesh)
        print("="*60)

    if plot:
        points = mesh.points[:, :2]
        cells = mesh.get_cells_type("triangle")
        plt.figure(figsize=(8, 4))
        plt.triplot(points[:, 0], points[:, 1], cells, color="black", linewidth=0.5)
        plt.gca().set_aspect("equal")
        plt.title("Triangular Mesh")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()


    return node_result, elements


# Enhanced main execution
if __name__ == "__main__":
    # nodes, elements = utils.load_mesh_meshio("mesh_with_hole")
    # nodes, elements = gen.mesh_with_perturbed_hole(lc=0.5, plot=True)
    bc_type = 'mixed'
    L=5.0
    b=1.0
    a=0.0
    xc=0.0
    yc=0.0
    r=0.3
    lc=0.1
    k_guess = 1.7
    N_modes = 20

    det_list = []
    kb_list = np.arange(1, 4, 0.05)

    nodes, elements = mesh_with_obstacle_center(L=L, b=b, a=a, xc=xc, yc=yc, r=r, lc=lc, plot=True)

    np.savez("square_mesh", nodes=nodes, elements=elements)
