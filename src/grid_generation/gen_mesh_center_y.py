"""
Mesh Generation with Pygmsh

@author: kssgarcia
"""
# %%
import pygmsh
from meshio import CellBlock
import meshio
import matplotlib.pyplot as plt
import numpy as np

def mesh_with_obstacle(L=100.0, b=5.0, xc=10.0, yc=0.0, r=2.0, lc=0.05, plot=False):
    #=========== Creating Mesh ===========
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = lc
        geom.characteristic_length_max = lc

        # Rectángulo: x de 0 a L, y de -b a b
        rectangle = geom.add_rectangle([0.0, -b, 0.0], L, 2*b)

        # Agujero: disco centrado en x medio del rectángulo
        hole = geom.add_disk([xc, yc, 0.0], r)

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

def mesh_with_obstacle_center(L=5.0, b=1.0, xc=0.0, yc=0.0, r=0.1, lc=0.05, plot=False):
    #=========== Creating Mesh ===========
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = lc
        geom.characteristic_length_max = lc

        # Rectángulo: x de 0 a L, y de -b a b
        rectangle = geom.add_rectangle([-L, -b, 0.0], 2 * L, 2 * b)

        # Agujero: disco centrado en x medio del rectángulo
        hole = geom.add_disk([xc, yc, 0.0], r)

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

def mesh_with_obstacle_center_translated(L=5.0, b=1.0, xc=0.0, yc=0.0, r=0.1, lc=0.05, plot=False):
    #=========== Creating Mesh ===========
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = lc
        geom.characteristic_length_max = lc

        # Rectángulo: x de 0 a L, y de -b a b
        rectangle = geom.add_rectangle([-L, 0.0, 0.0], 2 * L, 2 * b)

        # Agujero: disco centrado en x medio del rectángulo
        hole = geom.add_disk([xc, 1+yc, 0.0], r)

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

def mesh_with_two_obstacle(L=5.0, b=1.0, xc=0.0, yc=0.0, r=0.1, lc=0.05, plot=False):
    #=========== Creating Mesh ===========
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = lc
        geom.characteristic_length_max = lc

        # Rectángulo: x de 0 a L, y de -b a b
        rectangle = geom.add_rectangle([-L, -b, 0.0], 2 * L, 2 * b)

        # Agujero: disco centrado en x medio del rectángulo
        hole1 = geom.add_disk([xc, yc, 0.0], r)
        hole2 = geom.add_disk([-xc, yc, 0.0], r)

        # Región final = rectángulo menos agujero
        domain = geom.boolean_difference([rectangle], [hole1, hole2])

        # Etiquetas físicas
        geom.add_physical(domain, label="domain")
        geom.add_physical(hole1, label="hole1")
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

def mesh_without_obstacle(L=100.0, b=5.0, xc=10.0, yc=0.0, r=2.0, lc=0.05, plot=False):
    #=========== Creating Mesh ===========
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


def mesh_with_parametric_obstacle_even_y(L=5.0, b=1.0, beta=0.2, xc=0.0, yc=0.0, lc=0.05, n_points=200, scale=0.5, plot=False):
    """
    Create a triangular mesh for a rectangle with a parametric obstacle.
    
    scale : float
        Scaling factor for the obstacle size.
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = lc
        geom.characteristic_length_max = lc

        # Rectangle
        rectangle = geom.add_rectangle([-L, -b, 0.0], 2*L, 2*b)

        # Parametric obstacle
        t = np.linspace(0, 2*np.pi, n_points, endpoint=False)
        X = scale * (np.sin(t) - (beta/2)*np.sin(2*t)) + xc
        Y = scale * (-np.cos(t) + (beta/2)*np.cos(2*t)) + yc
        obstacle_points = np.column_stack([X, Y, np.zeros_like(X)])

        # Obstacle loop
        obstacle_loop = geom.add_polygon(obstacle_points.tolist(), make_surface=True)

        # Domain = rectangle - obstacle
        domain = geom.boolean_difference([rectangle], [obstacle_loop])

        # Physical tags
        geom.add_physical(domain, label="domain")
        geom.add_physical(obstacle_loop, label="obstacle")
        geom.add_physical(rectangle, label="outer")

        mesh = geom.generate_mesh()
        nodes = mesh.points
        elements = mesh.get_cells_type("triangle")
        node_result = nodes[:, :2]

    if plot:
        plt.figure(figsize=(8,4))
        plt.triplot(node_result[:,0], node_result[:,1], elements, color="black", linewidth=0.5)
        plt.gca().set_aspect("equal")
        plt.title("Triangular Mesh with Scaled Parametric Obstacle")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()

    return node_result, elements

def mesh_with_parametric_obstacle_even_x(L=5.0, b=1.0, beta=0.2, xc=0.0, yc=0.0, lc=0.05, n_points=200, scale=0.5, plot=False):
    """
    Create a triangular mesh for a rectangle with a parametric obstacle.
    
    scale : float
        Scaling factor for the obstacle size.
    """
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_min = lc
        geom.characteristic_length_max = lc
        # Rectangle
        rectangle = geom.add_rectangle([-L, -b, 0.0], 2*L, 2*b)
        # Parametric obstacle (even in X, odd in Y)
        t = np.linspace(0, 2*np.pi, n_points, endpoint=False)
        X = scale * (-np.cos(t) + (beta/2)*np.cos(2*t)) + xc
        Y = scale * (np.sin(t) - (beta/2)*np.sin(2*t)) + yc
        obstacle_points = np.column_stack([X, Y, np.zeros_like(X)])
        # Obstacle loop
        obstacle_loop = geom.add_polygon(obstacle_points.tolist(), make_surface=True)
        # Domain = rectangle - obstacle
        domain = geom.boolean_difference([rectangle], [obstacle_loop])
        # Physical tags
        geom.add_physical(domain, label="domain")
        geom.add_physical(obstacle_loop, label="obstacle")
        geom.add_physical(rectangle, label="outer")
        mesh = geom.generate_mesh()
        nodes = mesh.points
        elements = mesh.get_cells_type("triangle")
        node_result = nodes[:, :2]
    if plot:
        plt.figure(figsize=(8,4))
        plt.triplot(node_result[:,0], node_result[:,1], elements, color="black", linewidth=0.5)
        plt.gca().set_aspect("equal")
        plt.title("Triangular Mesh with Scaled Parametric Obstacle")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()
    return node_result, elements