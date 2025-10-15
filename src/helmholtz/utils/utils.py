import meshio

def load_mesh_meshio(fileName):
    folder = "../../meshes"
    # Read the original mesh
    mesh = meshio.read(f"{folder}/{fileName}.msh")

    print("="*60)
    print("Loading Mesh")
    print("="*60)
    print(mesh)
    print("="*60)

    # Extract node coordinates
    nodes = mesh.points  # Shape: (num_nodes, 3)

    # Extract element connectivity
    elements = mesh.get_cells_type("triangle")  # Shape: (num_elements, 3)

    # Extract physical group information
    print("Physical Groups:")
    if mesh.field_data:
        for name, (tag, dim) in mesh.field_data.items():
            print(f"  {name}: tag={tag}, dimension={dim}")

    # Get physical group tags for elements
    physical_tags = None
    if mesh.cell_data and "gmsh:physical" in mesh.cell_data:
        # Get physical tags for triangular elements
        cell_types = [cell.type for cell in mesh.cells]
        triangle_idx = cell_types.index("triangle") if "triangle" in cell_types else None

        if triangle_idx is not None:
            physical_tags = mesh.cell_data["gmsh:physical"][triangle_idx]
            print(f"\nPhysical tags for elements: {len(physical_tags)} elements")
            print(f"Unique physical tags: {set(physical_tags)}")

    # Alternative way to get cell data
    print(f"\nAvailable cell data keys: {list(mesh.cell_data.keys()) if mesh.cell_data else 'None'}")
    node_result = nodes[:, :2]

    return node_result, elements
