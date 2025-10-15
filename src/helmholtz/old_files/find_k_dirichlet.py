# %%
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import numpy as np
import matplotlib.pyplot as plt
import src.grid_generation.gen_mesh_center_y as gen
import src.helmholtz.full_dirichlet as dtn

bc_type = 'neuman_dirichlet'
L=5.0
b=1.0
xc=0.0
yc=0.0
r=0.1
lc=0.2

det_list = []
kb_list = np.arange(1, 4, 0.05)

nodes, elements = gen.mesh_with_obstacle_center(L=L, b=b, xc=xc, yc=yc, r=r, lc=lc)
for k in kb_list:
    u_full, interior_nodes, nodes, elements, boundary_nodes, A_bc, M_bc = \
        dtn.solve_helmholtz_eigenproblem(nodes, elements, bc_type=bc_type, b=b, k_guess=k)

    A_dense = A_bc.toarray()
    M_dense = M_bc.toarray()

    # Valor propio más pequeño de A - k^2 M
    eigvals = np.linalg.eigvals(A_dense - (k**2) * M_dense)
    eigvals_min = np.min(np.abs(eigvals))

    # LHS = (A_bc - k**2 * M_bc)
    # lu = splu(LHS)
    # logdet = np.sum(np.log(np.abs(lu.U.diagonal())))
    # sign = np.prod(np.sign(lu.U.diagonal()))
    # det = sign * np.exp(logdet)
    tolerance = 1e-6

    det_list.append(eigvals_min)
    if np.abs(eigvals_min) < tolerance:
        print("Matrix is near-singular (det ≈ 0)")
    else:
        print("Matrix is not singular")

y_mag = np.abs(det_list)

# Plot
plt.figure(figsize=(8,5))
plt.plot(kb_list, y_mag, marker='o')
plt.xlabel('k')
plt.ylabel('|det(LHS)|')
plt.title('Determinant magnitude vs k ')
plt.grid(True)
plt.show()

idx_eig = np.argmin(y_mag)
print(np.min(y_mag))
print(kb_list[idx_eig])

# %%
bc_type = 'mixed'
L = 5.0
b = 1.0
a = 0.0
xc = 0.0
yc = 0.0
r = 0.1
lc = 0.1

# Valores de k para barrido
kb_list = np.arange(1, 4, 0.05)  # Puedes ajustar el paso

# Generar malla
nodes, elements = gen.mesh_with_obstacle_center(L=L, b=b, xc=xc, yc=yc, r=r, lc=lc)

# Almacenar resultados
modes_trapped = []
k_trapped = []

# Barrido de k
for k_guess in kb_list:
    try:
        # Resolver Helmholtz
        u_full, interior_nodes, nodes_array, elements_array, boundary_nodes, A_bc, M_bc = \
            dtn.solve_helmholtz_eigenproblem(nodes, elements,
                                         bc_type=bc_type,
                                         b=b,
                                         k_guess=k_guess)

        # Criterio para modo atrapado:
        # Si la solución en los nodos de frontera DtN (left + right) es pequeña, puede ser atrapado
        top_nodes, bottom_nodes, left_nodes, right_nodes, obstacle_nodes = boundary_nodes
        edge_nodes = np.concatenate([left_nodes, right_nodes])
        edge_norm = np.linalg.norm(u_full[edge_nodes])

        # Umbral arbitrario: si es pequeño -> modo atrapado
        if edge_norm < 1e-2:
            modes_trapped.append(u_full)
            k_trapped.append(k_guess)
            print(f"Modo atrapado detectado en k ≈ {k_guess:.3f}, norma en bordes = {edge_norm:.3e}")

    except Exception as e:
        print(f"Error en k={k_guess:.3f}: {e}")
        continue

print(f"\nSe encontraron {len(modes_trapped)} modos atrapados en el barrido de k")
