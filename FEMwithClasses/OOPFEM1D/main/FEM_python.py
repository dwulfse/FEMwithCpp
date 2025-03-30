import numpy as np
import math

# === Mesh reading functions ===

def read_node_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    header = lines[0].split()
    n_nodes = int(header[0])
    nodes = np.zeros((n_nodes, 2))
    for line in lines[1:]:
        if line.strip() == "" or line[0] == '#':
            continue
        parts = line.split()
        idx = int(parts[0]) - 1  # convert to 0-index
        x = float(parts[1])
        y = float(parts[2])
        nodes[idx] = [x, y]
    return nodes

def read_ele_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    header = lines[0].split()
    n_elems = int(header[0])
    elems = []
    for line in lines[1:]:
        if line.strip() == "" or line[0] == '#':
            continue
        parts = line.split()
        # assuming triangles: id, n1, n2, n3
        node_ids = [int(parts[i]) - 1 for i in range(1, 4)]
        elems.append(node_ids)
    return elems

# === Geometry helper functions ===

def compute_triangle_area(P, Q, R):
    return 0.5 * abs(P[0]*(Q[1]-R[1]) + Q[0]*(R[1]-P[1]) + R[0]*(P[1]-Q[1]))

def compute_affine_matrix(P, Q, R):
    A = np.zeros((2,2))
    A[:,0] = Q - P
    A[:,1] = R - P
    return A

def invert_2x2(A):
    det = A[0,0]*A[1,1] - A[0,1]*A[1,0]
    if abs(det) < 1e-12:
        raise ValueError("Singular matrix")
    return np.array([[A[1,1], -A[0,1]],
                     [-A[1,0], A[0,0]]]) / det

# === FEM local functions (using the same reference element as C++ code) ===

# Reference element: vertices (-1,-1), (1,-1), (-1,1); area = 2.
# Shape functions: N1 = -0.5*x - 0.5*y, N2 = 0.5*x + 0.5, N3 = 0.5*y + 0.5.
def shape_functions(qp):
    x, y = qp
    N1 = -0.5 * x - 0.5 * y
    N2 =  0.5 * x + 0.5
    N3 =  0.5 * y + 0.5
    return np.array([N1, N2, N3])

# The reference gradients (constant)
ref_grads = np.array([[-0.5, -0.5],
                      [ 0.5,  0.0],
                      [ 0.0,  0.5]])

# A simple quadrature rule on T_ref (you may replace these with a proper rule)
# Here we provide one point (the centroid) with weight = area.
quad_points = np.array([[ -1/3, -1/3 ]])  # approximate centroid for T_ref (not exact)
quad_weights = np.array([2.0])             # because area(T_ref)=2

def map_to_physical(P, Q, R, qp):
    # For our reference element with vertices (-1,-1), (1,-1), (-1,1),
    # the mapping is:
    # x = P + 0.5*(xi+1)*(Q-P) + 0.5*(eta+1)*(R-P)
    xi, eta = qp
    return P + 0.5*(xi+1)*(Q-P) + 0.5*(eta+1)*(R-P)

def local_stiffness(P, Q, R):
    A = compute_affine_matrix(P, Q, R)
    A_inv = invert_2x2(A)
    # Transform reference gradients: phys_grad = (A_inv)^T * ref_grad
    phys_grads = np.zeros_like(ref_grads)
    for i in range(3):
        phys_grads[i] = A_inv.T @ ref_grads[i]
    area = compute_triangle_area(P, Q, R)
    scale = area / 2.0  # Correction factor
    K_local = scale * (phys_grads @ phys_grads.T)
    return K_local

def local_load(P, Q, R, f):
    area = compute_triangle_area(P, Q, R)
    J = area / 2.0
    F_local = np.zeros(3)
    for qp, w in zip(quad_points, quad_weights):
        phys = map_to_physical(P, Q, R, qp)
        N = shape_functions(qp)
        F_local += w * f(phys[0], phys[1]) * N
    F_local *= J
    return F_local

# Example right-hand side: f(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y)
def f_func(x, y):
    return 2 * math.pi**2 * math.sin(math.pi*x)*math.sin(math.pi*y)

# === Assembly of Global System ===

nodes = read_node_file("C:/Users/dylan/OneDrive/Documents/uni_work/FEMwithCPP/Code/FEMwithClasses/OOPFEM1D/main/domain.1.node")  # adjust file paths as needed
elems = read_ele_file("C:/Users/dylan/OneDrive/Documents/uni_work/FEMwithCPP/Code/FEMwithClasses/OOPFEM1D/main/domain.1.ele")
n_nodes = nodes.shape[0]

K_global = np.zeros((n_nodes, n_nodes))
F_global = np.zeros(n_nodes)

for elem in elems:
    P = nodes[elem[0]]
    Q = nodes[elem[1]]
    R = nodes[elem[2]]
    K_loc = local_stiffness(P, Q, R)
    F_loc = local_load(P, Q, R, f_func)
    for i_local, i_global in enumerate(elem):
        F_global[i_global] += F_loc[i_local]
        for j_local, j_global in enumerate(elem):
            K_global[i_global, j_global] += K_loc[i_local, j_local]

print("Global stiffness matrix:")
print(K_global)
for i in range(n_nodes):
    for j in range(n_nodes):
        if K_global[i][j] != 0:
            print(f"({j}, {K_global[i][j]})", end=" ")
    print()
print("\nGlobal load vector:")
print(F_global)

# === Solve the linear system ===

# Apply Dirichlet boundary conditions (for simplicity, set u=0 on all boundaries)
# Here we simply zero out rows and columns corresponding to boundary nodes.
# In a real application, you would need to modify the system more carefully.
boundary_nodes = [node_id for node_id, (x, y) in enumerate(nodes) if x == 0 or x == 1 or y == 0 or y == 1]
for node_id in boundary_nodes:
    K_global[node_id, :] = 0
    K_global[:, node_id] = 0
    K_global[node_id, node_id] = 1
    F_global[node_id] = 0

# Solve the linear system
u = np.linalg.solve(K_global, F_global)
print("\nSolution (u):")
print(u)
