import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def read_node_file(filename):
    """
    Reads a .node file from Triangle.
    The first non-comment line is the header:
      <# of nodes> <dimension> <# of attributes> <# of boundary markers>
    Each subsequent line gives:
      <node id> <x> <y> [attributes] [boundary marker]
    Returns a dictionary mapping node IDs to (x,y) tuples.
    """
    nodes = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Filter out blank lines and comments.
    lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
    if not lines:
        raise ValueError("No data in node file.")

    header = lines[0].split()
    num_nodes = int(header[0])
    # Process each node line.
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 3:
            continue
        node_id = int(parts[0])
        x = float(parts[1])
        y = float(parts[2])
        nodes[node_id] = (x, y)
    if len(nodes) != num_nodes:
        print(f"Warning: Expected {num_nodes} nodes, but got {len(nodes)}")
    return nodes

def read_ele_file(filename):
    """
    Reads a .ele file from Triangle.
    The first non-comment line is the header:
      <# of elements> <nodes per element> <# of attributes>
    Each subsequent line gives:
      <element id> <node id 1> <node id 2> ... <node id n> [attributes]
    Returns a list of elements, where each element is a list of node IDs.
    """
    elements = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Filter out blank lines and comments.
    lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
    if not lines:
        raise ValueError("No data in ele file.")

    header = lines[0].split()
    num_elements = int(header[0])
    nodes_per_elem = int(header[1])

    for line in lines[1:]:
        parts = line.split()
        if len(parts) < nodes_per_elem + 1:
            continue
        # Skip the element ID and get the node IDs.
        node_ids = list(map(int, parts[1:1+nodes_per_elem]))
        elements.append(node_ids)
    if len(elements) != num_elements:
        print(f"Warning: Expected {num_elements} elements, but got {len(elements)}")
    return elements

def plot_mesh(nodes, elements):
    """
    Plots the triangular mesh.
    nodes: dict mapping node ID -> (x,y)
    elements: list of elements (each a list of node IDs)
    """
    # Build a list of triangles (each triangle is a list of (x,y) pairs)
    triangles = []
    for elem in elements:
        try:
            triangle = [nodes[node_id] for node_id in elem]
            triangles.append(triangle)
        except KeyError as e:
            print(f"Missing node: {e}")

    _, ax = plt.subplots()
    # Create a PolyCollection from the list of triangles
    poly = PolyCollection(triangles, edgecolors='black', facecolors='none', linewidths=0.8)
    ax.add_collection(poly)

    # Set the plot limits based on node coordinates.
    xs = [coord[0] for coord in nodes.values()]
    ys = [coord[1] for coord in nodes.values()]
    margin = 0.05 * (max(xs) - min(xs))
    ax.set_xlim(min(xs)-margin, max(xs)+margin)
    ax.set_ylim(min(ys)-margin, max(ys)+margin)
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_title("Triangle Mesh")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.show()

def main():
    # Change these filenames to match your output files from Triangle.
    node_filename = "C:/Users/dylan/OneDrive/Documents/uni_work/FEMwithCPP/Code/FEMwithClasses/OOPFEM1D/main/domain.1.node"
    ele_filename = "C:/Users/dylan/OneDrive/Documents/uni_work/FEMwithCPP/Code/FEMwithClasses/OOPFEM1D/main/domain.1.ele"

    nodes = read_node_file(node_filename)
    elements = read_ele_file(ele_filename)
    plot_mesh(nodes, elements)

if __name__ == '__main__':
    main()
