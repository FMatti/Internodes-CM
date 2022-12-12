import numpy as np
import matplotlib.pyplot as plt

def get_theoretical_contact_radius(R, d):
    return (d*R)**0.5

def get_theoretical_pressure_amplitude(R, d, E, nu):
    a = get_theoretical_contact_radius(R, d)
    E_red = E/(2*(1-nu**2))
    F = 4/3 * E_red * R**0.5 * d**1.5
    return 3*F / (2*np.pi*a**2)

def get_theoretical_normal_displacement(R, d, E, nu):
    a = get_theoretical_contact_radius(R, d)
    p0 = get_theoretical_pressure_amplitude(R, d, E, nu)
    K = np.pi/2
    return - (1-nu**2)/E * K * p0 * a

def plot_mesh(positions, triangle_indices, nodes_interface=None, nodes_interpenetrating=None, nodes_tension=None):
    plt.figure()
    plt.triplot(positions[:, 0], positions[:, 1], triangles=triangle_indices)
    if nodes_interface is not None:
        plt.scatter(positions[nodes_interface, 0], positions[nodes_interface, 1], color="blue")
    if nodes_interpenetrating is not None and nodes_tension is not None:
        plt.scatter(positions[nodes_interpenetrating, 0], positions[nodes_interpenetrating, 1], color="red")
        plt.scatter(positions[nodes_tension, 0], positions[nodes_tension, 1], color="green")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('scaled')
    plt.show()

def write_solution(mesh, new_nodal_positions):
    with open(mesh, 'r+') as file:
        meshdata = file.read().split('\n')

    nodes_start = meshdata.index("$Nodes")
    nodes_end = meshdata.index("$EndNodes")
    nodal_positions = meshdata[nodes_start+1:nodes_end]
    n_nodes = int(meshdata[nodes_start+1].split(' ')[1])

    node_file_indices = np.empty(n_nodes, dtype=int)
    for i in range(n_nodes):
        node_file_indices[i] = nodal_positions.index(str(i+1))
    node_file_indices = np.append(node_file_indices, len(nodal_positions)+1)
    node_file_indices_spacing = np.diff(node_file_indices)

    if new_nodal_positions.shape[1] < 3:
        new_nodal_positions = np.c_[new_nodal_positions, np.zeros(n_nodes)]

    new_nodal_positions_formatted = ["{} {} {}".format(p[0], p[1], p[2]) for p in new_nodal_positions]

    i = 0
    k = 0
    for d in node_file_indices_spacing:
        if d > 1:
            nodal_positions[node_file_indices[i]+1:node_file_indices[i]+d-1] = new_nodal_positions_formatted[k:k+d-2]
            k += d-2
        i += 1

    with open(mesh.replace(".msh", "_solved.msh"), 'w') as file:
        # Calling the read() method for file objects, to read all its content.
        meshdata[nodes_start+1:nodes_end] = nodal_positions
        file.write("\n".join(meshdata))