import numpy as np

def get_theoretical_contact_radius(R, d):
    return (d * R) ** 0.5

def get_theoretical_pressure_amplitude(R, d, E, nu):
    a = get_theoretical_contact_radius(R, d)
    E_red = E / (2 * (1 - nu**2))
    F = 4 / 3 * E_red * R**0.5 * d**1.5
    return 3 * F / (2 * np.pi * a**2)

def get_theoretical_normal_displacement(R, d, E, nu):
    a = get_theoretical_contact_radius(R, d)
    p0 = get_theoretical_pressure_amplitude(R, d, E, nu)
    K = np.pi / 2
    return (1 - nu**2) / E * K * p0 * a

def plot_mesh(
    ax,
    nodal_positions,
    connectivity,
    nodes_primary=None,
    nodes_secondary=None,
    nodes_candidate_primary=None,
    nodes_candidate_secondary=None,
    nodes_interface_primary=None,
    nodes_interface_secondary=None,
    nodes_added_primary=None,
    nodes_added_secondary=None,
    nodes_dumped_primary=None,
    nodes_dumped_secondary=None
):
    """Plot a mesh."""

    if nodes_primary is None and nodes_secondary is None:
        ax.triplot(nodal_positions[:, 0], nodal_positions[:, 1], triangles=connectivity)
    else:
        connectivity_primary = connectivity[
            np.isin(connectivity, nodes_primary).any(axis=1)
        ]
        connectivity_secondary = connectivity[
            np.isin(connectivity, nodes_secondary).any(axis=1)
        ]
        ax.triplot(
            nodal_positions[:, 0],
            nodal_positions[:, 1],
            triangles=connectivity_primary,
            color="#0e437c",
            linewidth=0.5,
        )
        ax.triplot(
            nodal_positions[:, 0],
            nodal_positions[:, 1],
            triangles=connectivity_secondary,
            color="#d55209",
            linewidth=0.5,
        )
    if nodes_candidate_primary is not None:
        ax.scatter(
            nodal_positions[nodes_candidate_primary, 0],
            nodal_positions[nodes_candidate_primary, 1],
            s=20,
            facecolors="none",
            edgecolors="#0e437c",
            linewidth=0.5,
        )
    if nodes_candidate_secondary is not None:
        ax.scatter(
            nodal_positions[nodes_candidate_secondary, 0],
            nodal_positions[nodes_candidate_secondary, 1],
            s=20,
            facecolors="none",
            edgecolors="#d55209",
            linewidth=0.5,
        )
    if nodes_interface_primary is not None:
        ax.scatter(
            nodal_positions[nodes_interface_primary, 0],
            nodal_positions[nodes_interface_primary, 1],
            s=20,
            color="#0e437c",
        )
    if nodes_interface_secondary is not None:
        ax.scatter(
            nodal_positions[nodes_interface_secondary, 0],
            nodal_positions[nodes_interface_secondary, 1],
            s=20,
            color="#d55209",
        )
    if nodes_added_primary is not None:
        ax.scatter(
            nodal_positions[nodes_added_primary, 0],
            nodal_positions[nodes_added_primary, 1],
            s=20,
            facecolors="none",
            edgecolors="red",
        )
    if nodes_added_secondary is not None:
        ax.scatter(
            nodal_positions[nodes_added_secondary, 0],
            nodal_positions[nodes_added_secondary, 1],
            s=20,
            facecolors="none",
            edgecolors="red",
        )
    if nodes_dumped_primary is not None:
        ax.scatter(
            nodal_positions[nodes_dumped_primary, 0],
            nodal_positions[nodes_dumped_primary, 1],
            s=20,
            marker="x",
            color="green",
        )
    if nodes_dumped_secondary is not None:
        ax.scatter(
            nodal_positions[nodes_dumped_secondary, 0],
            nodal_positions[nodes_dumped_secondary, 1],
            s=20,
            marker="x",
            color="green",
        )

def write_mesh(mesh, new_nodal_positions):
    """Update nodal positions from a .msh mesh file."""
    with open(mesh, "r+") as file:
        meshdata = file.read().split("\n")

    nodes_start = meshdata.index("$Nodes")
    nodes_end = meshdata.index("$EndNodes")
    nodal_positions = meshdata[nodes_start + 1 : nodes_end]
    n_nodes = int(meshdata[nodes_start + 1].split(" ")[1])

    if new_nodal_positions.shape[1] < 3:
        new_nodal_positions = np.c_[new_nodal_positions, np.zeros(n_nodes)]

    new_nodal_positions_formatted = [
        "{} {} {}".format(p[0], p[1], p[2]) for p in new_nodal_positions
    ]

    i = 1
    k = 0
    while i < len(nodal_positions):
        n_nodes_in_group = int(nodal_positions[i].split(" ")[-1])
        if n_nodes_in_group == 0:
            i += 1
            continue
        nodal_positions[
            i + 1 + n_nodes_in_group : i + 1 + 2 * n_nodes_in_group
        ] = new_nodal_positions_formatted[k : k + n_nodes_in_group]
        i += 1 + 2 * n_nodes_in_group
        k += n_nodes_in_group

    with open(mesh.replace(".msh", "_solved.msh"), "w") as file:
        # Calling the read() method for file objects, to read all its content.
        meshdata[nodes_start + 1 : nodes_end] = nodal_positions
        file.write("\n".join(meshdata))
