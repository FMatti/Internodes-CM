"""
References
----------

[1] Y. Voet et. al.: The INTERNODES method for applications in contact mechanics
    and dedicated preconditioning techniques.
    Computers & Mathematics with Applications, vol. 127, 2022, pp. 48-64
    https://doi.org/10.1016/j.camwa.2022.09.019
[2] Y. Voet: On the preconditioning of the INTERNODES matrix for applications
    in contact mechanics.
    Master's thesis EPFL, 2021.
"""

import numpy as np
import scipy as sp
import akantu as aka

def nodes_to_dofs(nodes, dim):
    """Obtain the DOFs that correspond to the node indices

    Parameters
    ----------
    nodes : 1D list or np.ndarray
        List of node indices.
    dim : int
        Spatial dimension.

    Returns
    -------
    DOFs : np.ndarray
        Array with the degrees of freedom corresponding to the nodes.

    Example
    -------
    >>> nodes =  np.array([1, 2, 5])
    >>> nodes_to_dofs(nodes, dim=2)
    >>>   = array([2, 3, 4, 5, 10, 11])
    """
    return dim * np.repeat(nodes, dim) + np.tile(range(dim), len(nodes))

def expand_to_dim(matrix, dim):
    """Expand matrix entries to (dim x dim)-blocks corresponding to the DOFs.

    Parameters
    ----------
    matrix : 2D np.ndarray, sparse or dense scipy matrix
        A matrix.
    dim : int
        The spatial dimension of the problem.

    Returns
    -------
    matrix_expanded : sparse scipy matrix
        The expanded matrix.

    Example
    -------
    >>> matrix = np.array([[1., 2.],
    >>>                    [3., 4.]])
    >>> expand_to_dim(matrix, dim=2).todense()
    >>>   = matrix([[1., 0., 2., 0.],
    >>>             [0., 1., 0., 2.],
    >>>             [3., 0., 4., 0.],
    >>>             [0., 3., 0., 4.]])
    """
    return sp.sparse.kron(matrix, np.eye(dim), format="csr")

def remove_rows_without_all_items(matrix, items):
    """Remove rows from matrix that are not fully contained in items.

    Parameters
    ----------
    matrix : 2D np.ndarray
        A matrix.
    items : 1D np.ndarray or list
        A list with items.

    Returns
    -------
    matrix_removed : 2D np.ndarray
        The matrix where rows not fully contained in item have been eliminated.

    Example
    -------
    >>> matrix = np.array([[1, 2],
    >>>                    [3, 4],
    >>>                    [5, 6],
    >>>                    [7, 8]])
    >>> items = [3, 4, 8]
    >>> remove_rows_without_all_items(matrix, items)
    >>>   = array([[3, 4]])
    """
    return matrix[np.isin(matrix, items).all(axis=1)]

def wendland(delta):
    """Evaluate Wendland C2 function.

    Reference
    ---------
    [1], Page 49, Section 2.1, Table 1
    """
    return (1 - delta) ** 4 * (1 + 4 * delta) * (delta <= 1)

def wendland_rbf(distances, radiuses):
    """Evaluate Wendland radial basis function.

    Reference
    ---------
    [1], Page 49, Section 2.1, Bottom left
    """
    return wendland(distances / radiuses)

def compute_rbf_radius_parameters(positions, rbf=wendland_rbf, c=0.5):
    """Iteratively compute radius parameters for radial basis functions until
    invertibility conditions are satisfied (increase 'c' in every iteration).

    Parameters
    ----------
    positions : np.ndarray
        Positions of interpolation nodes.
        rbf : function, default is wendland
        The radial basis function being used.
    c : float in (0, 1), default is 0.5
        [1], Page 51, Section 2.3, Equation 2
        (The empirical default value is given at the right on same page)

    Returns
    -------
    rbf_radius_parameters : np.ndarray
        The radius parameters to use in the radial basis function.

    Reference
    ---------
    [1], Page 51, Section 2.3
    """
    # Compute distance matrices among nodes
    # distance_matrix(X, Y)[i, j] = ||x_i - y_j||
    distance_matrix_MM = sp.spatial.distance.cdist(positions, positions)

    # Minimum distance of each node from closest distinct interpolation node
    np.fill_diagonal(distance_matrix_MM, np.inf)
    min_distance_MM = np.min(distance_matrix_MM, axis=0)

    # Iteratively increase parameter 'c' if necessary
    while True:

        # Define radius parameters for [1], Page 51, Section 2.3, Equation 2
        # (Also mentioned in [1], Page 51, Section 2.3, right center)
        rbf_radius_parameters = min_distance_MM / c

        # Number radial basis functions in whose support interpolation nodes are
        n_supports_interpolation = np.sum(
            distance_matrix_MM < rbf_radius_parameters, axis=0
        )

        # Check if criterion [1], Page 51, Section 2.3, Equation 4 is satisfied
        # for all interpolation nodes
        if np.max(n_supports_interpolation) < 1 / rbf(c, 1):
            break

        # If criterion was not satisfied, increase c and reiterate the procedure
        c = (c + 1) / 2
        if c >= 1:
            ValueError("Tried to increase c to", c, "which is larger than 1.")

    return rbf_radius_parameters

def find_interpolation_nodes(
    positions_primary, positions_secondary, rbf=wendland_rbf, C=0.95
):
    """Find contact/interface nodes while trying to satisfy the constraints.

    Parameters
    ----------
    positions_primary: np.ndarray
        Positions of candidate primary nodes.
    positions_secondary: np.ndarray
        Positions of candidate secondary nodes.
    rbf : function, default is wendland
        The radial basis function being used.
    C : float in (c, 1), default is 0.5
        [1], Page 51, Section 2.3, Equation 3
        (The empirical default value is given at the bottom right on same page)

    Returns
    -------
    interface_primary_mask : np.ndarray
        Boolean array for nodes in primary that belong to the interface.
    interface_secondary_mask : np.ndarray
        Boolean array for nodes in secondary that belong to the interface.

    Reference
    ---------
    [1], Page 51, Section 2.3, Equation 2 and 3
    """

    interface_primary_mask = np.ones(len(positions_primary), dtype=bool)
    interface_secondary_mask = np.ones(len(positions_secondary), dtype=bool)

    while True:

        # Determine the radial basis function radius parameters and to how
        # many opposite radial basis function supports each node belongs
        rbf_radius_parameters_primary = compute_rbf_radius_parameters(
            positions_primary[interface_primary_mask], rbf=rbf
        )
        rbf_radius_parameters_secondary = compute_rbf_radius_parameters(
            positions_secondary[interface_secondary_mask], rbf=rbf
        )

        distance_matrix_MN = sp.spatial.distance.cdist(
            positions_primary[interface_primary_mask],
            positions_secondary[interface_secondary_mask],
        )

        # Determine isolated nodes (i.e. nodes outside support of all opposite rbf)
        primary_mask = (
            np.min(distance_matrix_MN / rbf_radius_parameters_secondary, axis=1) < C
        )
        secondary_mask = (
            np.min(distance_matrix_MN.T / rbf_radius_parameters_primary, axis=1) < C
        )
        interface_primary_mask[interface_primary_mask] = primary_mask
        interface_secondary_mask[interface_secondary_mask] = secondary_mask

        # Stop algorithm if for all secondary and primary interface nodes
        # [1], Page 51, Section 2.3, Equation 3 is satisfied
        if primary_mask.all() or secondary_mask.all():
            break

        # If all nodes of one interface are isolated, raise an error
        if not interface_primary_mask.any() or not interface_secondary_mask.any():
            raise RuntimeError("No contact nodes were detected.")

    # Remove radius parameters of nodes not belonging to interface
    rbf_radius_parameters_primary = rbf_radius_parameters_primary[primary_mask]
    rbf_radius_parameters_secondary = rbf_radius_parameters_secondary[secondary_mask]

    return (
        rbf_radius_parameters_primary,
        rbf_radius_parameters_secondary,
        interface_primary_mask,
        interface_secondary_mask,
    )

def construct_rbf_matrix(positions, positions_ref, radiuses_ref, rbf=wendland_rbf):
    """Construct radial basis matrix $\Phi_{NM}$.

    Reference
    ---------
    [1], Page 49, Section 2.1, Point 1
    """
    return rbf(
        sp.spatial.distance.cdist(
            positions.reshape(len(positions), -1),
            positions_ref.reshape(len(positions_ref), -1),
        ),
        radiuses_ref,
    )

def construct_interpolation_matrix(
    positions, positions_ref, radiuses_ref, rbf=wendland_rbf
):
    """Construct interpolation matrix $\R_{NM}$.

    Reference
    ---------
    [1], Page 49, Section 2.1, Bottom right
    """
    # Construct radial basis function matrices $\Phi_{MM}$ and $\Phi_{NM}$
    rbf_matrix_MM = construct_rbf_matrix(
        positions_ref, positions_ref, radiuses_ref, rbf
    )
    rbf_matrix_NM = construct_rbf_matrix(positions, positions_ref, radiuses_ref, rbf)

    # Compute raw interpolation matrix without rescaling
    interpolation_matrix = np.linalg.solve(rbf_matrix_MM.T, rbf_matrix_NM.T).T

    # Compute the diagonal rescaling factors (diagonal entries of $D_{NN}$)
    scaling_factors = np.sum(interpolation_matrix, axis=1)[:, np.newaxis]
    return sp.sparse.csr_matrix(interpolation_matrix / scaling_factors)

def compute_normals(nodal_positions, nodes, surface_connectivity, dim):
    """Compute the normals.

    Parameters
    ----------
    nodal_positions : np.ndarray
        Positions of nodes.
    nodes : np.ndarray
        Nodes for which the normals are computed for.
    surface_connectivity : np.ndarray
        Surface connectivity of the boundary elements in the mesh.

    Returns
    -------
    node_normals : np.ndarray
        Average of surface normals of elements adjacent to each node.
    """

    # Compute tangents corresponding to the elements at new positions
    if dim == 2:
        tangent1 = (
            nodal_positions[surface_connectivity[:, 1]]
            - nodal_positions[surface_connectivity[:, 0]]
        )
        tangent2 = [0, 0, 1]
    elif dim == 3:
        tangent1 = (
            nodal_positions[surface_connectivity[:, 1]]
            - nodal_positions[surface_connectivity[:, 0]]
        )
        tangent2 = (
            nodal_positions[surface_connectivity[:, 2]]
            - nodal_positions[surface_connectivity[:, 0]]
        )

    # Compute normal vectors to all the boundary elements
    surface_normals = np.cross(tangent1, tangent2)[:, :dim]

    node_normals = np.zeros((len(nodes), dim))
    for j, node in enumerate(nodes):
        id = np.isin(surface_connectivity, node).any(axis=1)

        # Compute average surface normal for each candidate node
        node_normals[j] = np.sum(
            surface_normals[id]
            / np.linalg.norm(surface_normals[id], axis=1)[:, np.newaxis],
            axis=0,
        )

    node_normals /= np.linalg.norm(node_normals, axis=1)[:, np.newaxis]

    return node_normals

def alternative_find_penetration_nodes(
    positions_primary, positions_secondary, normals_primary, normals_secondary, n=2
):
    """Find the nodes on secondary and primary interface that penetrate.
    TODO: Require reference and make mesh size configurable!!

    Parameters
    ----------
    positions_new : np.ndarray
        New positions of nodes after solving the INTERNODES system.
    normals_interface_primary : np.ndarray
        Normal vectors of primary interface.
    normals_interface_secondary : np.ndarray
        Normal vectors of secondary interface.
    tolerance : float in (0, 1), default is 0.9
        Tolerance for what counts as penetration or not.
    mesh_size : float, default is 0.05 (TODO: No default but adaptive)
        Representative size of mesh to determine when penetration happens.

    Returns
    -------
    nodes_penetration_primary : np.ndarray
        Nodes of primary interface where penetration is observed.
    nodes_penetration_secondary : np.ndarray
        Nodes of secondary interface where penetration is observed.
    """

    distance_matrix = sp.spatial.distance.cdist(positions_primary, positions_secondary)

    closest_secondary_to_primary_nodes = np.argsort(distance_matrix, axis=1)[:, :n]
    closest_primary_to_secondary_nodes = np.argsort(distance_matrix.T, axis=1)[:, :n]

    closest_secondary_to_primary_positions = positions_secondary[
        closest_secondary_to_primary_nodes
    ]
    closest_primary_to_secondary_positions = positions_primary[
        closest_primary_to_secondary_nodes
    ]

    gaps_primary_to_secondary = (
        closest_secondary_to_primary_positions - positions_primary[:, np.newaxis, :]
    )
    gaps_secondary_to_primary = (
        closest_primary_to_secondary_positions - positions_secondary[:, np.newaxis, :]
    )

    inner_primary = np.sum(
        gaps_primary_to_secondary
        * normals_secondary[closest_secondary_to_primary_nodes],
        axis=2,
    )
    inner_secondary = np.sum(
        gaps_secondary_to_primary * normals_primary[closest_primary_to_secondary_nodes],
        axis=2,
    )

    penetration_node_primary = ~np.any(inner_primary < -1e-3, axis=1)
    penetration_node_secondary = ~np.any(inner_secondary < -1e-3, axis=1)

    return penetration_node_primary, penetration_node_secondary

def find_penetration_nodes(
    positions_primary,
    positions_secondary,
    normals_candidates_primary,
    normals_candidates_secondary,
    rbf=wendland_rbf,
    tolerance=0.9,
    mesh_size=0.02,
):
    """Find the nodes on secondary and primary interface that penetrate.
    TODO: Require reference and make mesh size configurable!!

    Parameters
    ----------
    positions_new : np.ndarray
        New positions of nodes after solving the INTERNODES system.
    normals_interface_primary : np.ndarray
        Normal vectors of primary interface.
    normals_interface_secondary : np.ndarray
        Normal vectors of secondary interface.
    tolerance : float in (0, 1), default is 0.9
        Tolerance for what counts as penetration or not.
    mesh_size : float, default is 0.05 (TODO: No default but adaptive)
        Representative size of mesh to determine when penetration happens.

    Returns
    -------
    nodes_penetration_primary : np.ndarray
        Nodes of primary interface where penetration is observed.
    nodes_penetration_secondary : np.ndarray
        Nodes of secondary interface where penetration is observed.
    """

    # Find contact nodes with contact detection algorithm
    (
        rbf_radius_parameters_primary,
        rbf_radius_parameters_secondary,
        nodes_interface_primary_mask,
        nodes_interface_secondary_mask,
    ) = find_interpolation_nodes(positions_primary, positions_secondary)

    # Update positions of primary and secondary nodes along interface
    positions_interface_primary = positions_primary[nodes_interface_primary_mask]
    positions_interface_secondary = positions_secondary[nodes_interface_secondary_mask]

    # Update normals of primary and secondary nodes along interface
    normals_interface_primary = normals_candidates_primary[nodes_interface_primary_mask]
    normals_interface_secondary = normals_candidates_secondary[
        nodes_interface_secondary_mask
    ]

    R12 = construct_interpolation_matrix(
        positions_interface_primary,
        positions_interface_secondary,
        rbf_radius_parameters_secondary,
        rbf,
    )
    R21 = construct_interpolation_matrix(
        positions_interface_secondary,
        positions_interface_primary,
        rbf_radius_parameters_primary,
        rbf,
    )

    # Determine the size of the nodal gaps on interface
    nodal_gaps_primary = (
        R12 * positions_interface_secondary - positions_interface_primary
    )
    nodal_gaps_secondary = (
        R21 * positions_interface_primary - positions_interface_secondary
    )

    # Penetration if projected nodal gap onto normal is sufficiently small
    threshold = -tolerance * mesh_size
    penetrates_secondary = (
        np.sum(nodal_gaps_primary * normals_interface_primary, axis=1) < threshold
    )
    penetrates_primary = (
        np.sum(nodal_gaps_secondary * normals_interface_secondary, axis=1) < threshold
    )

    penetration_node_primary = nodes_interface_primary_mask
    penetration_node_secondary = nodes_interface_secondary_mask

    penetration_node_primary[nodes_interface_primary_mask] = penetrates_secondary
    penetration_node_secondary[nodes_interface_secondary_mask] = penetrates_primary

    return penetration_node_primary, penetration_node_secondary

class ContactMechanicsInternodes(object):
    def __init__(
        self,
        dim,
        nodal_positions,
        nodal_displacements,
        surface_connectivity,
        nodes_candidate_primary,
        nodes_candidate_secondary,
        nodes_blocked,
        external_force,
        M,
        K,
        E,
        rbf=wendland_rbf,
    ):
        self.dim = dim

        # Nodal positions of all nodes included in mesh
        self.nodal_positions = nodal_positions.copy()
        self.surface_connectivity = surface_connectivity.copy()

        # Candidate nodes and their connectivity (2D: segments, 3D: triangles)
        self.nodes_candidate_primary = nodes_candidate_primary.copy()
        self.nodes_candidate_secondary = nodes_candidate_secondary.copy()
        self.connectivity_candidate_primary = remove_rows_without_all_items(
            self.surface_connectivity, self.nodes_candidate_primary
        )
        self.connectivity_candidate_secondary = remove_rows_without_all_items(
            self.surface_connectivity, self.nodes_candidate_secondary
        )

        # Nodes, positions, and corresponding dofs of primary/secondary interface
        self.nodes_interface_primary = nodes_candidate_primary.copy()
        self.nodes_interface_secondary = nodes_candidate_secondary.copy()
        self.positions_interface_primary = self.nodal_positions[
            self.nodes_interface_primary
        ]
        self.positions_interface_secondary = self.nodal_positions[
            self.nodes_interface_secondary
        ]
        self.dofs_interface_primary = nodes_to_dofs(
            self.nodes_interface_primary, dim=self.dim
        )
        self.dofs_interface_secondary = nodes_to_dofs(
            self.nodes_interface_secondary, dim=self.dim
        )

        # All dofs, blocked dofs, and free dofs
        self.dofs = np.arange(len(self.nodal_positions) * self.dim)
        self.dofs_blocked = self.dofs[nodes_blocked.ravel()]
        self.dofs_free = self.dofs[~nodes_blocked.ravel()]

        # Radial basis function and the corresponding radius parameters
        self.rbf = rbf
        self.rbf_radius_parameters_primary = None
        self.rbf_radius_parameters_secondary = None

        # Objects for use in the solution process
        self.scaling_factor = E  # Young's modulus
        self.external_force = external_force
        self.nodal_displacements = nodal_displacements
        self.M = sp.sparse.csr_matrix(M)  # Mass matrix
        self.K = sp.sparse.csr_matrix(K)  # Stiffness matrix

        self.R12 = None  # Interpolation matrix from secondary to primary
        self.R21 = None  # Interpolation matrix from primary to secondary

        self.B = None  # Block component of INTERNODES matrix
        self.B_tilde = None  # Block component of INTERNODES matrix

        self.internodes_matrix = None  # INTERNODES matrix
        self.force_term = None  # Force term b

    def _assemble_interpolation_matrices(self):
        """Assemble the interpolation matrices $R_{12}$ and $R_{21}$.

        Reference
        ---------
        [1], Page 53, Section 3.4, Bottom right
        """
        self.R12 = construct_interpolation_matrix(
            self.positions_interface_primary,
            self.positions_interface_secondary,
            self.rbf_radius_parameters_secondary,
            self.rbf,
        )
        self.R21 = construct_interpolation_matrix(
            self.positions_interface_secondary,
            self.positions_interface_primary,
            self.rbf_radius_parameters_primary,
            self.rbf,
        )

    def _assemble_B_matrices(self):
        """Assemble block components of INTERNODES matrix.

        Reference
        ---------
        [1], Page 54, Section 4, Bottom left
        """
        # Find indices in 'dofs_free' corresponding to interface of primary and secondary
        idx_interface_primary = np.argwhere(
            np.in1d(self.dofs_free, self.dofs_interface_primary)
        ).squeeze()
        idx_interface_secondary = np.argwhere(
            np.in1d(self.dofs_free, self.dofs_interface_secondary)
        ).squeeze()

        # Initialize the matrices
        self.B = sp.sparse.lil_matrix(
            (len(self.dofs_free), len(self.dofs_interface_primary)), dtype=np.float64
        )
        self.B_tilde = sp.sparse.lil_matrix(
            (len(self.dofs_interface_primary), len(self.dofs_free)), dtype=np.float64
        )

        # Fill in the blocks as described in the reference [1]
        M1 = self.M[np.ix_(self.dofs_interface_primary, self.dofs_interface_primary)]
        M2 = self.M[
            np.ix_(self.dofs_interface_secondary, self.dofs_interface_secondary)
        ]
        self.B[idx_interface_primary, :] = -M1
        self.B[idx_interface_secondary, :] = M2 * expand_to_dim(self.R21, dim=self.dim)

        self.B_tilde[:, idx_interface_primary] = sp.sparse.eye(
            len(self.dofs_interface_primary), dtype=np.float64, format="lil"
        )
        self.B_tilde[:, idx_interface_secondary] = -expand_to_dim(
            self.R12, dim=self.dim
        )

        # Convert sparsity format to csr
        self.B = sp.sparse.csr_matrix(self.B)
        self.B_tilde = sp.sparse.csr_matrix(self.B_tilde)

    def _assemble_internodes_matrix(self):
        """Assemble the INTERNODES matrix.

        Reference
        ---------
        [1], Page 54, Section 4, Equation 13
        """
        # Rescale stiffness matrix restricted to free dofs by Young's modulus
        K_free = self.K[np.ix_(self.dofs_free, self.dofs_free)] / self.scaling_factor
        self.internodes_matrix = sp.sparse.vstack(
            [
                sp.sparse.hstack([K_free, self.B]),
                sp.sparse.hstack(
                    [
                        self.B_tilde,
                        sp.sparse.csr_matrix(
                            (self.B_tilde.shape[0], self.B.shape[1]), dtype=np.float64
                        ),
                    ]
                ),
            ]
        )

    def _assemble_force_term(self):
        """Assemble the force term.

        Reference
        ---------
        [1], Page 54, Section 4, Equation 13
        """
        # Compute displacements for Dirichlet boundary condition offset
        dirichlet_offset = (
            self.K[np.ix_(self.dofs_free, self.dofs_blocked)]
            * self.nodal_displacements.ravel()[self.dofs_blocked]
        )

        # First component $f$ of force term is a adjusted by offset and rescaled
        virtual_force = (
            self.external_force.ravel()[self.dofs_free] - dirichlet_offset
        ) / self.scaling_factor

        # Nodal gaps $d$ between interpolated and true positions of primary nodes
        nodal_gaps = (
            self.R12 * self.positions_interface_secondary
            - self.positions_interface_primary
        )

        self.force_term = np.concatenate([virtual_force, nodal_gaps.ravel()])

    def define_interface(self):
        """Find contact/interface nodes while trying to satisfy the constraints.

        Reference
        ---------
        [1], Page 51, Section 2.3, Equation 2, 3, and 4
        """
        # Find the contact nodes and corresponding radial basis parameters
        (
            rbf_radius_parameters_primary,
            rbf_radius_parameters_secondary,
            interface_primary_mask,
            interface_secondary_mask,
        ) = find_interpolation_nodes(
            self.positions_interface_primary,
            self.positions_interface_secondary,
            rbf=wendland_rbf,
        )

        # Update current radius parameters of radial basis functions
        self.rbf_radius_parameters_primary = rbf_radius_parameters_primary
        self.rbf_radius_parameters_secondary = rbf_radius_parameters_secondary

        # Update interface node sets
        self.nodes_interface_primary = self.nodes_interface_primary[
            interface_primary_mask
        ]
        self.nodes_interface_secondary = self.nodes_interface_secondary[
            interface_secondary_mask
        ]

        # Update positions of interface node sets
        self.positions_interface_primary = self.positions_interface_primary[
            interface_primary_mask
        ]
        self.positions_interface_secondary = self.positions_interface_secondary[
            interface_secondary_mask
        ]

        # Update degrees of freedom with new node indices obtained
        self.dofs_interface_primary = nodes_to_dofs(
            self.nodes_interface_primary, dim=self.dim
        )
        self.dofs_interface_secondary = nodes_to_dofs(
            self.nodes_interface_secondary, dim=self.dim
        )

    def assemble_full_model(self):
        self._assemble_interpolation_matrices()
        self._assemble_B_matrices()
        self._assemble_internodes_matrix()
        self._assemble_force_term()

    def solve_direct(self):
        """Solve the INTERNODES system of equations.

        Returns
        -------
        displacements : np.ndarray
            Displacements $u$ of the free DOFs after the solve step.
        lambdas : np.ndarray
            Lagrange multipliers $\lambda$ after the solve step.

        Reference
        ---------
        [1], Page 54, Section 4, Equation 13
        """
        # Solve the internodes system
        x = sp.sparse.linalg.spsolve(self.internodes_matrix, self.force_term)

        # Fill in the computed displacements at the free dofs
        n_dofs_free = len(self.dofs_free)
        displacements = self.nodal_displacements.ravel()
        displacements[self.dofs_free] = x[:n_dofs_free]
        displacements = displacements.reshape((-1, self.dim))

        # Reshape and revert the scaling of the Lagrange multipliers $\lambda$
        lambdas = x[n_dofs_free:].reshape((-1, self.dim)) * self.scaling_factor

        return displacements, lambdas

    def update_interface(self, displacements, lambdas, return_changes=False):
        """Update the interfaces according to penetrations and tension

        Parameters
        ----------
        displacements : np.ndarray
            Displacements $u$ of the free DOFs after the solve step.
        lambdas : np.ndarray
            Lagrange multipliers $\lambda$ after the solve step.

        Returns
        -------
        converged : bool
            If all Lagrange multipliers negative and no penetration is detected.

        Reference
        ---------
        [2], Page 23, Algorithm 2, Lines 5-16
        """
        # Get new positions
        positions_new = self.nodal_positions + displacements

        # Compute the normals
        normals_candidate_primary = compute_normals(
            positions_new,
            self.nodes_candidate_primary,
            self.connectivity_candidate_primary,
            self.dim,
        )
        normals_candidate_secondary = compute_normals(
            positions_new,
            self.nodes_candidate_secondary,
            self.connectivity_candidate_secondary,
            self.dim,
        )
        nodes_interface_primary_mask = np.in1d(
            self.nodes_candidate_primary, self.nodes_interface_primary
        )
        nodes_interface_secondary_mask = np.in1d(
            self.nodes_candidate_secondary, self.nodes_interface_secondary
        )
        normals_interface_primary = normals_candidate_primary[
            nodes_interface_primary_mask
        ]
        normals_interface_secondary = normals_candidate_secondary[
            nodes_interface_secondary_mask
        ]

        # Interpolate the Lagrange multipliers of the secondary
        # [1], Page 53, Section 3.2, Equation (11)
        lambdas_secondary = -self.R21 * lambdas

        # Mark nodes with positive projected Lagrange multipliers for dumping
        positive_lambda_primary = (
            np.sum(lambdas * normals_interface_primary, axis=1) > 0
        )
        positive_lambda_secondary = (
            np.sum(lambdas_secondary * normals_interface_secondary, axis=1) > 0
        )
        nodes_to_dump_primary = self.nodes_interface_primary[positive_lambda_primary]
        nodes_to_dump_secondary = self.nodes_interface_secondary[
            positive_lambda_secondary
        ]

        # Add new nodes to primary and secondary
        positions_candidate_primary = positions_new[self.nodes_candidate_primary]
        positions_candidate_secondary = positions_new[self.nodes_candidate_secondary]
        nodes_penetration_primary, nodes_penetration_secondary = find_penetration_nodes(
            positions_candidate_primary,
            positions_candidate_secondary,
            normals_candidate_primary,
            normals_candidate_secondary,
        )
        nodes_to_add_primary = self.nodes_candidate_primary[
            np.logical_and(nodes_penetration_primary, ~nodes_interface_primary_mask)
        ]
        nodes_to_add_secondary = self.nodes_candidate_secondary[
            np.logical_and(nodes_penetration_secondary, ~nodes_interface_secondary_mask)
        ]

        # Add penetration nodes
        self.nodes_interface_primary = np.union1d(
            self.nodes_interface_primary, nodes_to_add_primary
        )
        self.nodes_interface_secondary = np.union1d(
            self.nodes_interface_secondary, nodes_to_add_secondary
        )

        # If projected Lagrange multipliers are all negative
        if len(nodes_to_dump_primary) + len(nodes_to_dump_secondary) == 0:

            # Convergence if all penetration nodes are already in interface
            if len(nodes_to_add_primary) + len(nodes_to_add_secondary) == 0:
                if return_changes:
                    return True, *[np.array([], dtype=int)] * 4
                return True

        else:
            # Dump nodes from primary and secondary
            self.nodes_interface_primary = np.setdiff1d(
                self.nodes_interface_primary, nodes_to_dump_primary
            )
            self.nodes_interface_secondary = np.setdiff1d(
                self.nodes_interface_secondary, nodes_to_dump_secondary
            )

        # Update interface positions according to new set of interface nodes
        self.positions_interface_primary = self.nodal_positions[
            self.nodes_interface_primary
        ]
        self.positions_interface_secondary = self.nodal_positions[
            self.nodes_interface_secondary
        ]

        if return_changes:
            return (
                False,
                nodes_to_add_primary,
                nodes_to_add_secondary,
                nodes_to_dump_primary,
                nodes_to_dump_secondary,
            )
        return False
