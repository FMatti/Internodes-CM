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
    return dim*np.repeat(nodes, dim) + np.tile(range(dim), len(nodes))

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
    return sp.sparse.kron(matrix, np.eye(dim), format='csr')

def remove_rows_without_items(matrix, items):
    """Remove rows from matrix that are not contained in items.

    Parameters
    ----------
    matrix : 2D np.ndarray
        A matrix.
    items : 1D np.ndarray or list
        A list with items.

    Returns
    -------
    matrix_removed : 2D np.ndarray
        The matrix where all rows without an item have been eliminated.

    Example
    -------
    >>> matrix = np.array([[1, 2],
    >>>                    [3, 4],
    >>>                    [5, 6],
    >>>                    [7, 8]])
    >>> items = [3, 8]
    >>> remove_rows_without_items(matrix, items)
    >>>   = array([[3, 4],
                   [7, 8]])
    """
    return matrix[np.isin(matrix, items).any(axis=1)]

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
    return (1 - delta)**4 * (1 + 4*delta) * (delta <= 1)

def wendland_rbf(distances, radiuses):
    """Evaluate Wendland radial basis function.
    
    Reference
    ---------
    [1], Page 49, Section 2.1, Bottom left
    """
    return wendland(distances / radiuses)

def _check_rbf_radius_parameters(positions, positions_ref, rbf_radius_parameters, c, C):
    distance_matrix_MM = sp.spatial.distance.cdist(positions, positions)
    distance_matrix_NM = sp.spatial.distance.cdist(positions_ref, positions)

    # Assert condition: [1], Page 51, Section 2.3, Equation 2
    np.fill_diagonal(distance_matrix_MM, np.inf)
    is_not_close_MM = distance_matrix_MM >= (c - 1e-10)*rbf_radius_parameters
    assert is_not_close_MM.all()

    # Assert condition: [1], Page 51, Section 2.3, Equation 3
    is_close_NM = distance_matrix_NM <= C*rbf_radius_parameters
    assert is_close_NM.any(axis=1).all()

    # Assert condition: [1], Page 51, Section 2.3, Equation 4
    n_supports_interpolation = np.sum(distance_matrix_MM < rbf_radius_parameters, axis=0)
    assert np.max(n_supports_interpolation) < 1 / wendland(c)

def compute_rbf_radius_parameters(positions, positions_ref, c=0.5, C=0.95, rbf=wendland):
    """Iteratively compute radius parameters for radial basis functions until
    invertibility conditions are satisfied (increase 'c' in every iteration).

    Parameters
    ----------
    positions : np.ndarray
        Positions of interpolation nodes.
    positions_ref : np.ndarray
        Positions of reference/evaluation nodes.
    c : float in (0, 1), default is 0.5
        [1], Page 51, Section 2.3, Equation 2
        (The empirical default value is given at the right on same page)
    C : float in (c, 1), default is 0.5
        [1], Page 51, Section 2.3, Equation 3
        (The empirical default value is given at the bottom right on same page)
    max_iter : int, default is 10
        Maximum number of iterations for finding suitable radius parameters.
    rbf : function, default is wendland
        The radial basis function being used.

    Returns
    -------
    rbf_radius_parameters : np.ndarray
        The radius parameters to use in the radial basis function.
    n_supports_reference : np.ndarray
        Number of radial basis functions to whose "shrunk" support (by factor C)
        each reference node belongs ([1], Page 51, Section 2.3, Equation 3)

    Reference
    ---------
    [1], Page 51, Section 2.3
    """
    # Compute distance matrices among nodes and from nodes to reference nodes:
    # distance_matrix(X, Y)[i, j] = ||x_i - y_j||
    distance_matrix_MM = sp.spatial.distance.cdist(positions, positions)
    distance_matrix_NM = sp.spatial.distance.cdist(positions_ref, positions)

    # Minimum distance of each node from closest distinct interpolation node
    np.fill_diagonal(distance_matrix_MM, np.inf)
    min_distance_MM = np.min(distance_matrix_MM, axis=0)
    
    # Extremely heuristic: Estimate for the minimum distance of the nodes
    # from the reference nodes (min_distance_MM/2 is supposed to be
    # representative of the distance an orthogonal projection of the
    # reference nodes onto the interface is away from the closest node)
    #d0 = 0.05
    #min_distance_NM = (d0**2 + (min_distance_MM/2)**2)**0.5
    # Minimum distance of each reference point from interpolation node
    #min_distance_NM = np.min(distance_matrix_NM, axis=0)
    # Take maximum of the two minimum distances from the closest nodes
    #rbf_radius_parameters = np.maximum(min_distance_MM, min_distance_NM)

    #rbf_radius_parameters = min_distance_MM / c
    #n_supports_reference = np.sum(distance_matrix_NM < C*rbf_radius_parameters, axis=0)

    # Iteratively increase parameter 'c' if necessary
    while True:

        rbf_radius_parameters = min_distance_MM / c

        # Enforce condition [1], Page 51, Section 2.3, Equation 2
        # (Also mentioned in [1], Page 51, Section 2.3, right center)
        #rbf_radius_parameters = np.minimum(candidate_radiuses, min_distance_MM / c)

        # Number of radial basis functions in whose support reference nodes are
        # belong (support is shrunk by C: [1], Page 51, Section 2.3, Equation 3)
        n_supports_reference = np.sum(distance_matrix_NM < C*rbf_radius_parameters, axis=0)

        # Number radial basis functions in whose support interpolation nodes are
        n_supports_interpolation = np.sum(distance_matrix_MM < rbf_radius_parameters, axis=0)

        # Check if criterion [1], Page 51, Section 2.3, Equation 4 is satisfied
        # for all interpolation nodes
        if np.max(n_supports_interpolation) < 1 / rbf(c):
            break

        # If criterion was not satisfied, increase c and reiterate the procedure
        c = (c + 1) / 2
        if c >= C:
            ValueError("Tried to increase c to", c, "which is larger than C.")

    return rbf_radius_parameters, n_supports_reference

def find_interface_nodes(positions_primary, positions_secondary, nodes_primary, nodes_secondary, rbf=wendland):
    """Find contact/interface nodes while trying to satisfy the constraints.

    Parameters
    ----------
    positions_primary: np.ndarray
        Positions of primary nodes (possibly only interface or candidates).
    positions_secondary: np.ndarray
        Positions of secondary nodes (possibly only interface or candidates).
    nodes_primary : np.ndarray
        Node indices of primary (possibly only interface or candidates).
    nodes_primary : np.ndarray
        Node indices of secondary (possibly only interface or candidates). 
    rbf : function, default is wendland
        The radial basis function being used.

    Reference
    ---------
    [1], Page 51, Section 2.3, Equation 2 and 3
    """

    while True:
        # Determine the radial basis function radius parameters and to how
        # many opposite radial basis function supports each node belongs
        rbf_radius_parameters_primary, n_supports_in_secondary = compute_rbf_radius_parameters(positions_primary, positions_secondary, rbf=rbf)
        rbf_radius_parameters_secondary, n_supports_in_primary = compute_rbf_radius_parameters(positions_secondary, positions_primary, rbf=rbf)
        
        # Update interface nodes by removing isolated nodes
        nodes_primary = nodes_primary[n_supports_in_secondary > 0]
        nodes_secondary = nodes_secondary[n_supports_in_primary > 0]
        
        # Update interface positions by removing isolated nodes
        positions_primary = positions_primary[n_supports_in_secondary > 0]
        positions_secondary = positions_secondary[n_supports_in_primary > 0]

        # Stop algorithm if for all secondary and primary interface nodes
        # [1], Page 51, Section 2.3, Equation 3 is satisfied
        if np.all(n_supports_in_secondary > 0) and np.all(n_supports_in_primary > 0):
            break

    # Verify conditions [1], Page 51, Section 2.3, (2), (3), and (4)
    #_check_rbf_radius_parameters(positions_primary, positions_secondary, rbf_radius_parameters_primary, c, C)

    return rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary

def construct_rbf_matrix(positions, positions_ref, radiuses_ref, rbf):
    """Construct radial basis matrix $\Phi_{NM}$.
    
    Reference
    ---------
    [1], Page 49, Section 2.1, Point 1
    """
    return rbf(sp.spatial.distance.cdist(positions.reshape(len(positions), -1), positions_ref.reshape(len(positions_ref), -1)), radiuses_ref)

def construct_interpolation_matrix(positions, positions_ref, radiuses_ref, rbf):
    """Construct interpolation matrix $\R_{NM}$.

    Reference
    ---------
    [1], Page 49, Section 2.1, Bottom right
    """
    # Construct radial basis function matrices $\Phi_{MM}$ and $\Phi_{NM}$
    rbf_matrix_MM = construct_rbf_matrix(positions_ref, positions_ref, radiuses_ref, rbf)
    rbf_matrix_NM = construct_rbf_matrix(positions, positions_ref, radiuses_ref, rbf)

    # Compute raw interpolation matrix without rescaling
    interpolation_matrix = np.linalg.solve(rbf_matrix_MM.T, rbf_matrix_NM.T).T

    # Compute the diagonal rescaling factors (diagonal entries of $D_{NN}$)
    rescaling_factors = np.sum(interpolation_matrix, axis=1)[:, np.newaxis]
    return sp.sparse.csr_matrix(interpolation_matrix / rescaling_factors)

class ContactMechanicsInternodes(object):

    def __init__(self, dim, mesh, model, name_candidate_interface_primary, name_candidate_interface_secondary, rbf=wendland_rbf, blocked_nodes_name=None):
        self.dim = dim
        self.mesh = mesh
        self.model = model
    
        # Candidate nodes for contact interface
        self.nodes_candidate_primary = mesh.getElementGroup(name_candidate_interface_primary).getNodeGroup().getNodes().ravel()
        self.nodes_candidate_secondary = mesh.getElementGroup(name_candidate_interface_secondary).getNodeGroup().getNodes().ravel()
        
        # Nodes, positions, and corresponding dofs of primary/secondary interface
        self.nodes_interface_primary = mesh.getElementGroup(name_candidate_interface_primary).getNodeGroup().getNodes().ravel()
        self.nodes_interface_secondary = mesh.getElementGroup(name_candidate_interface_secondary).getNodeGroup().getNodes().ravel()
        self.positions_interface_primary = mesh.getNodes()[self.nodes_interface_primary]
        self.positions_interface_secondary = mesh.getNodes()[self.nodes_interface_secondary]
        self.dofs_interface_primary = nodes_to_dofs(self.nodes_interface_primary, dim=self.dim)
        self.dofs_interface_secondary = nodes_to_dofs(self.nodes_interface_secondary, dim=self.dim)

        # All dofs, blocked dofs, and free dofs
        self.dofs = np.arange(mesh.getNbNodes()*self.dim)
        blocked_dofs = model.getBlockedDOFs()
        if blocked_nodes_name:
            blocked_dofs[mesh.getElementGroup(blocked_nodes_name).getNodeGroup().getNodes()] = True
        self.dofs_blocked = self.dofs[blocked_dofs.ravel()]
        self.dofs_free = self.dofs[~blocked_dofs.ravel()]

        # Connectivity of the model (line segments and triangular elements)
        triangles = mesh.getConnectivity(aka._triangle_3)
        if self.dim == 2:
            self.connectivity = triangles[:, np.array([0, 1, 0, 2, 1, 2])].reshape(-1, 2)
        elif self.dim == 3:
            self.connectivity = triangles

        # Radial basis function and the corresponding radius parameters
        self.rbf = rbf
        self.rbf_radius_parameters_primary = None
        self.rbf_radius_parameters_secondary = None
        
        # Objects for use in the solution process
        self.rescaling_factor = 30e+9  # TODO: Tunable, Young's modulus
        self.K = None  # Stiffness matrix

        self.R12 = None  # Interpolation matrix from secondary to primary
        self.R21 = None  # Interpolation matrix from primary to secondary

        self.M1 = None  # Interface mass matrix of primary
        self.M2 = None  # Interface mass matrix of secondary

        self.B = None  # Block component of INTERNODES matrix
        self.B_tilde = None  # Block component of INTERNODES matrix

        self.internodes_matrix = None  # INTERNODES matrix
        self.force_term = None  # Force term b

    def find_interface_nodes(self):
        """Find contact/interface nodes while trying to satisfy the constraints.
       
        Reference
        ---------
        [1], Page 51, Section 2.3, Equation 2, 3, and 4
        """
        # Find the contact nodes and corresponding radial basis parameters
        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_interface_primary, nodes_interface_secondary = find_interface_nodes(self.positions_interface_primary, self.positions_interface_secondary, self.nodes_interface_primary, self.nodes_interface_secondary, rbf=wendland)
        
        # Update current radius parameters of radial basis functions
        self.rbf_radius_parameters_primary = rbf_radius_parameters_primary
        self.rbf_radius_parameters_secondary = rbf_radius_parameters_secondary
        
        # Update positions with new node indices obtained
        self.positions_interface_primary = self.positions_interface_primary[np.in1d(self.nodes_interface_primary, nodes_interface_primary)]
        self.positions_interface_secondary = self.positions_interface_secondary[np.in1d(self.nodes_interface_secondary, nodes_interface_secondary)]

        # Update node indices new node indices obtained
        self.nodes_interface_primary = nodes_interface_primary
        self.nodes_interface_secondary = nodes_interface_secondary

        # Update degrees of freedom with new node indices obtained
        self.dofs_interface_primary = nodes_to_dofs(self.nodes_interface_primary, dim=self.dim)
        self.dofs_interface_secondary = nodes_to_dofs(self.nodes_interface_secondary, dim=self.dim)

    def assemble_interpolation_matrices(self):
        """Assemble the interpolation matrices $R_{12}$ and $R_{21}$.
        
        Reference
        ---------
        [1], Page 53, Section 3.4, Bottom right
        """
        self.R12 = construct_interpolation_matrix(self.positions_interface_primary, self.positions_interface_secondary, self.rbf_radius_parameters_secondary, self.rbf)
        self.R21 = construct_interpolation_matrix(self.positions_interface_secondary, self.positions_interface_primary, self.rbf_radius_parameters_primary, self.rbf)

    def assemble_interface_mass_matrices(self):
        """Assemble the interface mass matrices $M_1, M_2$.

        Reference
        ---------
        [1], Page 53, Section 3.4, Equation 12
        """
        self.model.assembleMass()
        mass_matrix = self.model.getDOFManager().getMatrix('M')
        mass_matrix = aka.AkantuSparseMatrix(mass_matrix).toarray()

        # Slice global mass matrix into primary and secondary matrices by dofs
        self.M1 = mass_matrix[np.ix_(self.dofs_interface_primary, self.dofs_interface_primary)]
        self.M2 = mass_matrix[np.ix_(self.dofs_interface_secondary, self.dofs_interface_secondary)]

        # Convert matrices to sparse csr format
        self.M1 = sp.sparse.csr_matrix(self.M1)
        self.M2 = sp.sparse.csr_matrix(self.M2)

    def assemble_stiffness_matrix(self):
        """Assemble the global stiffness matrix $K$."""
        self.model.assembleStiffnessMatrix()
        self.K = self.model.getDOFManager().getMatrix('K')
        self.K = aka.AkantuSparseMatrix(self.K).toarray()
        self.K = sp.sparse.csr_matrix(self.K)

    def assemble_B_matrices(self):
        """Assemble block components of INTERNODES matrix.

        Reference
        ---------
        [1], Page 54, Section 4, Bottom left
        """
        # Find indices in 'dofs_free' corresponding to interface of primary and secondary
        idx_interface_primary = np.argwhere(np.in1d(self.dofs_free, self.dofs_interface_primary)).squeeze()
        idx_interface_secondary = np.argwhere(np.in1d(self.dofs_free, self.dofs_interface_secondary)).squeeze()

        # Initialize the matrices
        self.B = sp.sparse.lil_matrix((len(self.dofs_free), len(self.dofs_interface_primary)), dtype=np.float64)
        self.B_tilde = sp.sparse.lil_matrix((len(self.dofs_interface_primary), len(self.dofs_free)), dtype=np.float64)

        # Fill in the blocks as described in the reference [1]
        self.B[idx_interface_primary, :] = - self.M1
        self.B[idx_interface_secondary, :] = self.M2 * expand_to_dim(self.R21, dim=self.dim)

        self.B_tilde[:, idx_interface_primary] = sp.sparse.eye(len(self.dofs_interface_primary), dtype=np.float64, format="lil")
        self.B_tilde[:, idx_interface_secondary] = - expand_to_dim(self.R12, dim=self.dim)

        # Convert sparsity format to csr
        self.B = sp.sparse.csr_matrix(self.B)
        self.B_tilde = sp.sparse.csr_matrix(self.B_tilde)

    def assemble_internodes_matrix(self):
        """Assemble the INTERNODES matrix.

        Reference
        ---------
        [1], Page 54, Section 4, Equation 13
        """
        # Rescale stiffness matrix restricted to free dofs by Young's modulus
        K_free = self.K[np.ix_(self.dofs_free, self.dofs_free)] / self.rescaling_factor
        self.internodes_matrix = sp.sparse.vstack([
            sp.sparse.hstack([K_free, self.B]),
            sp.sparse.hstack([self.B_tilde, sp.sparse.csr_matrix(
                (self.B_tilde.shape[0], self.B.shape[1]), dtype=np.float64
            )])
        ])

    def assemble_force_term(self, f_free, displacements):
        """Assemble the force term.
        
        Parameters
        ----------
        f_free : np.ndarray
            Force applied to free DOFs.
        displacements : np.ndarray
            Displacements of the free DOFs.

        Reference
        ---------
        [1], Page 54, Section 4, Equation 13
        """
        # Compute displacements for Dirichlet boundary condition offset
        dirichlet_offset = self.K[np.ix_(self.dofs_free, self.dofs_blocked)] * displacements[self.dofs_blocked]
        
        # First component $f$ of force term is a adjusted by offset and rescaled
        virtual_force = (f_free - dirichlet_offset) / self.rescaling_factor

        # Nodal gaps $d$ between interpolated and true positions of primary nodes
        nodal_gaps = self.R12 * self.positions_interface_secondary - self.positions_interface_primary
    
        self.force_term = np.concatenate([virtual_force, nodal_gaps.ravel()])

    def solve_direct(self, displacements):
        """Solve the INTERNODES system of equations.
        
        Parameters
        ----------
        displacements : np.ndarray
            Displacements of the free DOFs before the solve step 

        Returns
        -------
        positions_new : np.ndarray
            Positions after the solve step 
        displacements : np.ndarray
            Displacements $u$ of the free DOFs after the solve step 
        lambdas : np.ndarray
            Lagrange multipliers $\lambda$ after the solve step
        
        Reference
        ---------
        [1], Page 54, Section 4, Equation 13
        """
        x = sp.sparse.linalg.spsolve(self.internodes_matrix, self.force_term)

        # FIll in the computed displacements at the free dofs
        displacements[self.dofs_free] = x[:len(self.dofs_free)]
        positions_new = self.mesh.getNodes() + displacements.reshape((-1, self.dim))

        # Reshape and revert the rescaling of the Lagrange multipliers $\lambda$
        lambdas = x[len(self.dofs_free):].reshape((-1, self.dim)) * self.rescaling_factor

        return positions_new, displacements, lambdas

    def update_interface(self, positions_new, lambdas):
        """Update the interfaces according to penetrations and tension
        
        Parameters
        ----------
        positions_new : np.ndarray
            New positions of nodes after solving the INTERNODES system.
        lambdas : np.ndarray
            Lagrange multipliers obtained from solving the INTERNODES system.

        Returns
        -------
        converged : bool
            If all Lagrange multipliers negative and no penetration is detected

        Reference
        ---------
        [2], Page 23, Algorithm 2, Lines 5-16
        """
        # Get connectivities (triangle/segment indices) between candidate nodes
        connectivity_candidate_primary = remove_rows_without_all_items(self.connectivity, self.nodes_candidate_primary)
        connectivity_candidate_secondary = remove_rows_without_all_items(self.connectivity, self.nodes_candidate_secondary)
        
        # Compute the normals
        normals_candidate_primary = self.compute_normals(positions_new, self.nodes_candidate_primary, connectivity_candidate_primary)
        normals_candidate_secondary = self.compute_normals(positions_new, self.nodes_candidate_secondary, connectivity_candidate_secondary)
        normals_interface_primary = normals_candidate_primary[np.in1d(self.nodes_candidate_primary, self.nodes_interface_primary)]
        normals_interface_secondary = normals_candidate_secondary[np.in1d(self.nodes_candidate_secondary, self.nodes_interface_secondary)]

        # Interpolate the Lagrange multipliers of the secondary
        # [1], Page 53, Section 3.2, Equation (11)
        lambdas_secondary = - self.R21 * lambdas

        # Mark nodes with positive projected Lagrange multipliers for dumping
        positive_lambda_proj_primary = np.sum(lambdas * normals_interface_primary, axis=1) > 0
        positive_lambda_proj_secondary = np.sum(lambdas_secondary * normals_interface_secondary, axis=1) > 0
        nodes_to_dump_primary = self.nodes_interface_primary[positive_lambda_proj_primary]
        nodes_to_dump_secondary = self.nodes_interface_secondary[positive_lambda_proj_secondary]

        # Add new nodes to primary and secondary
        nodes_to_add_primary, nodes_to_add_secondary = self.detect_penetration_nodes(positions_new, normals_candidate_primary, normals_candidate_secondary)
        self.nodes_interface_primary = np.union1d(self.nodes_interface_primary, nodes_to_add_primary)
        self.nodes_interface_secondary = np.union1d(self.nodes_interface_secondary, nodes_to_add_secondary)

        # If projected Lagrange multipliers are all negative
        if not (positive_lambda_proj_primary.any() or positive_lambda_proj_secondary.any()):

            # Convergence if all penetration nodes are already in interface
            if (np.all(np.in1d(nodes_to_add_primary, self.nodes_interface_primary))
             or np.all(np.in1d(nodes_to_add_secondary, self.nodes_interface_secondary))):
                return True

        else:
            # Dump nodes from primary and secondary
            self.nodes_interface_primary = np.setdiff1d(self.nodes_interface_primary, nodes_to_dump_primary)
            self.nodes_interface_secondary = np.setdiff1d(self.nodes_interface_secondary, nodes_to_dump_secondary)

        # Update interface positions according to new set of interface nodes
        self.positions_interface_primary = self.mesh.getNodes()[self.nodes_interface_primary]
        self.positions_interface_secondary = self.mesh.getNodes()[self.nodes_interface_secondary]

        return False

    def compute_normals(self, positions_new, nodes, connectivity):
        """Compute the normals.

        Parameters
        ----------
        positions_new : np.ndarray
            New positions of nodes after solving the INTERNODES system.
        nodes : np.ndarray 
            Nodes for which the normals are computed for.
        segments : np.ndarray
            Connectivity of the line segments in the mesh.
        triangles : np.ndarray
            Connectivity of the triangular elements in the mesh.
    
        Returns
        -------
        normals_avg : np.ndarray
        """

        # Determine boundary nodes that aren't interface elements
        #nodes_candidate_noninterface = np.setdiff1d(all_triangles_candidate, nodes_candidate)

        # Compute tangents corresponding to the elements at new positions
        if self.dim == 2:
            tangent1 = positions_new[connectivity[:, 1]] - positions_new[connectivity[:, 0]]
            tangent2 = [0, 0, 1]
        elif self.dim == 3:
            tangent1 = positions_new[connectivity[:, 1]] - positions_new[connectivity[:, 0]]
            tangent2 = positions_new[connectivity[:, 2]] - positions_new[connectivity[:, 0]]

        # Compute normal vectors
        normals = np.cross(tangent1, tangent2)[:, :self.dim]

        normals_avg = np.zeros((len(nodes), self.dim))
        for j, node in enumerate(nodes):
            id = np.isin(connectivity, node).any(axis=1)

            # Compute average surface normal for each candidate node
            normals_avg[j] = np.sum(normals[id] / np.linalg.norm(normals[id], axis=1)[:, np.newaxis], axis=0)

            """
            # This below steps are unnecessary if mesh oriented correctly!
            step_size=1e-3

            # Compute average normal on interface boundary
            normal_avg = np.sum(normals[id] / np.linalg.norm(normals[id], axis=1)[:, np.newaxis], axis=0)

            # Extend positions by +/- a step size in the normal direction
            positions_plus = positions_new[node] + step_size*normal_avg
            positions_minus = positions_new[node] - step_size*normal_avg

            # Find minimum distance of new positions of non-interface nodes to WHAT?!?!
            min_plus = sp.spatial.distance.cdist(positions_new[nodes_candidate_noninterface], positions_plus.reshape((-1, self.dim))).min()
            min_minus = sp.spatial.distance.cdist(positions_new[nodes_candidate_noninterface], positions_minus.reshape((-1, self.dim))).min()

            # Change sign of the average normal vector
            normals_avg[j] = normal_avg * (1 if min_plus > min_minus else -1)
            """

        normals_avg /= np.linalg.norm(normals_avg, axis=1)[:, np.newaxis]

        return normals_avg

    def detect_penetration_nodes(self, positions_new, normals_candidates_primary, normals_candidates_secondary, tolerance=0.9, mesh_size=0.1):
        """Detect the nodes on secondary and primary interface that penetrate.
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

        # Update positions of primary and secondary nodes along interface
        positions_candidate_primary = positions_new[self.nodes_candidate_primary]
        positions_candidate_secondary = positions_new[self.nodes_candidate_secondary]

        # Find contact nodes with contact detection algorithm
        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_interface_primary, nodes_interface_secondary = find_interface_nodes(positions_candidate_primary, positions_candidate_secondary, self.nodes_candidate_primary, self.nodes_candidate_secondary)

        # Update positions of primary and secondary nodes along interface
        positions_interface_primary = positions_new[nodes_interface_primary]
        positions_interface_secondary = positions_new[nodes_interface_secondary]

        # Update normals of primary and secondary nodes along interface
        normals_interface_primary = normals_candidates_primary[np.in1d(self.nodes_candidate_primary, nodes_interface_primary)]
        normals_interface_secondary = normals_candidates_secondary[np.in1d(self.nodes_candidate_secondary, nodes_interface_secondary)]

        R12 = construct_interpolation_matrix(positions_interface_primary, positions_interface_secondary, rbf_radius_parameters_secondary, self.rbf)
        R21 = construct_interpolation_matrix(positions_interface_secondary, positions_interface_primary, rbf_radius_parameters_primary, self.rbf)
  
        # Determine the size of the nodal gaps on interface
        nodal_gaps_primary = R12 * positions_interface_secondary - positions_interface_primary
        nodal_gaps_secondary = R21 * positions_interface_primary - positions_interface_secondary

        # Penetration if projected nodal gap onto normal is sufficiently small
        threshold = - tolerance * mesh_size
        penetrates_secondary = np.sum(nodal_gaps_primary * normals_interface_primary, axis=1) < threshold
        penetrates_primary = np.sum(nodal_gaps_secondary * normals_interface_secondary, axis=1) < threshold

        # Add the nodes where penetration is observed to the interface
        nodes_penetration_primary = nodes_interface_primary[penetrates_secondary]
        nodes_penetration_secondary = nodes_interface_secondary[penetrates_primary]

        return nodes_penetration_primary, nodes_penetration_secondary
