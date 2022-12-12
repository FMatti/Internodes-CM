import pytest
import numpy as np
import scipy as sp
import matplotlib.tri as tri

from helper import get_theoretical_normal_displacement, get_theoretical_contact_radius
from contact_mechanics_internodes import *

def _get_uniform_interface_grid(dim, graph, x_min=-1, x_max=1, n_grid=10):
    linspace = np.linspace(x_min, x_max, n_grid)
    meshgrid = np.meshgrid(*[linspace]*(dim-1))
    grid = np.c_[[meshgrid[i].ravel() for i in range(dim-1)]]
    mesh = np.c_[grid.T, graph(grid)]
    return mesh

def _get_random_interface_grid(dim, graph, x_min=-1, x_max=1, n_grid=10, seed=0):
    np.random.seed(0)
    grid = np.random.uniform(x_min, x_max, (dim-1, n_grid**(dim-1)))
    mesh = np.c_[grid.T, graph(grid)]
    return mesh

def _check_rbf_radius_conditions(positions, rbf_radius_parameters, c_min=0.49, c_max=0.95):
    # Check condition: [1], Page 51, Section 2.3, Equation 2
    dist_MM = sp.spatial.distance.cdist(positions, positions)
    np.fill_diagonal(dist_MM, np.inf)
    min_ratio = np.min(dist_MM / rbf_radius_parameters, axis=1)
    np.testing.assert_array_less(c_min*np.ones_like(min_ratio),  min_ratio)

    # Check condition: [1], Page 51, Section 2.3, Equation 4
    n_supports = np.sum(dist_MM < rbf_radius_parameters, axis=0)
    n_supports_max = np.ones_like(n_supports) / wendland_rbf(c_max, 1)
    np.testing.assert_array_less(n_supports,  n_supports_max)

def _check_node_contained_in_support(positions, positions_ref, rbf_radius_parameters, C_max=0.99):
    # Check condition: [1], Page 51, Section 2.3, Equation 3
    dist_MN = sp.spatial.distance.cdist(positions_ref, positions)
    min_ratio = np.min(dist_MN / rbf_radius_parameters, axis=1)
    np.testing.assert_array_less(min_ratio, C_max*np.ones_like(min_ratio))

def test_compute_rbf_radius_parameters_regular():

    for dim in [2, 3]:

        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
        positions_primary = _get_uniform_interface_grid(dim, g_zeros, n_grid=10)
        positions_secondary = _get_uniform_interface_grid(dim, g_parabolic, n_grid=10)

        rbf_radius_parameters_primary = compute_rbf_radius_parameters(positions_primary)
        _check_rbf_radius_conditions(positions_primary, rbf_radius_parameters_primary)

        rbf_radius_parameters_secondary = compute_rbf_radius_parameters(positions_secondary)
        _check_rbf_radius_conditions(positions_secondary, rbf_radius_parameters_secondary)

def test_compute_rbf_radius_parameters_random():

    for dim in [2, 3]:
        
        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
        positions_primary = _get_random_interface_grid(dim, g_zeros, n_grid=10)
        positions_secondary = _get_random_interface_grid(dim, g_parabolic, n_grid=10)

        rbf_radius_parameters_primary = compute_rbf_radius_parameters(positions_primary)
        _check_rbf_radius_conditions(positions_primary, rbf_radius_parameters_primary)

        rbf_radius_parameters_secondary = compute_rbf_radius_parameters(positions_secondary)
        _check_rbf_radius_conditions(positions_secondary, rbf_radius_parameters_secondary)

def test_find_interpolation_nodes_regular():
    
    for dim in [2, 3]:

        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
        positions_primary = _get_uniform_interface_grid(dim, g_zeros, n_grid=10)
        positions_secondary = _get_uniform_interface_grid(dim, g_parabolic, n_grid=10)

        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary = find_interpolation_nodes(positions_primary, positions_secondary, rbf=wendland_rbf, C=0.95)
        _check_node_contained_in_support(positions_primary[nodes_primary], positions_secondary[nodes_secondary], rbf_radius_parameters_primary)
        _check_node_contained_in_support(positions_secondary[nodes_secondary], positions_primary[nodes_primary], rbf_radius_parameters_secondary)

def test_find_interpolation_nodes_random():
    
    for dim in [2, 3]:

        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
        positions_primary = _get_random_interface_grid(dim, g_zeros, n_grid=10)
        positions_secondary = _get_random_interface_grid(dim, g_parabolic, n_grid=10)

        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary = find_interpolation_nodes(positions_primary, positions_secondary, rbf=wendland_rbf, C=0.95)
        _check_node_contained_in_support(positions_primary[nodes_primary], positions_secondary[nodes_secondary], rbf_radius_parameters_primary)
        _check_node_contained_in_support(positions_secondary[nodes_secondary], positions_primary[nodes_primary], rbf_radius_parameters_secondary)

def test_construct_interpolation_matrix_constant():

    # Test if constant function is interpolated exactly
    positions = np.array([1, 2, 3])
    positions_ref = np.array([1.5, 2.5])
    radiuses_ref = np.array([1, 1])
    rbf = wendland_rbf
    R = construct_interpolation_matrix(positions, positions_ref, radiuses_ref, rbf)
    f = lambda x : np.ones_like(x)
    np.testing.assert_allclose(R * f(positions_ref), f(positions), rtol=1e-10)

def test_construct_interpolation_matrix_polynomial():

    # Test if polynomial is interpolated well enough
    positions = np.linspace(0.05, 0.95, 10)
    positions_ref = np.linspace(0, 1, 11)
    radiuses_ref = 1e-1*np.ones_like(positions_ref)
    rbf = wendland_rbf
    R = construct_interpolation_matrix(positions, positions_ref, radiuses_ref, rbf)
    f = lambda x : x**3 + x**2 + x + 1
    np.testing.assert_allclose(R * f(positions_ref), f(positions), rtol=1e-2)

def test_construct_gap_function_interpolation_regular():

     for dim in [2, 3]:

        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
        positions_primary = _get_uniform_interface_grid(dim, g_zeros, n_grid=10)
        positions_secondary = _get_uniform_interface_grid(dim, g_parabolic, n_grid=10)

        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary = find_interpolation_nodes(positions_primary, positions_secondary, rbf=wendland_rbf, C=0.95)
        R12 = construct_interpolation_matrix(positions_primary[nodes_primary], positions_secondary[nodes_secondary], rbf_radius_parameters_secondary, rbf=wendland_rbf)
        R21 = construct_interpolation_matrix(positions_secondary[nodes_secondary], positions_primary[nodes_primary], rbf_radius_parameters_primary, rbf=wendland_rbf)
        positions_interpolated_primary = R12 * positions_secondary[nodes_secondary]
        positions_interpolated_secondary = R21 * positions_primary[nodes_primary]

        np.testing.assert_allclose(g_parabolic(positions_interpolated_primary[:, :-1].T), positions_interpolated_primary[:, -1], atol=2e-2)
        np.testing.assert_allclose(g_zeros(positions_interpolated_secondary[:, :-1].T), positions_interpolated_secondary[:, -1], atol=2e-2)

def test_construct_gap_function_interpolation_random():

     for dim in [2, 3]:

        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
        positions_primary = _get_random_interface_grid(dim, g_zeros, n_grid=10)
        positions_secondary = _get_random_interface_grid(dim, g_parabolic, n_grid=10)
        nodes_primary = np.arange(len(positions_primary))
        nodes_secondary = np.arange(len(positions_secondary))

        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary = find_interpolation_nodes(positions_primary, positions_secondary, rbf=wendland_rbf, C=0.95)
        R12 = construct_interpolation_matrix(positions_primary[nodes_primary], positions_secondary[nodes_secondary], rbf_radius_parameters_secondary, rbf=wendland_rbf)
        R21 = construct_interpolation_matrix(positions_secondary[nodes_secondary], positions_primary[nodes_primary], rbf_radius_parameters_primary, rbf=wendland_rbf)
        positions_interpolated_primary = R12 * positions_secondary[nodes_secondary]
        positions_interpolated_secondary = R21 * positions_primary[nodes_primary]

        np.testing.assert_allclose(g_parabolic(positions_interpolated_primary[:, :-1].T), positions_interpolated_primary[:, -1], atol=2e-2)
        np.testing.assert_allclose(g_zeros(positions_interpolated_secondary[:, :-1].T), positions_interpolated_secondary[:, -1], atol=2e-2)

def test_compute_normals_plane():

    n_grid = 50

    for dim in [2, 3]:

        connectivity = np.arange(n_grid)
        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic_grad = lambda x: -np.zeros_like(x)
        nodal_positions = _get_uniform_interface_grid(dim, g_zeros, n_grid=n_grid)

        if dim == 2:
            nodes = np.argsort(nodal_positions[:, 0])
            connectivity = np.c_[nodes[1:], nodes[:-1]]
        if dim == 3:
            connectivity = tri.Triangulation(nodal_positions[:, 0], nodal_positions[:, 1]).triangles

        normals = compute_normals(nodal_positions, np.arange(n_grid**(dim-1)), connectivity, dim)

        normals_ex = np.c_[g_parabolic_grad(nodal_positions[:, :-1]), np.ones_like(nodal_positions[:, -1])]
        normals_ex /= np.linalg.norm(normals_ex, axis=1)[:, np.newaxis]

        if dim == 2:
            inner_nodes = ~(np.abs(nodal_positions[:, 0]) == 1)
        elif dim == 3:
            inner_nodes = ~np.logical_or(np.abs(nodal_positions[:, 0]) == 1, np.abs(nodal_positions[:, 1] == 1))

        np.testing.assert_allclose(normals_ex[inner_nodes], normals[inner_nodes], atol=1e-10)

def test_compute_normals_parabolic():

    n_grid = 50

    for dim in [2, 3]:

        connectivity = np.arange(n_grid)
        g_parabolic = lambda x: np.sum(x**2, axis=0)
        g_parabolic_grad = lambda x: -2*x
        nodal_positions = _get_uniform_interface_grid(dim, g_parabolic, n_grid=n_grid)

        if dim == 2:
            nodes = np.argsort(nodal_positions[:, 0])
            connectivity = np.c_[nodes[1:], nodes[:-1]]
        if dim == 3:
            connectivity = tri.Triangulation(nodal_positions[:, 0], nodal_positions[:, 1]).triangles

        normals = compute_normals(nodal_positions, np.arange(n_grid**(dim-1)), connectivity, dim)

        normals_ex = np.c_[g_parabolic_grad(nodal_positions[:, :-1]), np.ones_like(nodal_positions[:, -1])]
        normals_ex /= np.linalg.norm(normals_ex, axis=1)[:, np.newaxis]

        if dim == 2:
            inner_nodes = ~(np.abs(nodal_positions[:, 0]) == 1)
        elif dim == 3:
            inner_nodes = ~np.logical_or(np.abs(nodal_positions[:, 0]) == 1, np.abs(nodal_positions[:, 1] == 1))

        np.testing.assert_allclose(normals_ex[inner_nodes], normals[inner_nodes], atol=2e-2)

def test_contact_problem_3d():

    # Set up contact problem
    mesh_file = 'mesh/contact3d_sphere.msh'
    material_file = 'material/material.dat'
    spatial_dimension = 3
    aka.parseInput(material_file)

    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._implicit_dynamic)

    # Apply boundary conditions
    model.applyBC(aka.FixedValue(0., aka._x), 'primary_fixed')
    model.applyBC(aka.FixedValue(0., aka._y), 'primary_fixed')
    model.applyBC(aka.FixedValue(0., aka._z), 'primary_fixed')
    model.applyBC(aka.FixedValue(0., aka._x), 'secondary_fixed')
    model.applyBC(aka.FixedValue(0., aka._z), 'secondary_fixed')

    # Get positions of all nodes, surface connectivity and candidate nodes
    nodal_positions = mesh.getNodes()
    surface_connectivity = mesh.getConnectivity(aka._triangle_3)
    nodes_candidate_primary = mesh.getElementGroup('primary_candidates').getNodeGroup().getNodes().ravel()
    nodes_candidate_secondary = mesh.getElementGroup('secondary_candidates').getNodeGroup().getNodes().ravel()
    external_force = model.getExternalForce()
    nodal_displacements = model.getDisplacement()
    nodes_blocked = model.getBlockedDOFs()

    model.assembleMass()
    M = aka.AkantuSparseMatrix(model.getDOFManager().getMatrix('M')).toarray()

    model.assembleStiffnessMatrix()
    K = aka.AkantuSparseMatrix(model.getDOFManager().getMatrix('K')).toarray()

    E = model.getMaterial(0).getReal("E")

    # Run internodes algorithm
    d0 = 0.05
    R = 0.5

    E = model.getMaterial(0).getReal("E")
    nu = model.getMaterial(0).getReal("nu")

    d_list = np.linspace(0.05, 0.2, 4)
    a_list = np.empty_like(d_list)
    u_list = np.empty_like(d_list)

    for j, d in enumerate(d_list):

        model.applyBC(aka.FixedValue(-d+d0, aka._y), 'secondary_fixed')
        nodal_displacements = model.getDisplacement()

        internodes_model = ContactMechanicsInternodes(spatial_dimension, nodal_positions, nodal_displacements, surface_connectivity, nodes_candidate_primary, nodes_candidate_secondary, nodes_blocked, external_force, M, K, E)

        max_iter = 10
        for i in range(max_iter):
            # Find the interface nodes
            internodes_model.define_interface()

            # Assemble model
            internodes_model.assemble_full_model()

            # Solve model
            displacements, lambdas = internodes_model.solve_direct()

            # Update the interface nodes and check if it converged
            converged = internodes_model.update_interface(displacements, lambdas)

            if converged:
                break

        assert i+1 < max_iter

        positions = internodes_model.nodal_positions + displacements
        positions_interface_secondary = positions[internodes_model.nodes_interface_secondary]

        a_list[j] = np.max(sp.spatial.distance.cdist(positions_interface_secondary, positions_interface_secondary))/2
        u_list[j] = np.min(positions_interface_secondary[:, 1])

    # Check if normal displacement corresponds to theoretical expectation
    np.testing.assert_allclose(get_theoretical_normal_displacement(R, d_list, E, nu), u_list, atol=5e-3)

    # Check if contact radius corresponds to theoretical expectation
    np.testing.assert_allclose(get_theoretical_contact_radius(R, d_list), a_list, atol=5e-2)

def test_contact_problem_2d():
    """
    # Set up contact problem
    mesh_file = 'mesh/contact2d_circle.msh'
    material_file = 'material/material.dat'
    spatial_dimension = 2
    aka.parseInput(material_file)

    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._implicit_dynamic)

    # Apply boundary conditions
    model.applyBC(aka.FixedValue(0., aka._x), 'primary_fixed')
    model.applyBC(aka.FixedValue(0., aka._y), 'primary_fixed')
    model.applyBC(aka.FixedValue(0., aka._z), 'primary_fixed')
    model.applyBC(aka.FixedValue(0., aka._x), 'secondary_fixed')
    model.applyBC(aka.FixedValue(0., aka._z), 'secondary_fixed')

    # Get positions of all nodes, surface connectivity and candidate nodes
    nodal_positions = mesh.getNodes()
    surface_connectivity = mesh.getConnectivity(aka._segment_2)
    nodes_candidate_primary = mesh.getElementGroup('primary_candidates').getNodeGroup().getNodes().ravel()
    nodes_candidate_secondary = mesh.getElementGroup('secondary_candidates').getNodeGroup().getNodes().ravel()
    external_force = model.getExternalForce()
    nodal_displacements = model.getDisplacement()
    nodes_blocked = model.getBlockedDOFs()

    model.assembleMass()
    M = aka.AkantuSparseMatrix(model.getDOFManager().getMatrix('M')).toarray()

    model.assembleStiffnessMatrix()
    K = aka.AkantuSparseMatrix(model.getDOFManager().getMatrix('K')).toarray()

    E = model.getMaterial(0).getReal("E")

    # Run internodes algorithm
    d0 = 0.05
    R = 0.5

    E = model.getMaterial(0).getReal("E")
    nu = model.getMaterial(0).getReal("nu")

    d_list = np.linspace(0.05, 0.2, 4)
    a_list = np.empty_like(d_list)
    u_list = np.empty_like(d_list)

    for j, d in enumerate(d_list):

        model.applyBC(aka.FixedValue(-d+d0, aka._y), 'secondary_fixed')
        nodal_displacements = model.getDisplacement()

        internodes_model = ContactMechanicsInternodes(spatial_dimension, nodal_positions, nodal_displacements, surface_connectivity, nodes_candidate_primary, nodes_candidate_secondary, nodes_blocked, external_force, M, K, E)

        max_iter = 10
        for i in range(max_iter):
            # Find the interface nodes
            internodes_model.find_interface_nodes()

            # Assemble model
            internodes_model.assemble_full_model()

            # Solve model
            displacements, lambdas = internodes_model.solve_direct()

            # Update the interface nodes and check if it converged
            converged = internodes_model.update_interface(displacements, lambdas)

            if converged:
                break

        assert i+1 < max_iter

        positions = internodes_model.nodal_positions + displacements
        positions_interface_secondary = positions[internodes_model.nodes_interface_secondary]

        a_list[j] = np.max(sp.spatial.distance.cdist(positions_interface_secondary, positions_interface_secondary))/2
        u_list[j] = np.min(positions_interface_secondary[:, 1])

    np.testing.assert_allclose(get_theoretical_normal_displacement(R, d_list, E, nu), u_list, atol=5e-3)
    """

if __name__ == '__main__':
    pytest.main()
