import pytest
import numpy as np
import scipy as sp

from helper import get_theoretical_normal_displacement
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
    n_supports_max = np.ones_like(n_supports) / wendland(c_max)
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

        rbf_radius_parameters_primary = compute_rbf_radius_parameters(positions_primary, positions_secondary)
        _check_rbf_radius_conditions(positions_primary, rbf_radius_parameters_primary)

        rbf_radius_parameters_secondary = compute_rbf_radius_parameters(positions_secondary, positions_primary)
        _check_rbf_radius_conditions(positions_secondary, rbf_radius_parameters_secondary)

def test_find_interface_nodes_regular():
    
    for dim in [2, 3]:

        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
        positions_primary = _get_uniform_interface_grid(dim, g_zeros, n_grid=10)
        positions_secondary = _get_uniform_interface_grid(dim, g_parabolic, n_grid=10)
        nodes_primary = np.arange(len(positions_primary))
        nodes_secondary = np.arange(len(positions_secondary))

        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary = find_interface_nodes(positions_primary, positions_secondary, nodes_primary, nodes_secondary, rbf=wendland, C=0.95)
        _check_node_contained_in_support(positions_primary[nodes_primary], positions_secondary[nodes_secondary], rbf_radius_parameters_primary)
        _check_node_contained_in_support(positions_secondary[nodes_secondary], positions_primary[nodes_primary], rbf_radius_parameters_secondary)

def test_compute_rbf_radius_parameters_random():

    for dim in [2, 3]:
        
        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
        positions_primary = _get_random_interface_grid(dim, g_zeros, n_grid=10)
        positions_secondary = _get_random_interface_grid(dim, g_parabolic, n_grid=10)

        rbf_radius_parameters_primary = compute_rbf_radius_parameters(positions_primary, positions_secondary)
        _check_rbf_radius_conditions(positions_primary, rbf_radius_parameters_primary)

        rbf_radius_parameters_secondary = compute_rbf_radius_parameters(positions_secondary, positions_primary)
        _check_rbf_radius_conditions(positions_secondary, rbf_radius_parameters_secondary)

def test_find_interface_nodes_random():
    
    for dim in [2, 3]:

        g_zeros = lambda x: np.zeros(x.shape[1])
        g_parabolic = lambda x: np.sum(x**2, axis=0) + 0.05
        positions_primary = _get_random_interface_grid(dim, g_zeros, n_grid=10)
        positions_secondary = _get_random_interface_grid(dim, g_parabolic, n_grid=10)
        nodes_primary = np.arange(len(positions_primary))
        nodes_secondary = np.arange(len(positions_secondary))

        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary = find_interface_nodes(positions_primary, positions_secondary, nodes_primary, nodes_secondary, rbf=wendland, C=0.95)
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
        nodes_primary = np.arange(len(positions_primary))
        nodes_secondary = np.arange(len(positions_secondary))

        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary = find_interface_nodes(positions_primary, positions_secondary, nodes_primary, nodes_secondary, rbf=wendland, C=0.95)
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

        rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary = find_interface_nodes(positions_primary, positions_secondary, nodes_primary, nodes_secondary, rbf=wendland, C=0.95)
        R12 = construct_interpolation_matrix(positions_primary[nodes_primary], positions_secondary[nodes_secondary], rbf_radius_parameters_secondary, rbf=wendland_rbf)
        R21 = construct_interpolation_matrix(positions_secondary[nodes_secondary], positions_primary[nodes_primary], rbf_radius_parameters_primary, rbf=wendland_rbf)
        positions_interpolated_primary = R12 * positions_secondary[nodes_secondary]
        positions_interpolated_secondary = R21 * positions_primary[nodes_primary]

        np.testing.assert_allclose(g_parabolic(positions_interpolated_primary[:, :-1].T), positions_interpolated_primary[:, -1], atol=2e-2)
        np.testing.assert_allclose(g_zeros(positions_interpolated_secondary[:, :-1].T), positions_interpolated_secondary[:, -1], atol=2e-2)

def test_contact3d_problem():

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
    model.applyBC(aka.FixedValue(-0.1, aka._y), 'secondary_fixed')
    model.applyBC(aka.FixedValue(0., aka._z), 'secondary_fixed')

    # Set initial conditions
    internodes_model = ContactMechanicsInternodes(spatial_dimension, mesh, model, 'primary_candidates', 'secondary_candidates')

    # Run internodes algorithm
    d0 = 0.05
    R = 0.5

    E = model.getMaterial(0).getReal("E")
    nu = model.getMaterial(0).getReal("nu")

    d_list = np.linspace(0.05, 0.1, 2)
    a_list = np.empty_like(d_list)
    u_list = np.empty_like(d_list)

    for j, d in enumerate(d_list):

        model.applyBC(aka.FixedValue(-d+d0, aka._y), 'secondary_fixed')

        internodes_model = ContactMechanicsInternodes(spatial_dimension, mesh, model, 'primary_candidates', 'secondary_candidates')

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

        positions = internodes_model.mesh.getNodes() + displacements
        positions_interface_secondary = positions[internodes_model.nodes_interface_secondary]

        a_list[j] = np.max(sp.spatial.distance.cdist(positions_interface_secondary, positions_interface_secondary))/2
        u_list[j] = np.min(positions_interface_secondary[:, 1])

    np.testing.assert_allclose(get_theoretical_normal_displacement(R, d_list, E, nu), u_list, atol=1e-2)

def test_contact2d_problem():
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
    model.applyBC(aka.FixedValue(0., aka._x), 'secondary_fixed')
    model.applyBC(aka.FixedValue(-0.1, aka._y), 'secondary_fixed')

    # Set initial conditions
    nodes_top = mesh.getElementGroup('secondary_fixed').getNodeGroup().getNodes().ravel()
    displacements = np.zeros((mesh.getNbNodes(), spatial_dimension))
    displacements[nodes_top, 1] = -0.1
    displacements = displacements.ravel()

    internodes_model = ContactMechanicsInternodes(spatial_dimension, mesh, model, 'primary_candidates', 'secondary_candidates')

    f_free = model.getExternalForce().ravel()
    f_free = f_free[internodes_model.dofs_free]
    
    # Run internodes algorithm
    for _ in range(10):
        # Find the interface nodes
        internodes_model.find_interface_nodes()

        # Assemble model
        internodes_model.assemble_interpolation_matrices()
        internodes_model.assemble_stiffness_matrix()
        internodes_model.assemble_interface_mass_matrices()
        internodes_model.assemble_B_matrices()
        internodes_model.assemble_internodes_matrix()
        internodes_model.assemble_force_term(f_free, displacements)

        # Solve model
        positions_new, displacements, lambdas = internodes_model.solve_direct(displacements)

        # Update the interface nodes and check if it converged
        converged = internodes_model.update_interface(positions_new, lambdas)

        if converged:
            break

    positions_interface_secondary = positions_new[internodes_model.nodes_interface_secondary]

    # Model parameters
    R = 0.5
    d = 0.15
    nu = model.getMaterial(0).getReal("nu")
    E = model.getMaterial(0).getReal("E")

    E_red = 0.5 * E/(1-nu**2)
    F = 4/3 * E_red * R**0.5 * d**1.5

    # Test contact circle radius
    x_min = np.min(positions_interface_secondary[:, 0])
    x_max = np.max(positions_interface_secondary[:, 0])
    a_eff_x = (x_max - x_min) / 2
    a_exp = (d*R)**0.5
    np.testing.assert_allclose(a_eff_x, a_exp, rtol=2e-1)

    # Test normal displacement
    p0 = 3*F / (2*np.pi*a_exp**2)
    u_exp = - 1/E * np.pi / 2 * p0 * a_exp
    u_eff = np.min(positions_interface_secondary[:, 1])
    np.testing.assert_allclose(u_eff, u_exp, rtol=2e-1)
    """

if __name__ == '__main__':
    pytest.main()
