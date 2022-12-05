import pytest
import numpy as np

from contact_mechanics_internodes import *

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

def test_find_interface_nodes():

    # Randomly generate points on a paraboloid
    n_points = 100
    random_positions = np.random.rand(n_points, 2) - 0.5
    random_positions = random_positions[np.sum(random_positions**2, axis=1) < 0.5**2]
    random_positions = np.c_[random_positions, np.sum(random_positions**2, axis=1)]

    # Find the interface nodes
    nodes_primary = np.arange(0, len(random_positions) // 2)
    nodes_secondary = np.arange(len(random_positions) // 2, len(random_positions))
    positions_primary = random_positions[nodes_primary]
    positions_secondary = random_positions[nodes_secondary]
    rbf = wendland
    rbf_radius_parameters_primary, rbf_radius_parameters_secondary, nodes_primary, nodes_secondary = find_interface_nodes(positions_primary, positions_secondary, nodes_primary, nodes_secondary, rbf)

    # Compute the distance matrices for the checks
    positions_primary = random_positions[nodes_primary]
    positions_secondary = random_positions[nodes_secondary]
    distance_matrix_MM = sp.spatial.distance.cdist(positions_primary, positions_primary)
    distance_matrix_NM = sp.spatial.distance.cdist(positions_secondary, positions_primary)
    distance_matrix_MN = sp.spatial.distance.cdist(positions_primary, positions_secondary)
    distance_matrix_NN = sp.spatial.distance.cdist(positions_secondary, positions_secondary)

    # Assert condition: [1], Page 51, Section 2.3, Equation 2
    np.fill_diagonal(distance_matrix_MM, np.inf)
    c_MM_max = np.min(distance_matrix_MM / rbf_radius_parameters_primary, axis=0)
    np.testing.assert_array_less(np.zeros_like(c_MM_max), c_MM_max)
    np.testing.assert_array_less(c_MM_max, np.ones_like(c_MM_max))
    np.fill_diagonal(distance_matrix_NN, np.inf)
    c_NN_max = np.min(distance_matrix_NN / rbf_radius_parameters_secondary, axis=0)
    np.testing.assert_array_less(np.zeros_like(c_NN_max), c_NN_max)
    np.testing.assert_array_less(c_NN_max, np.ones_like(c_NN_max))

    # Assert condition: [1], Page 51, Section 2.3, Equation 3
    C_NM_min = np.min(distance_matrix_NM / rbf_radius_parameters_primary, axis=0)
    #np.testing.assert_array_less(c_MM_max, C_NM_min)
    np.testing.assert_array_less(C_NM_min, np.ones_like(C_NM_min))
    C_MN_min = np.min(distance_matrix_MN / rbf_radius_parameters_secondary, axis=0)
    #np.testing.assert_array_less(c_NN_max, C_MN_min)
    np.testing.assert_array_less(C_MN_min, np.ones_like(C_MN_min)) 

    # Assert condition: [1], Page 51, Section 2.3, Equation 4
    n_supports_interpolation_MM = np.sum(distance_matrix_MM < rbf_radius_parameters_primary, axis=0)
    np.testing.assert_array_less(n_supports_interpolation_MM, 1 / wendland(c_MM_max))
    n_supports_interpolation_NN = np.sum(distance_matrix_NN < rbf_radius_parameters_secondary, axis=0)
    np.testing.assert_array_less(n_supports_interpolation_NN, 1 / wendland(c_NN_max))

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
    z_min = np.min(positions_interface_secondary[:, 2])
    z_max = np.max(positions_interface_secondary[:, 2])
    a_eff_x = (x_max - x_min) / 2
    a_eff_z = (z_max - z_min) / 2
    a_exp = (d*R)**0.5
    np.testing.assert_allclose(a_eff_x, a_exp, rtol=2e-1)
    np.testing.assert_allclose(a_eff_z, a_exp, rtol=2e-1)

    # Test normal displacement
    p0 = 3*F / (2*np.pi*a_exp**2)
    u_exp = - 1/E * np.pi / 2 * p0 * a_exp
    u_eff = np.min(positions_interface_secondary[:, 1])
    np.testing.assert_allclose(u_eff, u_exp, rtol=2e-1)

def test_contact2d_problem():
    
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

if __name__ == '__main__':
    pytest.main()
