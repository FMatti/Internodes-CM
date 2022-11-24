import pytest
import numpy as np

from contact_mechanics_internodes import *

def test_construct_interpolation_matrix():

    # Test if constant function is interpolated exactly
    positions = np.array([1, 2, 3])
    positions_ref = np.array([1.5, 2.5])
    radiuses_ref = np.array([1, 1])
    rbf = wendland_rbf
    R = construct_interpolation_matrix(positions, positions_ref, radiuses_ref, rbf)
    f = lambda x : np.ones_like(x)
    np.testing.assert_allclose(R * f(positions_ref), f(positions), rtol=1e-10)

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

if __name__ == '__main__':
    pytest.main()
