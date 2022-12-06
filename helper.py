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
    return - (1-nu**2)/E * np.pi/2 * p0 * a

def plot_mesh(positions, triangle_indices, nodes_interface=None, nodes_interpenetrating=None, nodes_tension=None):
    plt.figure()
    plt.triplot(positions[:, 0], positions[:, 1], triangle_indices)
    if nodes_interface is not None:
        plt.scatter(positions[nodes_interface, 0], positions[nodes_interface, 1], color="blue")
    if nodes_interpenetrating is not None and nodes_tension is not None:
        plt.scatter(positions[nodes_interpenetrating, 0], positions[nodes_interpenetrating, 1], color="red")
        plt.scatter(positions[nodes_tension, 0], positions[nodes_tension, 1], color="green")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('scaled')
    plt.show()