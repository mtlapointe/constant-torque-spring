""" Example script on how to use the TorqueSpring object to create springs and plot them.

"""


from torque_spring import TorqueSpring
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn')


def plot_spring(ax, spring, x_shift=0.0):
    """ Plot spring given x, y and theta arrays for spring curve

    Args:
        ax: MatPlotLib axis
        spring: TorqueSpring object
        x_shift: Shift plot on x-axis (plot defaults to roller center at 0)

    Returns:
        ax

    """


    x = spring.spline_x_vals
    y = spring.spline_y_vals
    theta = spring.spline_theta
    Rd = spring.Rd
    Rn = spring.Rn
    Rr = spring.Rr

    # Get drum x, y location
    x_drum = x[-1] + np.sin(theta[-1]) * Rd
    y_drum = y[-1] - np.cos(theta[-1]) * Rd

    ang = np.linspace(0, 2 * np.pi, 100)  # Angle segments for circle plot
    ax.plot(Rd * np.cos(ang) + x_drum + x_shift, Rd * np.sin(ang) + y_drum, 'k:')  # Plot drum circle
    ax.plot((x[-1] + x_shift, x_drum + x_shift), (y[-1], y_drum), 'k-')  # Plot drum to spring tangent line
    ax.plot(Rr * np.cos(ang) + x_shift, Rr * np.sin(ang), 'k:')  # Plot roller circle
    ax.plot((x[0] + x_shift, 0 + x_shift), (y[0], 0), 'k-')  # Plot roller to spring tangent line

    ax.plot(np.array(x) + x_shift, y)  # Plot spring curve

    return ax


fig, ax = plt.subplots(figsize=(6, 6))
ax.axis('equal')

spring = TorqueSpring(3.0, 1.0, 2.25, .011, 2.0, roller_angle=110)
plot_spring(ax, spring, -spring.dist_drum)
spring = TorqueSpring(3.0, 1.0, 2.25, .011, 2.0, roller_angle=90)
plot_spring(ax, spring, -spring.dist_drum)
spring = TorqueSpring(3.0, 1.0, 2.25, .011, 2.0, roller_angle=100)
plot_spring(ax, spring, -spring.dist_drum)




