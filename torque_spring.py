"""Constant torque spring shape generator

This module will create a TorqueSpring object given basic design information that calculates the shape of the
transition zone between the Roller and Output Spool. Uses methods presented in NASA paper
"Advances in Analysis and Design of Constant-Torque Springs by McGuire and Yura.

Example Use:
    spring = TorqueSpring(dia_spool=3.0, dia_roller=1.0, dia_coil=2.25, thk=.011, wd=2.5, roller_angle=55)

Todo:
    * Add more design guideline warnings per the NASA paper

"""

import numpy as np
from numpy.linalg import norm
from scipy.optimize import least_squares


class TorqueSpring:
    """ TorqueSpring class contains spring design and shape information

    Attributes:
        force_roller: Normal force on the roller applied by the spring in lbs
        spline_x_vals: List of x points on spring curve shape
        spline_y_vals: List of y points on spring curve shape
        spline_theta: List of slope for points on curve
        spline_len: Length of curve
        dist_drum: Drum center from roller
        theta_p: Spring output spool tangent angle
        theta_t: Spring roller tangent angle

    """




    def __init__(self, dia_spool, dia_roller, dia_coil, thk, wd, modulus=28e6, roller_angle=110):
        """ Constructor for TorqueSpring object. Generates curve geometry based on design inputs.

        Args:
            dia_spool: Output spool diameter
            dia_roller: Roller diameter
            dia_coil: Spring natural coil diameter
            thk: Spring thickness
            wd: Spring flat width
            modulus (optional): Material modulus of elasticity in PSI. Defaults to 28e6 psi.
            roller_angle: Angle between spring normal force on Roller and line between Roller and Spool centers
                Defaults to 110, which is the optimal roller angle per the NASA paper.
        """

        self.thk = thk
        self.wd = wd
        self.E = modulus
        self.Rn = dia_coil / 2
        self.Rd = dia_spool / 2
        self.Rr = dia_roller / 2
        self.Ix = wd * (thk ** 3) / 12
        self.theta_t = roller_angle

        # Initialize spring shape variables then calculate
        self.force_roller = 0
        self.spline_x_vals, self.spline_y_vals, self.spline_theta, \
            self.spline_len, self.dist_drum, self.theta_p, self.theta_t = self.get_spring_shape(roller_angle)

    def curve_slope(self, force_roller, x):
        # Given a roller force, determine slope of spring curve for given x location
        # Raises InvalidSpringGeometry if roller force is too low (drum distance too far apart)
        with np.errstate(all='raise'):
            try:
                return np.arcsin((1 / self.Rn) * x - force_roller * (x ** 2) / (2 * self.E * self.Ix))
            except FloatingPointError:
                raise InvalidSpringGeometry

    def calc_spring_curve(self, force_roller, prec=.001):
        # Generate spring curve shape using derivative method
        # Returns dict:
        #   x: list of x points on curve
        #   y: list of corresponding y points on curve
        #   theta: list of slope for points on curve
        #   len: length of curve
        #   drum_dist: drum center from roller
        #   theta_p: spring drum tangent angle
        #   theta_t: spring roller tangent angle
        # Raises InvalidSpringGeometry if roller force is too low (drum distance too far apart)

        # Uses derivative method to iterate over spring curve for a given roller force
        s = x = theta = 0
        y = -self.Rr
        phi_init = 1 / self.Rn
        spring_seg = [[s, theta, x, y, phi_init]]
        ds = prec * self.Rd  # Precision for derivative steps - too fine will run slow so adjust for spool size

        phi = phi_init
        while phi > (-1 / self.Rd):
            dx = ds * np.cos(theta)
            dy = ds * np.sin(theta)
            s += ds
            x += dx
            y += dy
            theta += phi * ds
            phi = phi_init - force_roller * x / (self.E * self.Ix)
            self.curve_slope(force_roller, x + dx)  # Check if slope is invalid
            spring_seg.append([s, theta, x, y, phi])

        # Unpack results
        s, theta, x, y, phi = list(zip(*spring_seg))
        x_drum = x[-1] + np.sin(theta[-1]) * self.Rd
        y_drum = y[-1] - np.cos(theta[-1]) * self.Rd
        dist_drum = np.sqrt(x_drum ** 2 + y_drum ** 2)

        # Rotate xy coordinates so that the roller and drum centers are on the x-axis
        x_rot = []
        y_rot = []
        rot_angle = -np.arctan(y_drum / x_drum)
        rot_matrix = np.array(((np.cos(rot_angle), -np.sin(rot_angle)),
                               (np.sin(rot_angle), np.cos(rot_angle))))
        for xy in zip(x, y):
            xx, yy = rot_matrix.dot(xy)
            x_rot.append(xx)
            y_rot.append(yy)

        theta = list(np.array(theta) + rot_angle)  # Rotate all the spring slope values

        # Determine spring drum tangent rotation
        v1 = (x[-1] - x_drum, y[-1] - y_drum)  # Vector for drum center to spring tangent
        v2 = (-x_drum, -y_drum)  # Vector for drum to roller
        c = np.dot(v1, v2) / norm(v1) / norm(v2)  # Cosine of the angle
        theta_p = np.arccos(np.clip(c, -1, 1))  # Drum angle

        # Determine spring roller tangent rotation
        v1 = (-x[0], -y[0])  # Vector for roller center to spring tangent
        v2 = (x_drum, y_drum)  # Vector for drum to roller
        c = np.dot(v1, v2) / norm(v1) / norm(v2)  # Cosine of the angle
        theta_t = np.arccos(np.clip(c, -1, 1))  # Roller angle

        spring_len = s[-1]

        return x_rot, y_rot, theta, spring_len, dist_drum, theta_p, theta_t

    def get_roller_angle(self, force_roller):
        # Returns roller angle given a roller force. If force is too low, return an angle of 0.
        try:
            return self.calc_spring_curve(force_roller)[6]
        except InvalidSpringGeometry:
            return 0  # Make roller angle get worse with lower force

    def get_spring_shape(self, roller_angle):

        self.force_roller = \
            least_squares(lambda x: self.get_roller_angle(x[0]) - np.radians(roller_angle), np.array(20),
                          bounds=[0, 'inf']).get('x')[0]

        return self.calc_spring_curve(self.force_roller)


class InvalidSpringGeometry(Exception):
    # Exception for when spring roller force is too low
    pass


