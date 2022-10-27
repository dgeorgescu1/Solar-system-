"""
 CompMod Ex2: Particle3D, a class to describe point particles in 3D space

 An instance describes a particle in Euclidean 3D space: 
 velocity and position are [3] arrays

 Includes time integrator methods to update the position of the particle using 
 the velocity (1st order) and second order using an inputted timestep and force. 
 Another time integrator method to update the velocity using an inputted force 
 and timestep

author: Denis Georgescu
s1974479

"""
import math
import numpy as np

# Constants should go here
G0 = 6.67384E-11            # From CODATA 2010
ASTRO_U = 149597870700.0    # From IAU resolution 2012/08 B2
YEAR = float(3.15576e7)     # Julian year = 365.25 days


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the particle
    mass: mass of the particle
    pos: position of the particle
    vel: velocity of the particle

        Methods:
    __init__  - initialises the particle
    __str__  - returns a string of the particle label and postition 
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos_1st - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_particle - initializes a P3D instance from a file handle
    sys_kinetic - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    """

    def __init__(self, label, mass, pos, vel):
        """
        Initialises a particle in 3D space

        :param label: String w/ the name of the particle
        :param mass: float, mass of the particle
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """
        self.label = str(label)
        self.mass = float(mass)
        self.position = pos 
        self.velocity = vel


    @staticmethod
    def new_particle(data):
        """
        Initialises a Particle3D instance given an input file handle.
        
        The input file should contain one per planet in the following format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>
        
        :param inputFile: Readable file handle in the above format

        :return Particle3D instance
        """
        line = data.readline()
        token = line.split()
        label = str(token[0])
        mass = float(token[1])
        position = np.array([float(token[2]), float(token[3]), float(token[4])])
        velocity = np.array([float(token[5]), float(token[6]), float(token[7])])   
        
        return Particle3D(label, mass, position, velocity)
        

    def __str__(self):
        """
        XYZ-compliant string. The format is
        <label>    <x>  <y>  <z>
        """
        label = self.label
        x = self.position[0]
        y = self.position[1]
        z = self.position[2]
        xyz_string = str(label) + "    " + str(x) + " " + str(y) + " " + str(z)
        return xyz_string


    def kinetic_e(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2 m v**2
        """
        return 0.5*self.mass*(np.linalg.norm(self.velocity))**2


    def momentum(self):
        """
        Returns the linear momentum of a Particle3D instance
        :return p: (3) float np.array, m*v
        """
        return self.mass*self.velocity


    def update_pos_1st(self, dt):
        """
        1st order position update

        :param dt: timestep
        """
        self.position += dt*self.velocity


    def update_pos_2nd(self, dt, force):
        """
        2nd order position update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """
        self.position += dt*self.velocity + (dt**2)*force/(2*self.mass)


    def update_vel(self, dt, force):
        """
        Velocity update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """
        self.velocity += dt*force/self.mass


    @staticmethod
    def sys_kinetic(p3d_list):
        """
        Returns the kinetic energy of the whole system

        :param p3d_list: list in which each item is a P3D instance

        :return sys_ke: \sum 1/2 m_i v_i^2 
        """
        sys_ke = 0.0
        for p in p3d_list:
            sys_ke += p.kinetic_e()
            
        return sys_ke


    @staticmethod
    def com_velocity(p3d_list):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return total_mass: The total mass of the system 
        :return com_vel: Centre-of-mass velocity
        """
        total_mass = 0.0
        total_momentum = np.zeros(3)
        for p in p3d_list:
            total_mass += p.mass
            total_momentum += p.momentum()
        
        com_vel = total_momentum/total_mass
        return total_mass, com_vel






