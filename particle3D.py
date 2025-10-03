"""
CompMod Ex2: Particle3D, a class to describe point particles in 3D space

An instance describes a particle in Euclidean 3D space: 
velocity and position are [3] arrays

Author: Elizabeth Williams
Number: s2137788

"""
import numpy as np


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

    Attributes
    ----------
    label: name of the particle
    mass: mass of the particle
    position: position of the particle
    velocity: velocity of the particle

    Methods
    -------
    __init__
    __str__
    kinetic_energy: computes the kinetic energy
    momentum: computes the linear momentum
    update_position_1st: updates the position to 1st order
    update_position_2nd: updates the position to 2nd order
    update_velocity: updates the velocity

    Static Methods
    --------------
    read_file: initializes a P3D instance from a file handle
    total_kinetic_energy: computes total K.E. of a list of particles
    com_velocity: computes centre-of-mass velocity of a list of particles
    """

    def __init__(self, label, mass, position, velocity):
        """
        Initialises a particle in 3D space.

        Parameters
        ----------
        label: str
            name of the particle
        mass: float
            mass of the particle
        position: [3] float array
            position vector
        velocity: [3] float array
            velocity vector
        """
        label = str(label)
        mass = float(mass)
        for i in range (0,3):
            position[i] = float(position[i])
        position = np.array([position[0], position[1], position[2]]) 
        for i in range (0,3):
            velocity[i] = float(velocity[i])      
        velocity = np.array([velocity[0], velocity[1], velocity[2]])           
        self.label = label
        self.mass = mass
        self.position = position
        self.velocity = velocity 
        
    def __str__(self):
        """
        Return an XYZ-format string. The format is
        label    x  y  z

        Returns
        -------
        string: str
            XYZ-format string of a Particle3D instance
        """
        string = str(self.label) + " " + str(self.position[0]) + " " + str(self.position[1]) + " " + str(self.position[2])
        return string

    def kinetic_energy(self):
        """
        Returns the kinetic energy of a Particle3D instance

        Returns
        -------
        KE: float
            ke = 1/2 m v**2 (kinetic energy of particle)
        """
        velocity_modulus_squared = (self.velocity[0])**2 + (self.velocity[1])**2 + (self.velocity[2])**2
        KE = 0.5*self.mass*velocity_modulus_squared
        return  KE

    def momentum(self):
        """
        Returns the linear momentum of a Particle3D instance

        Returns
        -------
        momentum: [3] float array
            p = mv (linear momentum of particle)
        """
        momentum = self.mass*self.velocity
        return  momentum

    def update_position_1st(self, dt):
        """
        Updates the position (to the 1st order) of a Particle3D instance
        
        Parameters
        ----------
        dt: float
            how much time has passed
        """
        self.position = self.position + dt*self.velocity
       
    def update_position_2nd(self, dt, force):
        """
        Updates the position (to the 2nd order) of a Particle3D instance
        
        Parameters
        ----------
        dt: float
            how much time has passed
        force: float
            forced used on particle
        """
        
        self.position = self.position + dt*self.velocity + 0.5*(dt**2)*force/self.mass

    def update_velocity(self, dt, force):
        """
        Updates the velocity of a Particle3D instance
        
        Parameters
        ----------
        dt: float
            how much time has passed
        force: float
            forced used on particle
        """
        force = np.array([force[0], force[1], force[2]])
        self.velocity = self.velocity + dt*force/self.mass
        

    @staticmethod
    def read_line(line):
        """
        Creates a Particle3D instance given a line of text.

        The input line should be in the format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>

        Parameters
        ----------
        filename: line
            Readable file handle in the above format

        Returns
        -------
        p: Particle3D
        """
        particle_data = line.split()
        label = str(particle_data[0])
        mass = float(particle_data[1])
        position = np.array([float(particle_data[2]), float(particle_data[3]), float(particle_data[4])])
        velocity = np.array([float(particle_data[5]), float(particle_data[6]), float(particle_data[7])])
        p = Particle3D(label, mass, position, velocity)
        return  p

    @staticmethod
    def total_kinetic_energy(particles):
        """
        Computes the total kinetic energy of a list of P3D's

        Parameters
        ----------
        particles: list
            A list of Particle3D instances

        Returns
        -------
        total_KE: float
            total kinetic energy of all particles in list
        """
        total_KE = 0
        for p in particles:
            KE = p.kinetic_energy()
            total_KE += KE
        
        return  total_KE

    @staticmethod
    def com_velocity(particles):
        """
        Computes the CoM velocity of a list of P3D's

        Parameters
        ----------
        particles: list
            A list of Particle3D instances

        Returns
        -------
        com_vel: array
            Centre-of-mass velocity of all particles in list
            com_vel = (m1v1 + m2v2 +...)/(m1 + m2 +...)
        """
        numerator = 0
        denomenator = 0
        for p in particles:
            numerator += p.velocity*p.mass
            denomenator += p.mass
        com_vel = numerator/denomenator
        return  com_vel