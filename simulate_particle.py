"""
Symplectic Euler and Velocity Verlet time integration of
two particles moving interacting via the morse potential.

Produces plots of the separations of the particles
and their total energy, both as function of time. Also
saves both to file.

Author: Elizabeth Williams
Number: S2137788

"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from particle3D import Particle3D
from scipy.signal import find_peaks

def force_morse(p1, p2, alpha, D_e, r_e):
    """
    Return the force between two particles interacting via the morse potential.

    The force is given by
        F(r1, r2) = 2*alpha*D_e*(1-exp(-alpha(r12-r_e)))*exp(-alpha(r12-r_e)
                    in direction of r12 unit vector where r12 = r2 - r1

    Parameters
    ----------
    p1: Particle3D
    p2: particle3D
    alpha: float
    D_e: float
    r_e: float

    Returns
    -------
    force: array
    """
    
    r12_vector = np.subtract(p2.position, p1.position) # vector between the two particle positions
    r12_norm = np.linalg.norm(r12_vector) # distance between the particles (separation)
    r12_unit = r12_vector/r12_norm # unit vector of r12_vector
    exponential = np.exp(-alpha * (r12_norm - r_e))
    force = 2 * alpha * D_e * (1-exponential) * exponential * r12_unit
    return force


def potential_morse(p1, p2, alpha, D_e, r_e):
    """
    Return the potential energy of two particles interacting via the morse potential.
    
    The potential is given by
        D_e(((1-exp(-alpha(r12-r_e)))**2)-1)

    Parameters
    ----------
    p1: Particle3D
    p2: particle3D
    alpha: float
    D_e: float
    r_e: float

    Returns
    -------
    potential: array
    """
    r12_vector = np.subtract(p2.position, p1.position) # vector between the two particle positions
    r12_norm = np.linalg.norm(r12_vector) # distance between the particles (separation)
    exponential = np.exp(-alpha * (r12_norm - r_e))
    potential = D_e * (((1-exponential)**2) - 1)
    return potential

def wavenumber(separations, dt):
    """
    Return the wavenumber of two particles interacting via the morse potential.
    
    The wavenumber is given by
        wavenumber = 1/period*c (where frequency = 1/period)

    Parameters
    ----------
    separations: array
    dt: float

    Returns
    -------
    wavenumber: float
    """
    c = 299792458 # speed of light
    peaks = find_peaks(separations)[0] # creating an array of the peak positions
    n = len(peaks) # number of peaks
    period = (peaks[n-1] - peaks[0])/(n-1) # average distance between peaks gives the average period
    period_secs = period*1.018050571*(10**(-14))*dt # change units to seconds
    wavenumber = 0.01/(period_secs*c) # find wavenumber and convert to cm^-1
    return wavenumber


def main():
    if len(sys.argv) != 4 : # Checking that the length of the user input is correct
    # Printing out a helpful error message if the user is missing one of four things:
    # the name of the file, the name of the particle file, euler or verlet, the name of the output file
        print("You left out the name of the output file when running.")
        print("In spyder, run like this instead:")
        print(f"    %run {sys.argv[0]} {sys.argv[1]} <euler or verlet> <desired output file>")
        sys.exit(1)
    else:
        inputfile_name = sys.argv[1] # name of the particle data file 
        mode = sys.argv[2] # time integration mode
        outfile_name = sys.argv[3] # name of the output file

    # Open the output file for writing
    outfile = open(outfile_name, "w") 
    # Open the input file for reading
    inputfile = open(inputfile_name, "r") 
    
    # Read the simulation parameters from the input file
    lines = inputfile.readlines() 
    
    dt = float(lines[1].split()[0]) # time step parameter
    numstep = int(lines[1].split()[1]) # the number of steps the simulation will run through
    time = 0.0
    
    # morse potential parameters dependent on the particle type
    D_e = float(lines[3].split()[0]) 
    r_e = float(lines[3].split()[1])
    alpha = float(lines[3].split()[2])
    
    # use Particle3D to create 2 particles with the specified initial conditions
    p1 = Particle3D.read_line(lines[4]) 
    p2 = Particle3D.read_line(lines[5])

    # Get initial force
    force = force_morse(p1, p2, alpha, D_e, r_e)

    # Write out starting time, separation, and total energy values to the output file.
    inital_energy = p1.kinetic_energy() + p2.kinetic_energy() + potential_morse(p1, p2, alpha, D_e, r_e)
    separation = np.linalg.norm(p2.position - p1.position)
    outfile.write(f"{time}    {separation}    {inital_energy}\n")
    
    # Initialise numpy arrays used to record the trajectories of the particles.
    times = np.zeros(numstep)
    separations = np.zeros(numstep)
    energies = np.zeros(numstep)
    positionsx = np.zeros(numstep)
    positionsy = np.zeros(numstep)
    positionsx2 = np.zeros(numstep)
    positionsy2 = np.zeros(numstep)

    # Start the time integration loop
    for i in range(numstep):

        # Update the positions and velocities
        # This will depend on the time integration scheme used
        if mode == "euler":
            # Update particle positions
            p1.update_position_1st(dt)
            p2.update_position_1st(dt)
            
            # Calculate force
            force = force_morse(p1, p2, alpha, D_e, r_e)

            # Update particle velocities
            p1.update_velocity(dt, force)
            p2.update_velocity(dt, -force)

        elif mode == "verlet":
            # Update particle positions using previous force
            p1.update_position_2nd(dt, force)
            p2.update_position_2nd(dt, -force)
            
            # Get the force value for the new positions
            force_new = force_morse(p1, p2, alpha, D_e, r_e)

            # Update particle velocity by averaging current and new forces
            p1.update_velocity(dt, 0.5*(force+force_new))
            p2.update_velocity(dt, -0.5*(force+force_new))
            
            # Re-define force value for the next iteration
            force = force_new
        else:
            raise ValueError(f"Unknown mode {mode} - should be euler or verlet")

        # Increase time
        time += dt
        
        # Output particle information
        energy = p1.kinetic_energy() + p2.kinetic_energy() + potential_morse(p1, p2, alpha, D_e, r_e)
        separation = np.linalg.norm(p2.position - p1.position)
        outfile.write(f"{time} {separation} {energy}\n")

        # Store time, separation and energy in numpy arrays to plot later
        times[i] = time
        separations[i] = separation
        energies[i] = energy
        positionsx[i] = p1.position[0]
        positionsy[i] = p1.position[1]
        positionsx2[i] = p2.position[0]
        positionsy2[i] = p2.position[1]

    # Close the output file
    outfile.close()

    
    if mode == "euler":
        name = "Symplectic Euler"
    elif mode == "verlet":
        name = "Velocity Verlet"
        
    # Plot particle separations
    pyplot.figure()
    pyplot.title(name + ": separation vs time")
    pyplot.xlabel("Time [T]")
    pyplot.ylabel("Separation (angstroms)")
    pyplot.plot(times, separations)
    pyplot.show()
 
    # Plot total energy of particles
    pyplot.figure()
    pyplot.title(name + ": total energy vs time")
    pyplot.xlabel("Time [T]")
    pyplot.ylabel("Energy (eV)")
    pyplot.plot(times, energies)
    pyplot.show()

    # Plot x and y positions
    
    pyplot.figure()
    pyplot.title(name + ": x position vs time of " + p1.label + " spin particle")
    pyplot.xlabel("Time [T]")
    pyplot.ylabel("Separation (angstroms)")
    pyplot.plot(times, positionsx)
    pyplot.plot(times, positionsx2)
    pyplot.show()
    
    pyplot.figure()
    pyplot.title(name + ": y position vs time of " + p1.label + " spin particle")
    pyplot.xlabel("Time [T]")
    pyplot.ylabel("Separation (angstroms)")
    pyplot.plot(times, positionsy)
    pyplot.plot(times, positionsy2)
    pyplot.show()


    # calls wavenumber function to find and print the wavenumber value of the simulation
    wavenumber_value = wavenumber(separations, dt)
    print("Wave number: " + str(wavenumber_value) + " cm^-1")
    
    #find and print energy inaccuracy
    energy_inaccuracy = (max(energies) - min(energies))/inital_energy
    print("Energy inaccuracy: " + str(energy_inaccuracy))

# the main function will only be run if we run this file
# i.e. not if we just import it from another python file
if __name__ == "__main__":
    main()

