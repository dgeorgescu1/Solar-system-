"""
Compomod project: Solar system
A program to simulate an N body system from given initial conditions using a Velocity Verlet time integrator method
Writes a .xyz file containing positions of N particles, writes observables to files

Denis Georgescu s1974479
"""

#Importing desired modules
import numpy as np
from Particle3D import Particle3D
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import sys


#Setting constant required for calculations (CODATA)
G_SI = 6.67430E-11                  #Gravitational constant in SI units
AU = 149597870700                   #Astronomical unit in metres
days = 24*60*60                     #Number of seconds in a day
G = G_SI/((AU**3/days**2)*(1E-24))  #Gives converted gravitational constant in AU, Days, and 10^24kg
YEAR = 365.24219                    #Earth year in earth days


def gravity_calculations(m1, m2, sep, r):
    """
    

    Parameters
    ----------
    m1 : float mass of first particle.
    m2 : float mass of second particle.
    sep : 3 array containing separtions in x y z vector.
    r : magnitude of separations vector sep.

    Returns
    -------
    TYPE
        force on particle as a 3 array.
    TYPE
        potential energy float.
        
    The function is used to calculate gravitational force as a vector and potential energy as a scaler between two particles
    """
    if r == 0:  #avoids dividing by 0 for separations between the same particles
        return 0, 0
    else:
        force = -G*m1*m2*sep/r**3
        energy = -G*m1*m2/r
        return force, energy


def sep_matrix(p_list, length):
    """
    

    Parameters
    ----------
    p_list : list containing particle instances.
    length : float length of the particle list (the number of particles).

    Returns
    -------
    seps : An (n x n x 3) matrix as numpy arrays containing the pair separations between all N particles .

    The function is used to calculate the separations between particles and store them in a matrix of the required size
    """
    #set matrix to zero each time calculating
    #i != j leaves diagonal elements as zero
    seps = np.zeros((length, length, 3))
    for i in range (length):
            for j in range (length):
                if i != j:
                    seps[i, j] = p_list[i].position - p_list[j].position
    return seps


def moduli_matrix(seps, length):
    """
    

    Parameters
    ----------
    seps : the (n x n x 3) matrix calculated using sep_matrix function.
    length : TYPE
        DESCRIPTION.

    Returns
    -------
    moduli : n x n matrix where each element now is the magnitude of the separation vector.

    The function calculates scaler separations between particles and stores them in another matrix
    """
    moduli = np.zeros((length, length))
    for i in range (length):
        for j in range(length):
            moduli[i, j] = np.linalg.norm(seps[i, j])
    return moduli


def com_vel_correction(p_list, length):
    """
    

    Parameters
    ----------
    p_list : list containing particle instances.
    length : float number of particles.

    Returns
    -------
    p_list : Returns the list of particles 
    
    The function is used to update the velocities of each particle to avoid a centre of mass shift of the whole system over time
    """
    sys_mass, sys_vel = Particle3D.com_velocity(p_list)
    
    for i in range(length):
        p_list[i].velocity = p_list[i].velocity - sys_vel
    return p_list  
    

def main():
    """
    #Quits programme if number of input arguments are incorrect
    if len(sys.argv) != 4:
        print ("Wrong number of arguments")
        print ("Correct form is <particle setup file> <time step (in earth days)> <total run time (in earth years)>")
        sys.exit()
    
    #Otherwise initialises system conditions
    else:
        initialising_file = sys.argv[1]
        dt = float(sys.argv[2])
        runtime = float(sys.argv[3])*YEAR
        steps = int(runtime/dt)
    """
    
    initialising_file = "solar_inputs.txt"
    dt = 1
    runtime = 600*YEAR
    steps = int(runtime/dt)
    
    #Opens particle starting conditions file 
    file = open(initialising_file, "r")
    data = open(initialising_file, "r")
    length = len(file.readlines())          #Determines number of particles
    
    #Creates particle list and appends all instances to it
    p_list = []
    for i in range(length):
        p_list.append(Particle3D.new_particle(data))
    
    #Finds arguments for sun and earth bodies which allow the particle initialiser file to be in any order
    list_sun = [x for x in p_list if x.label == "sun"]
    index_sun = p_list.index(list_sun[0])
    list_earth = [x for x in p_list if x.label == "earth"]
    index_earth = p_list.index(list_earth[0])    
    
    #Corrects for centre of mass shift
    com_vel_correction(p_list, length)
    
    #Creates .xyz outfile simply named Trajectories
    outfile = open("Trajectories.xyz", "w")
    
    #Calculates initial separations and moduli
    seps = sep_matrix(p_list, length)
    moduli = moduli_matrix(seps, length)
     
    #Creates two arrays containing zeros to contain total forces on each particle which will be used in integrator loop
    forces = np.zeros((length, 3))
    forces_new = np.zeros((length, 3))
    
    #Creates lists for energy outfiles and plots
    kinetic_energies = []
    potential_energies = []
    Total_energies = []
    
    #Sets time and potential energy to zero and appends time to list
    time = 0.0
    potential_energy = 0.0
    time_list = [time]     
    
    #Creates all separtion files required to calculate orbital periods and apsides
    #Moon file will contain separations from earth rather than from sun
    files = []
    for i in range(length):
        if p_list[i].label == "moon":
            file = open("moon_separation.txt", "w")
            file.write("{0:f} {1:f}\n".format(time, moduli[index_earth][i]))
        else:
            filename = p_list[i].label + "_separation.txt"
            file = open(filename, "w")
            file.write("{0:f} {1:f}\n".format(time, moduli[index_sun][i]))
            
        files.append(file)      #For ease later to append files to a list
    
    #Calculates the initial total forces acting on each particle preparing for integrator loop
    #Changes ith element of forces array to ith particle force
    #Calculates initial system potential energy
    for i in range (length):
            force_gravity = 0
            for j in range (length):
                force_gravity += gravity_calculations(p_list[i].mass, p_list[j].mass, seps[i,j], moduli[i,j])[0]
                potential_energy += gravity_calculations(p_list[i].mass, p_list[j].mass, seps[i,j], moduli[i,j])[1]
    
            forces[i] = force_gravity
            
    #Appends initial system energies to lists    
    kinetic_energies.append(Particle3D.sys_kinetic(p_list))
    potential_energies.append(potential_energy/2)
    Total_energies.append(potential_energy/2 + Particle3D.sys_kinetic(p_list))
    
    
    #Verlet integrator loop
    for t in range (steps):
        
        #Loop to write .xyz file with each line containing each particle location in correct format
        outfile.write(str(length) + "\n")
        outfile.write("Point" + str(t) + "\n")
        for s in p_list:
            outfile.write(s.__str__() + "\n")
        
        #Updates positions using second order method           
        for i in range(length):
            p_list[i].update_pos_2nd(dt, forces[i])
            
        #Updates time and appends to list
        time += dt
        time_list.append(time)         
        
        #Updates particle separations and moduli
        seps = sep_matrix(p_list, length)
        moduli = moduli_matrix(seps, length)
        
        #Writes required particle moduli to outfiles to be used to calculate orbit observables
        #Ignores sun file (needed sun file to keep particle indexes consistent)
        for i in range(length):
            file = files[i]
            if i == index_sun:
                continue
            if p_list[i].label == "moon":
                file.write("{0:f} {1:f}\n".format(time, moduli[index_earth][i]))
            else:
                file.write("{0:f} {1:f}\n".format(time, moduli[index_sun][i]))
                 
        #Calculates the new forces and adds values to numpy array
        #Updates velocities using Verlet method
        #Sets old forces to new forces 
        potential_new = 0
        for i in range (length):
            force_new = 0
            for j in range (length):
                force_new += gravity_calculations(p_list[i].mass, p_list[j].mass, seps[i,j], moduli[i,j])[0]
                potential_new += gravity_calculations(p_list[i].mass, p_list[j].mass, seps[i,j], moduli[i,j])[1]
                
            forces_new[i] = force_new
            
            p_list[i].update_vel(dt, 0.5*(forces[i] + force_new))
        
            forces[i] = forces_new[i]
            
        #Appends system energies to lists
        kinetic_energies.append(Particle3D.sys_kinetic(p_list))
        potential_energies.append(potential_new/2)
        Total_energies.append(potential_new/2 + Particle3D.sys_kinetic(p_list))  
        
    #Plots three energy graphs as indicated by the titles
    plt.plot(time_list, kinetic_energies)
    plt.title("Total system kinetic energy plotted against time")
    plt.xlabel("Time (earth days)")
    plt.ylabel("System kinetic energy (10^24kg AU^2 Days^-2)")
    plt.show()
    
    plt.plot(time_list, potential_energies)
    plt.title("Total system potential energy plotted against time")
    plt.xlabel("Time (earth days)")
    plt.ylabel("System potential energy (10^24kg AU^2 Days^-2)")
    plt.show()

    plt.plot(time_list, Total_energies) 
    plt.title("Total system energy plotted against time")
    plt.xlabel("Time (earth days)")
    plt.ylabel("System energy (10^24kg AU^2 Days^-2)")
    plt.show()
    
    #Creates energy files
    kinetic_file = open("kinetic_energies.txt", "w")
    potential_file = open("potential_energies.txt", "w")
    Total_energies_file = open("Total_energies.txt", "w")
    
    #Adds energies to files 
    for i in range(len(time_list)):
        kinetic_file.write("{0:f} {1:f}\n".format(time_list[i], kinetic_energies[i]))                             
        potential_file.write("{0:f} {1:f}\n".format(time_list[i], potential_energies[i]))
        Total_energies_file.write("{0:f} {1:f}\n".format(time_list[i], Total_energies[i]))
        
    #Closes files
    kinetic_file.close()
    potential_file.close()
    Total_energies_file.close()
    
    #Also calculates maximum energy variation and prints value
    start_energy = Total_energies[0]
    max_energy = max(Total_energies)
    min_energy = min(Total_energies)
    diff_in_energy = (max_energy - min_energy)/start_energy
    print ("The variation of the energy with respect to the starting energy (E(max) - E(min))/E(start), where E is the total energy of the system: " + str(diff_in_energy))
    
    #Closes separation files
    for file in files:
        file.close()
    
    length = 12
    #Writes new file which will contain all required apsides and periods
    astro_data = open("Astrophysical_data.txt", "w")
    #Opens the individual files 
    for i in range(length):
        separations = []
        filename = p_list[i].label + "_separation.txt"
        
        for line in open(filename, "r"):
            token = line.split()
            separations.append(token[1])
        
        #Again ignores sun file
        if filename == "sun_separation.txt":
            continue
        
        #Moon apsides have slightly different names and so this elif just changes the sring for the moon
        elif filename == "moon_separation.txt":
            #Max and min used to find apogee and perigee
            apogee = max(separations) 
            perigee = min(separations)
            #Find peaks to determine when moon at closest to earth
            peaks = find_peaks(separations)[0]
            #Average the peaks to find average orbit period (for moon in days)
            average_period = np.average(np.diff(peaks))*dt
            #Adds line in astro data file
            string = p_list[i].label + " " + "Orbital period in Earth days: " + str(average_period) + " " + "Apogee in AU: " + str(apogee) + " " + "Perigee in AU: " + str(perigee)
            astro_data.write(str(string) + "\n")
                
        #similar process for neptune
        elif filename == "neptune_separation.txt":
            aphelia = max(separations)
            perihelia = min(separations)
            #Neptune has very stable orbit and so prominence was needed to exagerate peaks 
            peaks = find_peaks(separations, prominence = (float(aphelia)/100))[0]
            average_period = np.average(np.diff(peaks))*dt/YEAR
            string = p_list[i].label + " " + "Orbital period in Earth years: " + str(average_period) + " " + "Aphelia in AU: " + str(aphelia) + " " + "Perihelia in AU: " + str(perihelia)
            astro_data.write(str(string) + "\n")
          
        #Same process for remaining planets   
        else:
            aphelia = max(separations)
            perihelia = min(separations)
            peaks = find_peaks(separations)[0]
            average_period = np.average(np.diff(peaks))*dt/YEAR
            string = p_list[i].label + " " + "Orbital period in Earth years: " + str(average_period) + " " + "Aphelia in AU: " + str(aphelia) + " " + "Perihelia in AU: " + str(perihelia)
            astro_data.write(str(string) + "\n")
            
            #Periods will not be calculated if objects havent made full orbits and in the outfile will be listed as nan    
            
    #Closes file
    astro_data.close()

#Runs function if condition met          
if __name__ == "__main__":
    main()