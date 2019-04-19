# Ising-Model
Simulating the 2D Ising Model for a ferromagnetic material


#README

Numerical Monte Carlo simulation of the 2D Ising Model based on the Metropolis algorithm – Lattice points are spins of ferromagnetic material. See https://en.wikipedia.org/wiki/Ising_model#Metropolis_algorithm for more information.

Code:

1) Ising_Class.py is the OOPclass which holds all the methods used in this simulation. Instance of a class is a square lattice of spins to simulate ferromagnetic material.

2) Graphs.py is the code that takes measurements of the system and generates all graphs. User is required to input in the command line: Graphs.py <temperature> <dimension of square lattice> <glauber/kawasaki (choice of dynamical algorithm)> 

3) Visualization.py is the code that visualises the system with a matplotlib animation. User is required to input in the command line: Visualization.py <temperature> <dimension of square lattice> <glauber/kawasaki>
Recommended: temp=1, dimension=50 or 100, glauber

—————————
2 Data files and 6 graphs are generated for glauber and kawasaki dynamics run for 3000 iterations for 30 increments of temperature.

Data files contain temperature, avg Energy, avg Heat Capacity, avg Magnetisation, avg Susceptibility all separated by a single space.
