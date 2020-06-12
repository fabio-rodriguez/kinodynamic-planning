# kinodynamic-planning
Kinodynamic planning algorithm to compute ornithopters trajectories that minimize energy consumption. A nonlinear dynamic model is used together with a tree-based heuristic search for trajectory planning.

Details of the algorithm and the benchmark experiments can be seen in the publication: 

*Kinodynamic Planning for an Energy-Efficient Autonomous Ornithopter*. F. Rodriguez, J. M. Diaz-Bañez, E. Sanchez-Laulhe, J. Capitan and A. Ollero. 

Requirements:
- Python 3.X
- Numpy/Scipy
- Treelib
- Matplotlib (optional)

Usage:

Run the main program with the "-default" option to perform a random simulation of the planning algorithm.

Run the main program without any options to perform a simulation with a specific start and final state. The program will ask these and additional parameters. 
