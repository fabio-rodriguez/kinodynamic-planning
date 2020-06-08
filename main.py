import sys, os
from experiments import *
from objects import FlightState
from random import uniform

if __name__ == "__main__":

    path = "images/"

    if len(sys.argv)>1 and sys.argv[1] == "-default":

        try:
            print("***Planning random simulation example***")
            print()

            ts = [12 / tc]  ## list of time steps in seconds (adimensional)
            criteria = compare_states  ## function for comparing states
            dw, dk = 20, 2 * 15 / c  ## Kd and Kw heuristic parameters

            transitions = transitions_reduced()
            s_0 = FlightState(0, 0, 1, 0, 0, 0, 0, 0)
            s_f = FlightState(0, 0, 1, 0, 0, 0, 2 * uniform(200,250) / c, 2 * uniform(-100,20) / c)
            nearest = find_nearest_state_whithin_square ## nearest state in tree

            result = simulate_experiment(s_0, s_f, criteria, ts, transitions, witness_nodes=dw, k_prunning=dk, find_nearest_function=nearest)
            save_img(s_0, s_f, result[0], result[1], path, "main default experiment")

            print()
            print("***************************")
            print()
            print("energy (W):", result[-2])
            print("Precision (Delta):", result[-3])
            print("computation time (s):", result[-1])
            print()
            print("***************************")
            print()

        except:

            print("ERROR: Some error occur, please, check that the simulation parameters are within the bounds")

    else:

        try:

            while True:
                print("*** Run the program with (-default) option for an example simulation")
                print()
                print("Select the simulations parameters")
                print()

                initial = input("Insert initial state position X and Z separated by space (meters):").split()
                s0x, s0z = 2*float(initial[0])/c, 2*float(initial[1])/c
                print()

                final = input("Insert final state position X and Z separated by space (meters):").split()
                sfx, sfz = 2*float(final[0])/c, 2*float(final[1])/c
                print()

                timesteps = input("Insert time steps in seconds separated by space:")
                ts = [float(ti)/tc for ti in timesteps.split()]
                print()

                dw = int(input("Insert Kw parameter:"))
                print()

                dk = 2*float(input("Insert Kd parameter (meters):"))/c
                print()

                criteria = compare_states  ## function for comparing states
                transitions = transitions_reduced()

                s_0 = FlightState(0, 0, 1, 0, 0, 0, s0x, s0z)
                s_f = FlightState(0, 0, 1, 0, 0, 0, sfx, sfz)
                nearest = find_nearest_state_whithin_square  ## nearest state in tree

                result = simulate_experiment(s_0, s_f, criteria, ts, transitions, witness_nodes=dw, k_prunning=dk,
                                             find_nearest_function=nearest)
                save_img(s_0, s_f, result[0], result[1], path, "main default experiment")

                print()
                print("***************************")
                print()
                print("energy (W):", result[-2])
                print("Precision (Delta):", result[-3])
                print("computation time (s):", result[-1])
                print()
                print("***************************")
                print()

        except:

            print("ERROR: Some error occur, please, check that the simulation parameters are within the bounds")

