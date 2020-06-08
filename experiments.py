import time
from simulate import simulate, find_nearest_state_in_tree, prune_by_ideal_curve, find_nearest_state_whithin_square
from math import sqrt, radians
import matplotlib.pyplot as plt
from functions import path_root_to_node, all_transitions, transitions_reduced, transitions_preperching
from model import c, tc
from objects import FlightState

def simulate_experiment(initial_state, final_state, criteria, ts, transitions, witness_nodes = 30, k_prunning = 20, find_nearest_function=None):
    '''
        Given:
            :param initial_state: tree of flight states
            :param final_state: final state
            :param criteria: function for comparing two states
            :param ts: time steps (adimensional)
            :param transitions: function for compare two states
            :param witness_nodes: function for compare two states
            :param k_prunning: function for compare two states
            :param find_nearest_function: function for compare two states

        Return the nearest state to the final state (S_f) given by planning algorithm statrting from the initial state
        (initial_state), using the time steos in (ts), the flight maneuvers in (transitions), the heuristic parameters
        Kd (k_pruning) and Kw (witness_nodes), and the function (find_nearest_function) for selecting states
    '''

    def prune_function(node):
        return prune_by_ideal_curve(initial_state, final_state, node, k_prunning)

    def condition(Y):

        if (Y[-2] <= final_state.x) and (0 <= Y[0] <= 20) and (-10 <= Y[1] <= 10) and \
                (-10 <= Y[2] <= 10) and (radians(-60) <= Y[3] <= radians(60)):
            return True

        return False


    t = time.time()
    tree = simulate(initial_state, final_state, ts, transitions, prune_function, witness_nodes, condition)
    # nearest = find_nearest_state_in_tree(tree, final_state, criteria)
    if find_nearest_function == None:
        nearest = find_nearest_state_in_tree(tree, final_state, criteria)
    else:
        nearest = find_nearest_function(tree, final_state)

    path = path_root_to_node(tree, nearest.tag)
    t = time.time() - t
    print("cant de nodos:", len(tree))

    actions = [s.split(",") for s in nearest.tag.split("_")[1:]]
    actions = [(float(a), float(f), float(t)) for a, f, t in actions]
    cost = nearest.data.cost
    distance = criteria(nearest.data, final_state)

    return tree, path, actions, distance, cost, t

def save_img(S_0, S_f, tree, branch, path, id):
    '''
        Save image of the tree of states and the planed path between the initial and the final state
    '''

    path_list = tree.paths_to_leaves()

    for b in path_list:

        xs, zs = [], []
        for index, node_id in enumerate(b):
            data = tree.get_node(node_id).data
            xs.append(data.x*c/2)
            zs.append(data.z*c/2)

            plt.plot(xs, zs, "-k")

    Xs = []
    Zs = []
    for node in branch:
        Xs.append(node.data.x*c/2)
        Zs.append(node.data.z*c/2)

    plt.plot(Xs, Zs, "yo-")
    plt.plot([S_0.x*c/2], [S_0.z*c/2], "ro")
    plt.plot([S_f.x*c/2], [S_f.z*c/2], "ro")

    plt.xlabel("x(m)")
    plt.ylabel("z(m)")
    plt.savefig(path+str(id))
    plt.close()

## Weigths = [X, Z, Velocity,Pitch]
def compare_states(S_0, S_1, weights=(1, 1, 1, 1)):
    '''
        Weighted distance between two states using variables X, Z, speed and pitch parameters
    '''

    square_diff = [weights[0]*(S_0.x - S_1.x)*c/2, weights[1]*(S_0.z - S_1.z)*c/2,
                   weights[2]*(S_0.Ub_value() - S_1.Ub_value()), weights[3]*(S_0.theta - S_1.theta)]
    return sqrt(sum([d**2 for d, w in zip(square_diff, weights)]))


if __name__ == '__main__':

    ## Planing example

    path = "images/"

    ts = [12 / tc]              ## list of time steps in seconds (adimensional)
    criteria = compare_states   ## function for comparing states
    dw, dk = 20, 2 * 15 / c     ## Kd and Kw heuristic parameters

    transitions = transitions_reduced()
    s_0 = FlightState(0, 0, 1, 0, 0, 0, 0, 0)
    s_f = FlightState(0, 0, 1, 0, 0, 0, 2 * 200 / c, 2 * -50 / c)
    nearest = find_nearest_state_whithin_square

    try:
        result = simulate_experiment(s_0, s_f, criteria, ts, transitions, witness_nodes=dw, k_prunning=dk, find_nearest_function=nearest)
        save_img(s_0, s_f, result[0], result[1], path, "example experiment")

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
