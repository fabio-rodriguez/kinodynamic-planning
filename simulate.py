from math import pi, cos, sin, sqrt
from objects import FlightState
from treelib import Tree
from model import ode_dydt, k_aero, tc, c_e, c
from scipy.integrate import odeint
import warnings
from functions import euclidean_dist
warnings.filterwarnings("error")

def simulate(S_0, S_f, ts, transitions, prune_function, witness_nodes, condition=None, h=0.03):
    '''
        Given:
            :param S_0: initial state
            :param S_f: final state
            :param ts: time step list (adimensional)
            :param transitions: dictionary of possible maneuvers to choose
            :param prune_function: function for pruning states out of the corridor (Kd)
            :param witness_nodes: number of witness nodes (Kw)
            :param condition: function for pruning states out of range
            :param h: integration step

        Return the tree of states for the planning algorithm between the initial state (S_0) and the target
        state (S_f) using the flight mmaneuvers in (transitions), the time steps in (ts), and the other pruning
        paramenters
    '''

    steps_tree = Tree()
    steps_tree.create_node("root", identifier="root", data=S_0)
    return get_tree(steps_tree, ["root"], transitions, ts, prune_function, witness_nodes, condition, h, S_f)

def get_tree(steps_tree, parents_id, transitions, ts, prune_function, witness_nodes, condition, h, S_f):
    '''
        Auxiliar recursive function for "simulate", return the tree of states for the planning algorithm between
        the initial state in the root of the tree (steps_tree) and the target state (S_f)
    '''

    if not parents_id:
        return steps_tree

    new_level = []
    for t in ts:
        time_nodes = []
        tSpan = [h * i for i in range(int(t/h))]
        for id in parents_id:
            node = steps_tree.get_node(id).data
            y0, cost0 = [node.u, node.v, node.omega, node.theta, node.x, node.z], node.cost

            for angle, fa in transitions[(node.tail_angle, node.fa)]:

                model = ode_dydt(angle, fa)
                try:
                    Y = odeint(model, y0[:], tSpan, args=(0,))
                except Exception as e:
                    continue

                if euclidean_dist(Y[-1][-2:], Y[0][-2:])*c/2 < 1*t*tc:
                    continue

                aux = False
                for i in range(int((len(Y)-1)/10)-1):
                    if Y[i*10][-2] >= Y[10*(i+1)][-2]:
                        aux = True
                        break

                if aux: continue

                tim = 1
                if not condition(Y[-1]):
                    index = index_in_range(Y, condition)
                    tim = index/len(Y)
                    Y = Y[:index]

                P_t, P_w = c_e, k_aero*(fa/tc)**3
                power = (P_t + P_w)*(tim*t*tc)
                u, v, omega, theta, x, z = Y[-1]

                new_node = FlightState(angle, fa, u, v, theta, omega, x, z)
                new_node.increment_cost(cost0, power)
                if not prune_function(new_node):
                    id_node = id + "_" + str(angle) + "," + str(fa) + "," + str(t)
                    time_nodes.append((new_node, id_node, id))

        final_nodes, cluster_nodes = [], []
        for node, x, y in time_nodes:
            if euclidean_dist([node.x], [S_f.x]) < 2*(S_f.x/200)/c:
                final_nodes.append((node,x,y))
            else:
                cluster_nodes.append((node,x,y))

        if witness_nodes and len(cluster_nodes) > witness_nodes/len(ts):
            cluster_nodes.sort(key=lambda x: x[0].z)

            range_height = (cluster_nodes[-1][0].z - cluster_nodes[0][0].z)/ witness_nodes
            dict_nodes = {num:[] for num in range(witness_nodes+1)}
            for node in cluster_nodes:
                height = node[0].z - cluster_nodes[0][0].z
                try:
                    dict_nodes[int(height/range_height)].append(node)
                except:
                    continue

            cluster_nodes = [min(node_list, key=lambda x: x[0].cost) for node_list in dict_nodes.values() if node_list]

        for n, child_id, id in final_nodes:
            steps_tree.create_node(child_id, child_id, parent=id, data=n)

        for n, child_id, id in cluster_nodes:
            new_level.append(child_id)
            steps_tree.create_node(child_id, child_id, parent=id, data=n)

    return get_tree(steps_tree, new_level, transitions, ts, prune_function, witness_nodes, condition, h, S_f)

def find_nearest_state_in_tree(tree, S_f, criteria):
    '''
        Given:
            :param tree: tree of flight states
            :param S_f: final state
            :param criteria: function for compare two states

        Return the tree nearest state to the final state (S_f) in the tree of flight states (tree) using the function
        (criteria) for comparing two states
    '''

    feasible_nodes = list(tree.all_nodes())
    feasible_nodes.sort(key=lambda x: criteria(x.data, S_f))

    return feasible_nodes[0]

def find_nearest_state_whithin_square(tree, S_f, error = 2*3/c, select_function = None):
    '''
        Given:
            :param tree: tree of flight states
            :param S_f: final state
            :param error: measure of distance
            :param select_function: function for selecting the best state within a list of states

        Return the best staes in the tree (tree) with maximum distance (error) from the final state (S_f)
    '''

    feasible_nodes = list(tree.all_nodes())
    result = []
    for node in feasible_nodes:
        if abs(node.data.x - S_f.x) < 2*error and abs(node.data.z - S_f.z) <= error:
            result.append(node)

    if result:
        if select_function == None:
            result.sort(key=lambda x: x.data.cost)
        else:
            result.sort(key=lambda x: select_function(S_f,x.data))

        return result[0]
    else:
        print("***NOT FOUNDED FEASIBLE NODE***")
        print()
        print("Nearest positional node result:")

        feasible_nodes.sort(key=lambda x: euclidean_dist([x.data.x, x.data.z], [S_f.x, S_f.z]))
        return feasible_nodes[0]

def prune_by_ideal_curve(state_0, state_f, new_state, k_prunning):
    '''
        Given:
            :param state_0: initial state
            :param state_f: final state
            :param new_state: flight state
            :param k_prunning: pruning parameter (Kd)

        Return True if the flight state (new state) is not within the scaled cosine corridor between the initial state
        (state_0) and the final state (state_f), or False in other case. The parameter (k_pruning) define the width of
        the corridor
    '''


    x0, h0 = state_0.x, state_0.z
    xf, hf = state_f.x, state_f.z
    x_new, h_new = new_state.x, new_state.z

    if x_new > xf or x_new < x0:
        return True

    if not k_prunning:
        return False

    xd, hd = abs(xf - x0), abs(hf - h0)
    if hf == h0:
        return abs(h_new - h0) > k_prunning
    elif hf > h0:
        f = lambda x, k: hd/2 + hd/2 * cos(pi + x * pi / xd) + k
    else:
        f = lambda x, k: hd/2 * sin(pi / 2 + x * pi / xd) - hd/2 + k

    h1, h2 = f(x_new, k_prunning), f(x_new, -k_prunning)
    h1, h2 = min(h1, h2), max(h1, h2)
    return not (h1 <= h_new <= h2)

def index_in_range(Y, condition, k_importance=1/10):
    '''
        Auxiliar function for pruning the states within the set of integration steps
    '''

    index = int(len(Y)*k_importance)
    if not condition(Y[index]):
        return len(Y)

    for i in range(index, len(Y)):
        if not condition(Y[i]):
            return i







