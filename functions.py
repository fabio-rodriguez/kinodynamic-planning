from math import *
from model import tc, ode_dydt, c

def path_root_to_node(tree, node_tag):
    '''
        Given:
            :param tree: tree of states
            :param node_tag: tag of the target node in the tree

        Return the branch of the tree, from the root until a target node with tag "node_tag"
    '''

    node_tag = "root" if node_tag == "Root" else node_tag
    seq = node_tag.split("_")
    seq.pop()
    node = tree.get_node(node_tag)

    if not len(seq):
        return [node]
    else:
        return path_root_to_node(tree, "_".join(seq)) + [node]

def line_point_distance(line_point1, line_point2, external_point):
    '''
        Given:
            :param line_point1: 1st  point from a straight line r
            :param line_point2: 2nd point from the straight line r
            :param node_tag: point external to the straight line r

        Return the distance between a straight line defined by "line_point1" and "line_point2",
        and the external point "external_point"
    '''

    m = (line_point1[1] - line_point2[1])/ (line_point1[0] - line_point2[0])
    A, B, C = m, -1, line_point1[1] - m*line_point1[0]

    return abs(A*external_point[0] + B*external_point[1] + C)/sqrt(A**2+B**2)

def all_transitions():
    '''
        Set M of 28 transitions considered combining tail angles and frequency values
    '''

    actions = []
    actions += [(round(-radians(x), 5), 0) for x in range(0, 7)]
    actions += [(round(-radians(x), 5), 4*tc) for x in range(0, 7)]
    actions += [(round(-radians(x), 5), 5*tc) for x in range(0, 7)]
    actions += [(round(-radians(x), 5), 6*tc) for x in range(0, 7)]

    return {a: actions for a in actions}

def transitions_reduced():
    '''
        Reduced set Mr of 17 transitions combining tail angles and frequency values
    '''

    actions = [(-0.0, 0), (-0.01745, 0), (-0.03491, 0), (-0.05236, 0), (-0.06981, 0), (-0.08727, 0),
               (-0.10472, 0), (-0.0, 0.12684496748616697), (-0.05236, 0.12684496748616697),
               (-0.06981, 0.12684496748616697), (-0.08727, 0.12684496748616697), (-0.10472, 0.12684496748616697),
               (-0.0, 0.1585562093577087), (-0.05236, 0.1585562093577087), (-0.06981, 0.1585562093577087),
               (-0.0, 0.19026745122925046), (-0.03491, 0.19026745122925046)]
    return {a:actions for a in actions}

def transitions_preperching():
    '''
        Reduced set Mr of 10 transitions for perching, combining tail angles and frequency values
    '''

    actions = [(radians(-1), 0.0), (radians(-2), 0.0), (radians(-3), 0.0), (radians(-4), 0.0), (radians(-5), 0.0),
               (radians(-6), 0.0), (0, 4.0*tc), (0, 5.0*tc), (0, 6.0*tc)]
    return {a:actions for a in actions}

def euclidean_dist(v1, v2):
    '''
        Euclidian distance between vectors v1 and v2
    '''
    return sqrt(sum([(xi0-xi1)**2 for xi0,xi1 in zip(v1, v2)]))
