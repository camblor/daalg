from queue import PriorityQueue
import time
import numpy as np

"""
II. MULTI-GRAPHS
"""


class AdjNode:
    """
    Class for adjacency nodes implementation
    """

    def __init__(self, value):
        """
        Constructor method for adjacency nodes
        :param value: Value stored in the node
        """
        self.vertex = value
        self.next = None


class Multigraph:
    """
    Class for multigraph representation
    """

    def __init__(self, n_nodes, decimals, num_max_multiple_edges, max_weight):
        """
        Constructor method for a multigraph with given parameters.
        :param n_nodes: Total number of nodes of the graph.
        :param decimals: Decimal digits in weight.
        :param num_max_multiple_edges: Maximum number of edges for each node.
        :param max_weight: Maximum weight of an edge
        """

        # Multigraph Creation
        self.mg = {}
        for node in range(n_nodes):
            self.mg[node] = {}

        # Setting of parameters
        self.num_max_multiple_edges = num_max_multiple_edges
        self.max_weight = max_weight

        # Adjacency list creation
        self.adjlist = {}

        # Decimal digits
        self.decimals = decimals

    def get_nodes(self):
        """
        Get list of nodes
        :return: list of nodes
        """
        return self.mg

    def get_edges(self, node):
        """
        Gets the edges from a given node
        :param node: given node
        :return: edges from given node
        """
        return self.mg[node]

    def set_edge(self, node1, node2, cost):
        """
        Creates an edge between two nodes with a given cost.
        :param node1: Origin node.
        :param node2: Destination node.
        :param cost: Edge cost.
        """

        if len(self.mg[node1]) < 1 or not node2 in self.mg[node1]:
            self.mg[node1][node2] = {0: cost}
        else:
            self.mg[node1][node2][len(self.mg[node1][node2])] = cost

    def get_num_max_multiple_edges(self):
        """
        Gets the maximum number of edges for each node in the current multigraph
        :return: Maximum number of edges for each node.
        """
        return self.num_max_multiple_edges

    def get_decimals(self):
        """
        Gets the decimal digits number
        :return: Decimal digits number
        """
        return self.decimals

    def get_max_weight(self):
        """
        Gets the maximum weight of the graph.
        :return: Maximum Weight.
        """
        return self.max_weight

    def generate_adjlist(self):
        """
        Generates the adjacency list
        :return: adjacency list
        """

        keys = self.get_nodes().keys()

        # Iterates for every node
        for key in keys:
            # Creates the adjacency node list
            self.adjlist[key] = {}
            # Appends the connections into the list
            for vertex in set(self.get_edges(key).keys()):
                # Base case
                if not self.adjlist[key]:
                    temp1 = AdjNode(vertex)
                    self.adjlist[key] = temp1
                # General case
                else:
                    temp2 = AdjNode(vertex)
                    temp1.next = temp2
                    temp1 = temp2

        return self.adjlist

    def __str__(self):
        return str(self.mg)


def edgeGeneration(mg, fl_diag, probability, fl_unweighted, fl_directed):
    """
    Edge generation method, applying probability for computations, for a given multigraph.
    :param mg: Multigraph structure
    :param fl_diag: Allow/Deny self-connected edges
    :param probability: Edge probability
    :param fl_unweighted: Flag that signals if a graph is or isn't unweighted
    :param fl_directed: Flag that signals if a graph is or isn't undirected
    """

    # Nodes of the graph
    nodes = mg.get_nodes().keys()

    for node in nodes:
        # Enabled/Disabled Self-linked edges
        if fl_diag:
            # NumPy chooses nodes to link from current node
            linked = np.random.choice(list(nodes), mg.get_num_max_multiple_edges())
        else:
            availablenodes = [x for x in nodes if x != node]
            # NumPy chooses nodes to link from current node
            linked = np.random.choice(availablenodes, mg.get_num_max_multiple_edges())

        # Edges creation
        for lnode in linked:

            # Applying probability to edge generation
            if np.random.uniform(0, 1) < probability:

                # Setting weight of the edge
                if fl_unweighted:
                    weight = 1
                else:
                    step = 1 / (10 ** mg.get_decimals())
                    decimal_range = np.arange(0, mg.get_max_weight(), step)
                    weight = round(np.random.choice(decimal_range), mg.get_decimals())

                # Creates the edges
                mg.set_edge(node, lnode, weight)

                # If undirected graph, connect both nodes.
                if not fl_directed:
                    mg.set_edge(lnode, node, weight)


def generate_adjlist(mg):
    """
    Generation for the adjacency list for a given multigraph.
    :param mg: Multigraph.
    :return: Adjacency list of the multigraph.
    """
    adjlist = {}
    # Iterates for every node
    for key in mg:
        # Creates the adjacency node list
        adjlist[key] = {}
        # Appends the connections into the list
        for vertex in set(mg[key]):
            # Base case
            if not adjlist[key]:
                temp1 = AdjNode(vertex)
                adjlist[key] = temp1
            # General case
            else:
                temp2 = AdjNode(vertex)
                temp1.next = temp2
                temp1 = temp2
    return adjlist


def rand_weighted_multigraph(n_nodes,
                             probability=0.2,
                             num_max_multiple_edges=1,
                             max_weight=50,
                             decimals=0,
                             fl_unweighted=False,
                             fl_diag=True):
    """
    Function that generates a directed multigraph with given parameters.
    :param n_nodes: Number of nodes of the graph.
    :param probability: Probability of edges between nodes.
    :param num_max_multiple_edges: Maximum number of edges for each node.
    :param max_weight: Maximum weight of an edge.
    :param decimals: Number of decimal digits in weights.
    :param fl_unweighted: Flag for classifying weighted/unweighted graph.
    :param fl_diag: Flag for enabling/disabling self-linked nodes.
    :return: Directed multigraph constructed with given parameters.
    """
    fl_directed = True
    mg = Multigraph(n_nodes, decimals, num_max_multiple_edges, max_weight)
    edgeGeneration(mg, fl_diag, probability, fl_unweighted, fl_directed)
    return mg.get_nodes()


def rand_weighted_undirected_multigraph(n_nodes,
                                        probability=0.2,
                                        num_max_multiple_edges=1,
                                        max_weight=50,
                                        decimals=0,
                                        fl_unweighted=False,
                                        fl_diag=True):
    """
    Function that generates an undirected multigraph with given parameters.
    :param n_nodes: Number of nodes of the graph.
    :param probability: Probability of edges between nodes.
    :param num_max_multiple_edges: Maximum number of edges for each node.
    :param max_weight: Maximum weight of an edge.
    :param decimals: Number of decimal digits in weights.
    :param fl_unweighted: Flag for classifying weighted/unweighted graph.
    :param fl_diag: Flag for enabling/disabling self-linked nodes.
    :return: Undirected multigraph constructed with given parameters.
    """
    mg = Multigraph(n_nodes, decimals, num_max_multiple_edges, max_weight)
    fl_directed = False
    edgeGeneration(mg, fl_diag, probability, fl_unweighted, fl_directed)
    return mg.get_nodes()


def print_adj_list_mg(mg):
    """
    Print the adjacency list for a given multigraph
    :param mg: Multigraph to print
    """

    adjlist = generate_adjlist(mg)

    for key in adjlist:
        print(str(key) + ":", end="")
        temp = adjlist[key]
        while temp:
            print(" -> {}".format(temp.vertex), end="")
            temp = temp.next
        print("")


"""
III. Distancias Mínimas en Multigrafos.
"""


# TODO queue = PriorityQueue()

def dijkstra_mg(mg, u):
    """
    Function that applies Dijkstra algorithm for minimum distance computation from a node.
    :param mg: Multigraph (Dict of Dict of Dict).
    :param u: Selected node.
    :return: Distance and Previous node dictionaries
    """
    previous = {u: u}

    distance = {v: np.inf for v in mg.keys()}
    distance[u] = 0
    visited = {v: False for v in mg.keys()}

    q = PriorityQueue()
    q.put((0, u))

    while not q.empty():
        _, v = q.get()
        if not visited[v]:
            visited[v] = True
            for z in mg[v].keys():
                for k in mg[v][z]:
                    if visited[z] is False and distance[z] > distance[v] + mg[v][z][k]:
                        distance[z] = distance[v] + mg[v][z][k]
                        previous[z] = v
                        q.put((distance[z], z))

    return distance, previous


def min_paths(d_prev):
    """
    Calculation of the path to traverse from the first node to every other node in the graph.
    :param d_prev: List with previous nodes required to traverse to get to selected node.
    :return: Dictionary with lists of shortest path to every node from the origin.
    """
    d_path = {}

    for vertex in d_prev:
        previous = [vertex, d_prev[vertex]]
        i = 1

        while previous[i] is not 0:
            previous.append(d_prev[previous[i]])
            i += 1

        previous.reverse()
        d_path[vertex] = previous

    return d_path


def time_dijkstra_mg(n_graphs,
                     n_nodes_ini,
                     n_nodes_fin,
                     step,
                     prob=0.2):
    # Storage for times
    countertime = {}

    # Generate n_graphs to measure
    for graph_number in range(n_graphs):

        # Storage of measurements for every node count.
        nodes_time = {}
        for n_nodes in range(n_nodes_ini, n_nodes_fin, step):
            # Graph generation
            mg = rand_weighted_multigraph(n_nodes, probability=prob)

            # Time measurement for every vertex
            counter = 0
            for vertex in mg.keys():
                start = time.time()
                dijkstra_mg(mg, vertex)
                end = time.time()
                counter += (end-start)

            # Mean of Dijkstra time for each vertex
            nodes_time[n_nodes] = counter / len(mg)
        # Save for this graph size
        countertime[graph_number] = nodes_time

    print(countertime)


if __name__ == '__main__':
    # GRAFOS
    """
    graph = rand_weighted_multigraph(10)
    print_adj_list_mg(graph)
    print("________________")
    graph = rand_weighted_multigraph(10, fl_diag=False)
    print_adj_list_mg(graph)
    print("________________")
    graph = rand_weighted_undirected_multigraph(10)
    print_adj_list_mg(graph)
    print("________________")"""
    graph = rand_weighted_undirected_multigraph(10, fl_diag=False)
    # print_adj_list_mg(graph)

    # DISTANCIAS MÍNIMAS
    dist, prev = dijkstra_mg(graph, 0)

    paths = min_paths(prev)
    # print(paths)
    time_dijkstra_mg(5, 1000, 1300, 100)
