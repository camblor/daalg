from queue import PriorityQueue
import time
import numpy as np

"""
II. MULTI-GRAPHS
"""


class Multigraph:
    """
    Class for multigraph representation.
    Its main purpose is to store all the variables needed during multigraph generation.
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

    def __str__(self):
        return str(self.mg)


def edgeGeneration(mg, fl_diag, probability, fl_unweighted, fl_directed):
    """
    Edge generation method, applying probability to computations for a given multigraph.
    :param mg: Multigraph structure
    :param fl_diag: Allow/Deny self-connected edges
    :param probability: Edge probability
    :param fl_unweighted: Flag that signals if a graph is or isn't unweighted
    :param fl_directed: Flag that signals if a graph is or isn't undirected
    """

    nodes = list(mg.get_nodes().keys())  # Nodes of the graph
    for node in nodes:
        # EDGE SELECTION

        if not fl_diag:  # When disabled self-linked-edges
            nodes.remove(node)

        linked = np.random.choice(list(nodes), mg.get_num_max_multiple_edges())  # NumPy chooses nodes to link

        if not fl_diag:  # When disabled self-linked-edges
            nodes.append(node)  # Return the node to list
            nodes.sort()  # Sort the list

        # EDGE CREATION

        for lnode in linked:  # For every selected node

            if np.random.uniform() < probability:  # Numpy randomly generates edges according to probability
                # Setting weight of the edge

                if fl_unweighted:  # When disabled graph weight
                    weight = 1
                else:
                    # Select a number between 0 and maximum weight with uniform probability. Truncate to decimals.
                    weight = round(np.random.uniform(0, mg.get_max_weight()), mg.get_decimals())

                mg.set_edge(node, lnode, weight)  # Creates each individual edge

                if not fl_directed:  # If undirected graph, connect both nodes.
                    mg.set_edge(lnode, node, weight)


def rand_weighted_multigraph(n_nodes, probability=0.6, num_max_multiple_edges=2, max_weight=50, decimals=0,
                             fl_unweighted=False, fl_diag=True):
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

    mg = Multigraph(n_nodes, decimals, num_max_multiple_edges, max_weight)  # Creation of the Multigraph Object

    edgeGeneration(mg, fl_diag, probability, fl_unweighted, fl_directed)  # Multigraph Random Generation

    mg = mg.get_nodes()  # Gets Multigraph dictionary structure

    return mg


def rand_weighted_undirected_multigraph(n_nodes, probability=0.2, num_max_multiple_edges=1, max_weight=50,
                                        decimals=0, fl_unweighted=False, fl_diag=True):
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

    for key, value in mg.items():
        print(str(key), end="")
        for destiny in sorted(value.keys()):
            print(" - {}".format(destiny), end="")
        print("")
    print("")


"""
III. Minimum distances in Multigraphs.
"""


def dijkstra_mg(mg, u):
    """
    Function that applies Dijkstra algorithm for minimum distance computation from a node.
    :param mg: Multigraph (Dict of Dict of Dict).
    :param u: Selected node.
    :return: Distance and Previous node dictionaries
    """

    # Previous node initialization
    previous = {u: u}

    # Distance dictionary creation and initialization
    distance = {v: np.inf for v in mg}
    distance[u] = 0

    # Visited dictionary creation and initialization
    visited = {v: False for v in mg}

    # Priority queue object creation.
    q = PriorityQueue()
    q.put((0, u))

    # Dijkstra algorithm implementation.
    while not q.empty():

        _, current = q.get()  # Obtains next node from the queue.
        if not visited[current]:  # Checks if we need to visit this node.

            visited[current] = True  # Marks the node as visited.

            for dest, routes in mg[current].items():  # Iterates through every destiny from current node.

                for route in routes:  # Iterates through every edge to destiny.

                    # Stores direct and traveling distance for comparison
                    current_distance = distance[dest]
                    traveling_distance = distance[current] + mg[current][dest][route]

                    # If (destiny node has not been visited + current distance is higher than the traveling distance)
                    if visited[dest] is False and current_distance > traveling_distance:
                        distance[dest] = traveling_distance  # Update distance to the traveling distance

                        previous[dest] = current  # Update previous node to destination as current node.

                        q.put((distance[dest], dest))  # Inserts next node for Dijkstra.

    return distance, previous


def min_paths(d_prev):
    """
    Calculation of the path to traverse from the first node to every other node in the graph.
    :param d_prev: Dictionary with previous nodes required to traverse to get to selected node.
    :return: Dictionary with lists of shortest path to every node from the origin.
    """
    d_path = {}  # Dictionary Initialization

    for vertex in d_prev:
        previous = [vertex, d_prev[vertex]]  # Get the previous node in the path to the selected node.
        i = 1
        while previous[i] is not 0:  # While we haven't arrived the origin.
            previous.append(d_prev[previous[i]])  # Add the node we are traversing to the path.
            i += 1

        previous.reverse()  # Reverse the path for data representation.
        d_path[vertex] = previous  # Store the data into our dictionary

    return d_path


def time_dijkstra_mg(n_graphs, n_nodes_ini, n_nodes_fin, step, prob=0.2):
    """
    Generated a given number of described graphs and computes Dijkstra algorithm for every node for every graph.
    :param n_graphs: Number of graphs used in the measurement.
    :param n_nodes_ini: Initial value for the number of nodes interval.
    :param n_nodes_fin: Final value for the number of nodes interval.
    :param step: Step between numbers of nodes in the interval.
    :param prob: Probability to generate edges between two given nodes.
    :return: Mean of the time elapsed in Dijkstra computations.
    """
    countertime = []  # Storage for times

    for n_nodes in np.arange(n_nodes_ini, n_nodes_fin + 1, step):  # Simulation with given interval
        counter = 0
        for graph_number in range(n_graphs):  # Generate n_graphs to measure
            mg = rand_weighted_multigraph(n_nodes, probability=prob)  # Graph generation

            start = time.time()  # Start time counter
            for node in mg.keys():  # Compute Dijkstra for every node in the graph.
                dijkstra_mg(mg, node)  # Compute Dijkstra for this node.
            end = time.time()  # End time counter

            counter += end - start  # Store result

        countertime.append(counter / (n_graphs * n_nodes))  # Save for this number of nodes

    return countertime


"""
IV. Dijkstra vs Floyd-Warshall
"""


def dijkstra_all_pairs(g):
    """
    Generates a matrix with all the results from Dijkstra Algorithm for every node in the given graph
    :param g: Given graph
    :return: Matrix with distances between two nodes represented in cells.
    """
    n = len(g)  # Get |g|
    result = np.full((n, n), np.inf)  # Generates the NumPy Matrix

    for node in np.arange(n):  # Traverses node index

        distance, _ = dijkstra_mg(g, node)  # Gets the distance

        for destination in np.arange(n):  # Traverses possible destinations
            result[node][destination] = distance[destination]  # Saves the distance in the correct cell.

    return result


def dg_2_ma(g):
    """
    Computes the adjacency matrix of a given graph.
    :param g: Given graph
    :return: Adjacency matrix of the graph
    """

    n = len(g)  # Get |g|

    adjmat = np.full((n, n), np.inf)  # Prepare matrix with NumPy

    for node in np.arange(n):  # Iterates node index
        adjmat[node][node] = 0  # Sets the value to same node to 0.

        for destination, cost in g[node].items():  # Iterates through destinations
            adjmat[node][destination] = cost[0]  # Updates the matrix with destination cost.

    return adjmat


def floyd_warshall(ma_g):
    """
    Computes Floyd-Warshall algorithm for a given adjacency matrix of a graph.
    :param ma_g: Adjacency Matrix
    :return: Matrix with minimum costs.
    """

    result = np.copy(ma_g)  # Save matrix variable
    n = len(ma_g)  # Store number of rows (square matrix)

    for k in range(n):  # For every 1D
        for i in range(n):  # For every 2D
            for j in range(n):  # For every 3D
                result[i][j] = min(result[i][j], result[i][k] + result[k][j])  # Compare the cost with other paths

    return result


def time_dijkstra_mg_all_pairs(n_graphs, n_nodes_ini, n_nodes_fin, step, num_max_multiple_edges=1, probability=0.5):
    """
    Generated a given number of described graphs and computes Dijkstra All-Pairs for every node in every graph.
    :param n_graphs: Number of graphs used in the measurement.
    :param n_nodes_ini: Initial value for the number of nodes interval.
    :param n_nodes_fin: Final value for the number of nodes interval.
    :param step: Step between numbers of nodes in the interval.
    :param num_max_multiple_edges: Maximum number of edges per node.
    :param probability: Probability to generate edges between two given nodes.
    :return: Mean of the time elapsed in Dijkstra computations.
    """
    countertime = []  # Storage for times
    for n_nodes in np.arange(n_nodes_ini, n_nodes_fin + 1, step):  # Simulation with given interval

        counter = 0  # Time counter variable
        for graph_number in range(n_graphs):  # Generate n_graphs to measure
            mg = rand_weighted_multigraph(n_nodes, num_max_multiple_edges=num_max_multiple_edges,
                                          probability=probability)  # Graph generation

            start = time.time()  # Start time counter
            dijkstra_all_pairs(mg)  # Compute Dijkstra for this node.
            end = time.time()  # End time counter

            counter += end - start  # Store time difference

        countertime.append(counter / n_graphs)  # Append to given size times.

    return countertime


def time_floyd_warshall(n_graphs, n_nodes_ini, n_nodes_fin, step, probability=0.5):
    """
    Generated a given number of described graphs and computes Dijkstra All-Pairs for every node in every graph.
    :param n_graphs: Number of graphs used in the measurement.
    :param n_nodes_ini: Initial value for the number of nodes interval.
    :param n_nodes_fin: Final value for the number of nodes interval.
    :param step: Step between numbers of nodes in the interval.
    :param probability: Probability to generate edges between two given nodes.
    :return: Mean of the time elapsed in Dijkstra computations.
    """
    countertime = []  # Storage for times
    for n_nodes in np.arange(n_nodes_ini, n_nodes_fin + 1, step):  # Simulation with given interval

        counter = 0  # Time counter variable
        for graph_number in range(n_graphs):  # Generate n_graphs to measure
            mg = rand_weighted_multigraph(n_nodes, probability=probability)  # Graph generation

            adjmatrix = dg_2_ma(mg)  # Obtain the adjacency matrix.

            start = time.time()  # Start time counter
            floyd_warshall(adjmatrix)  # Execute Floyd Warshall to the adjacency matrix.
            end = time.time()  # End time counter

            counter += end - start  # Store time difference

        countertime.append(counter / n_graphs)  # Append to given size times.

    return countertime


if __name__ == '__main__':
    # GRAFOS
    # graph = rand_weighted_multigraph(10)
    # print_adj_list_mg(graph)
    graph = rand_weighted_multigraph(10)
    matrix = dg_2_ma(graph)
    print(matrix)
    res = floyd_warshall(matrix)
    print(res)
    """
    
    print("________________")
    graph = rand_weighted_multigraph(10, fl_diag=False)
    print_adj_list_mg(graph)
    print("________________")
    graph = rand_weighted_undirected_multigraph(10)
    print_adj_list_mg(graph)
    print("________________")
    graph = rand_weighted_undirected_multigraph(10, fl_diag=False)
    # print_adj_list_mg(graph)

    # DISTANCIAS MÍNIMAS"""
    # dist, prev = dijkstra_mg(graph, 0)
    """
    paths = min_paths(prev)
    # print(paths)
    time_dijkstra_mg(1, 500, 1000, 100)

    # IV. Dijkstra vs Floyd-Warshall
    matrix = dijkstra_all_pairs(graph)
    # print(matrix)
    adjacencymat = dg_2_ma(graph)
    print(adjacencymat)
    testing = floyd_warshall(adjacencymat)
    print(testing)
    """
