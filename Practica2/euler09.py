import numpy as np
import sys
import copy
import random

"""
0.0 Auxiliary Class: Multigraph.
"""


class Multigraph:
    """
    Class for multigraph representation.
    ___________________________________
    Its main purpose is to store all the variables needed
    during multigraph generation.
    """

    def __init__(self, n_nodes, decimals, num_max_multiple_edges, max_weight):
        """
        Constructor method for a multigraph with given parameters.

        Parameters
        ----------
            n_nodes: Total number of nodes of the graph.
            decimals: Decimal digits in weight.
            num_max_multiple_edges: Maximum number of edges for each node.
            max_weight: Maximum weight of an edge
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
        """
        return self.mg

    def get_edges(self, node):
        """
        Gets the edges from a given node

        Parameters
        ----------
             node: given node

        Returns edges from given node.
        """
        return self.mg[node]

    def set_edge(self, node1, node2, cost):
        """
        Creates an edge between two nodes with a given cost.

        Parameters
        ----------
             node1: Origin node.
             node2: Destination node.
             cost: Edge cost.
        """

        if (len(self.mg[node1]) < 1) or (node2 not in self.mg[node1]):
            self.mg[node1][node2] = {0: cost}
        else:
            self.mg[node1][node2][len(self.mg[node1][node2])] = cost

    def get_num_max_multiple_edges(self):
        """
        Gets the maximum number of edges for each node
        in the current multigraph
        """
        return self.num_max_multiple_edges

    def get_decimals(self):
        """
        Gets the decimal digits number
        """
        return self.decimals

    def get_max_weight(self):
        """
        Gets the maximum weight of the graph.
        """
        return self.max_weight

    def __str__(self):
        return str(self.mg)


"""
0.1 Auxiliary Functions: Random Multigraph Generation.
"""


def edgeGeneration(mg, fl_diag, probability, fl_unweighted, fl_directed):
    """
    Edge generation method, applying probability to computations
    for a given multigraph.

    Parameters
    ----------
        mg: Multigraph structure
        fl_diag: Allow/Deny self-connected edges
        probability: Edge probability
        fl_unweighted: Flag that signals if a graph is or isn't unweighted
    """

    nodes = list(mg.get_nodes().keys())  # Nodes of the graph
    for node in nodes:
        # EDGE SELECTION

        if not fl_diag:  # When disabled self-linked-edges
            nodes.remove(node)

        # NumPy chooses nodes to link
        linked = np.random.choice(list(nodes),
                                  mg.get_num_max_multiple_edges())

        if not fl_diag:  # When disabled self-linked-edges
            nodes.append(node)  # Return the node to list
            nodes.sort()  # Sort the list

        # EDGE CREATION

        for lnode in linked:  # For every selected node
            # Numpy randomly generates edges according to probability
            if np.random.uniform() < probability:
                # Setting weight of the edge

                if fl_unweighted:  # When disabled graph weight
                    weight = 1
                else:
                    # Select a number between 0 and maximum weight with
                    # uniform probability. Truncate to decimals.
                    weight = round(np.random.uniform(0,
                                                     mg.get_max_weight()),
                                   mg.get_decimals())

                mg.set_edge(node, lnode, weight)  # Creates each edge

                if not fl_directed:  # If undirected graph, connect both nodes
                    mg.set_edge(lnode, node, weight)


def rand_weighted_multigraph(n_nodes, probability=0.6,
                             num_max_multiple_edges=2,
                             max_weight=50,
                             decimals=0,
                             fl_unweighted=False,
                             fl_diag=True):
    """
    Function that generates a directed multigraph with given parameters.

    Parameters
    ----------
        n_nodes: Number of nodes of the graph.
        probability: Probability of edges between nodes.
        num_max_multiple_edges: Maximum number of edges for each node.
        max_weight: Maximum weight of an edge.
        decimals: Number of decimal digits in weights.
        fl_unweighted: Flag for classifying weighted/unweighted graph.
        fl_diag: Flag for enabling/disabling self-linked nodes.
    Returns Directed multigraph constructed with given parameters.
    """

    fl_directed = True
    # Creation of the Multigraph Object
    mg = Multigraph(n_nodes, decimals, num_max_multiple_edges, max_weight)

    # Multigraph Random Generation
    edgeGeneration(mg, fl_diag, probability, fl_unweighted, fl_directed)

    mg = mg.get_nodes()  # Gets Multigraph dictionary structure

    return mg


def rand_weighted_undirected_multigraph(n_nodes,
                                        probability=0.2,
                                        num_max_multiple_edges=1,
                                        max_weight=50,
                                        decimals=0,
                                        fl_unweighted=False,
                                        fl_diag=True):
    """
    Function that generates an undirected multigraph with given parameters.

    Parameters
    ----------
        n_nodes: Number of nodes of the graph.
        probability: Probability of edges between nodes.
        num_max_multiple_edges: Maximum number of edges for each node.
        max_weight: Maximum weight of an edge.
        decimals: Number of decimal digits in weights.
        fl_unweighted: Flag for classifying weighted/unweighted graph.
        fl_diag: Flag for enabling/disabling self-linked nodes.
    Returns Undirected multigraph constructed with given parameters.
    """
    mg = Multigraph(n_nodes, decimals, num_max_multiple_edges, max_weight)
    fl_directed = False
    edgeGeneration(mg, fl_diag, probability, fl_unweighted, fl_directed)
    return mg.get_nodes()


"""
0.2 Auxiliary Functions: Adjacency and Incidency Lists.
"""


def print_adj_list_mg(mg):
    """
    Print the adjacency list for a given multigraph

    Parameters
    ----------
        mg: Multigraph to print
    """

    for key, value in mg.items():
        print(str(key), end="")
        for destiny in sorted(value.keys()):
            print(" - {}".format(destiny), end="")
        print("")
    print("")


"""
0.3 Auxiliary Functions: DFS, Connectivity and Adjacency/Incidency lists.
"""


def dfs(mg, visited, node):
    """
    Depth-first search for a given graph and initial node

    Parameters
    ----------
        mg: Multigraph.
        visited: List for visited nodes.
        node: Initial node.
    Mark as visited all nodes traversed during DFS.
    """
    visited[node] = True
    for i in mg[node]:
        if visited[i] == False:
            dfs(mg, visited, i)


def connected_component(mg, connected_nodes):
    """
    Function that checks if a given graph is connected

    Parameters
    ----------
        mg: Multigraph.
    Returns true or false depending on connected condition of the graph.
    """
    length = len(mg)
    visited = [False] * length

    # Component traversal.
    dfs(mg, visited, connected_nodes[0])

    # Connected component check.
    for node in connected_nodes:
        if not visited[node]:
            return False

    # All nonzero in-out nodes verified as connected component.
    return True


def inc_adj_analysis(gsize, adj, inc, nonzero=None):
    """
    Mathematical check function for eulerian paths.

    Parameters
    ----------
        gsize: Graph size.
        adj: Adjacency list of the graph.
        inc: Incidency list of the graph.
    Returns true or false depending on equalities in .
    """
    nodeA = None
    nodeB = None

    # Adjacency/Indidency vertex analysis
    for i in range(gsize):
        adji = len(adj[i])
        inci = len(inc[i])

        # Connected component check storage
        if nonzero != None:
            if inci != 0 or adji != 0:
                nonzero.append(i)

        # Every other vertex in == out
        if nodeA and nodeB and inci != adji:
            return False

        # 1 more in than out
        if inci+1 == adji:
            if nodeA:
                return False
            nodeA = i

        # 1 more out than in
        if inci == adji+1:
            if nodeB:
                return False
            nodeB = i

    return True


def dfs_lastnode(mg, visited, node, path):
    """
    Depth-first search for a given graph and initial node
    This version gets the last visited node

    Parameters
    ----------
        mg: Multigraph.
        visited: List for visited nodes.
        node: Initial node.
    Mark as visited all nodes traversed during DFS.
    """
    visited[node] = True
    path.append(node)
    for i in mg[node]:
        if visited[i] == False:
            dfs_lastnode(mg, visited, i, path)


""" TODO I. Guardando y Leyendo grafos """


def d_g_2_TGF(d_g, f_name):
    """
    Function that stores a given multigraph into a file.

    Parameters
    ----------
        d_g: Multigraph.
        f_name: Destination file name.
    Returns the multigraph into the file.
    """
    n_nodes = range(len(d_g))
    print(d_g)
    original_stdout = sys.stdout

    with open(f_name, 'w+') as f:
        sys.stdout = f
        for node in n_nodes:
            print(node)
        print("#")
        for adjacent_nodes, i in zip(d_g, n_nodes):
            for node in adjacent_nodes:
                print
                print(i, node)
        sys.stdout = original_stdout


def TGF_2_d_g(f_name):

    original_stdout = sys.stdout

    with open(f_name, 'r') as f:
        text = f.read()

    text = text.replace("\n", "")
    text = text.replace(" ", "")
    text = text.replace("\n", "")

    i = 0
    d_g = []
    n_nodes = 0
    nodeflag = True

    for line in text:

        if(line == '#'):
            nodeflag = False

            for _ in range(n_nodes):
                d_g.append([])

            print(d_g)

        if nodeflag:
            n_nodes += 1

        elif line != '#':
            print(line)


"""
II. Caminos eulerianos en multigrafos dirigidos
"""


def adj_inc_directed_multigraph(mg):
    """
    Function that returns adjacency and incidency lists for every vertex in mg

    Parameters
    ----------
        mg: Multigraph.
    Returns Adjacency and Incidency lists for every vertex.
    """

    adj = []
    inc = [[] for i in mg]

    for origin, edges in mg.items():
        adj_temp = []
        for destiny in sorted(edges):
            adj_temp.append(destiny)
            inc[destiny].append(origin)
        adj.append(adj_temp)

    return inc, adj


def isthere_euler_path_directed_multigraph(d_mg):
    """
    Function that checks if there is an eulerian path in the multigraph.

    Parameters
    ----------
        d_mg: Multigraph.
    Returns True or False depending on eulerian path existance in mg.
    """
    gsize = len(d_mg)
    inc, adj = adj_inc_directed_multigraph(d_mg)
    nonzero = []

    # Incidency and Adjacency list analysis
    if not inc_adj_analysis(gsize, adj, inc, nonzero):
        return False

    # Finding connected components: Component connectivity
    if not connected_component(d_mg, nonzero):
        return False

    # All conditions satisfied
    return True


def first_last_euler_path_directed_multigraph(d_mg):
    """
    Function that finds initial and final node of an eulerian path in mg.

    Parameters
    ----------
        d_mg: Multigraph to analyze.
    Returns initial and final nodes as tuple of an eulerian path.
    """
    nonzero = None
    last = None

    # Check if there is any eulerian path
    if not isthere_euler_path_directed_multigraph(d_mg):
        return ()

    # Get Incidency and Adjacency Lists
    inc, adj = adj_inc_directed_multigraph(d_mg)

    # Get first item in the connected component
    # Also check if we can find any last vertex
    for i in range(len(d_mg)):
        inci = len(inc[i])
        adji = len(adj[i])

        # Any initial in connected component
        if (inci != 0 or adji != 0) and nonzero == None:
            nonzero = i

        # Neccesary last vertex
        if inci == adji + 1 and last == None:
            last = i

        if nonzero != None and last != None:
            return (nonzero, last)

    # Traverse as Eulerian Walk
    mg = copy.deepcopy(d_mg)
    last = euler_walk_directed_multigraph(nonzero, mg)[-1]

    # Return Initial and Final nodes
    return (nonzero, last)


def euler_walk_directed_multigraph(u, d_mg):
    """
    Function that eulerianly traverses multigraphs.

    Parameters
    ----------
        u: Initial node to begin the path
        d_mg: Multigraph to work with.
    Returns traversed nodes of the eulerianly traversed path.
    """

    current_node = u
    visited = []

    while True:
        edges = list(d_mg[current_node])
        if not edges:
            break
        d_mg[current_node].pop(edges[0])
        visited.append(current_node)
        current_node = edges[0]

    visited.append(current_node)
    return visited


def next_first_node(l_path, d_mg):
    """
    Function that finds next node to reset eulerian traversal.

    Parameters
    ----------
        l_path: Traversed nodes during eulerian walk.
        d_mg: Multigraph to work with.
    Returns node to reset eulerian traversal.
    """

    for node in l_path:
        if d_mg[node]:
            return node

    return None


def path_stitch(path1, path2):
    """
    Function that pastes paths .

    Parameters
    ----------
        path1: Main path to be pasted into.
        path2: Additional path to be pasted into the main one.
    Returns joined paths.
    """
    stick = path2[0]
    ind = path1.index(stick)
    result = path1[:ind+1] + path2[1:-1] + path1[ind:]
    return result


def remaining_edges(mg):
    """
    Function that checks if there is any remaining edge in the multigraph.

    Parameters
    ----------
        mg: Multigraph.
    Returns the boolean evaluation of edge existance condition.
    """
    for edges in mg.values():
        if edges:
            return True

    return False


def euler_path_directed_multigraph(d_mg):
    """
    Function that tries to find an euler path in a given directed multigraph.

    Parameters
    ----------
        d_mg: Directed multigraph.
    Returns an euler path for the given multigraph if possible.
    """
    mg = copy.deepcopy(d_mg)
    visited = euler_walk_directed_multigraph(0, mg)
    origin = 0

    while remaining_edges(mg):
        origin = next_first_node(visited, mg)
        tmp = euler_walk_directed_multigraph(origin, mg)
        visited = path_stitch(visited, tmp)

    return visited


def vertex_degree_check(adj, inc, nonzero):
    """
    Function that checks adjacency and incidency equality as stores component.

    Parameters
    ----------
        d_mg: Directed multigraph.
    Returns boolean evaluation for the equality condition.
    """
    for i in range(len(adj)):
        inci = len(inc[i])
        adji = len(adj[i])

        # Connected component check storage
        if adji != 0 or inci != 0:
            nonzero.append(i)

        # Incidency must be equal to Adjacency
        if adji != inci:
            return False

    return True


def isthere_euler_circuit_directed_multigraph(d_mg):
    """
    Function that tries to find an euler circuit in a given multigraph.

    Parameters
    ----------
        d_mg: Directed multigraph.
    Returns boolean evaluation for euler path existance for given multigraph.
    """
    adj, inc = adj_inc_directed_multigraph(d_mg)
    visited = [False] * len(d_mg)
    nonzero = []

    # Incidency and Adjacency Equality
    if not vertex_degree_check(adj, inc, nonzero):
        return False

    # Connected Component: Depth first traversal
    dfs(d_mg, visited, visited[nonzero[0]])

    # Connected Component: Verification
    for node in nonzero:
        if not visited[node]:
            return False

    # Conditions satisfied
    return True


def euler_circuit_directed_multigraph(d_mg, u=0):
    """
    Function that tries to find an eulerian circuit in given multigraph.

    Parameters
    ----------
        d_mg: Size of the desired sequence
        u: Initial circuit vertex
    Returns euler circuit in given multigraph if existing, otherwise None.
    """
    if not isthere_euler_circuit_directed_multigraph(d_mg):
        return None

    circuit = euler_path_directed_multigraph(d_mg)

    return circuit


def random_sequence(len_seq):
    """
    Function that generates a protein sequence of selected length.

    Parameters
    ----------
        len_seq: Size of the desired sequence
    Returns random protein sequence.
    """
    proteins = ['A', 'C', 'G', 'T']
    sequence = []
    for _ in range(len_seq):
        sequence.append(random.choice(proteins))
    return sequence


def spectrum(sequence, len_read):
    """
    Function that returns unordered spectrum of a sequence.

    Parameters
    ----------
        sequence: Input sequence
        len_read: Read size
    Returns unordered spectrum of a sequence.
    """
    spectrum = []

    for i in range(len(sequence) - len_read + 1):
        spectrum.append(sequence[i:i+len_read])

    return spectrum


def spectrum_2(spectr):
    """
    Function that returns the l-1 spectrum associated to the l-spectrum.

    Parameters
    ----------
        spectr: l-spectrum 
    Returns Associated (l-1) spectrum.
    """
    spectrum = []
    length = len(spectr[0])
    for proteins in spectr:
        for i in range(length-1):
            prot = proteins[i:i+length-1]
            if prot not in spectrum:
                spectrum.append(prot)
    return spectrum


def spectrum_2_graph(spectr):
    """
    Function that returns a multigraph associated to l-spectrum.

    Parameters
    ----------
        spectr: l-spectrum 
    Returns Multigraph associated to the l-spectrum.
    """
    n_nodes = len(spectr)
    mg = Multigraph(n_nodes, 0, n_nodes, 1)

    # Iterates every item in l-spectrum as node
    for i in range(n_nodes):
        protein = spectr[i]
        # Iterates every possible connection
        for j in range(n_nodes):
            destiny = spectr[j]
            # Connects
            if protein != destiny and protein[-1:] == destiny[:1]:
                mg.set_edge(i, j, 1)

    # Multigraph generated
    return mg.get_nodes()


def spectrum_2_sequence(spectr):
    mg = spectrum_2_graph(spectr)
    path = euler_path_directed_multigraph(mg)
    return path


def path_2_sequence(l_path, spectrum_2):
    pass

# Si devuelve booleano, devuelve siempre verdadero, nunca falso por propiedades de las secuencias
# Coherencia tienen que coincidir los espectros:
# 1. generar secuencia aleatoria
# 2. generar su espectro
# 3. transformar espectro en grafo
# 4. en ese grafo aplicar funciones eulerianas para ver si hay circuito o camino
# 5. si lo hay, aplicar funciones de la segunda parte de la práctica
# 6. la reconstruccion es coherente si coincide con el espectro de la SECUENCIA INICIAL.


def check_sequencing(len_seq, len_read):
    pass


g = {
    0: {1: {0: 49.0}},
    1: {2: {0: 17.0}},
    2: {3: {0: 20.0}},
    3: {0: {0: 30.0}, 4: {0: 10.0}},
    4: {3: {0: 10.0}}
}


# EULER PATH CHECK
booleaneuler = isthere_euler_path_directed_multigraph(g)
# print(booleaneuler)
inc, adj = adj_inc_directed_multigraph(g)

# FILE SAVE TODO
#d_g_2_TGF(adj, "save.d")
# TGF_2_d_g("save.d")

# euler init and end
test = first_last_euler_path_directed_multigraph(g)
# print(test)

# euler walk
g_copy = copy.deepcopy(g)
variable = euler_walk_directed_multigraph(0, g_copy)
# print(variable)

var2 = next_first_node(variable, g_copy)
# print(var2)

# path stitch
list1 = [0, 1, 2, 3, 0]
list2 = [3, 4, 3]
resultado = path_stitch(list1, list2)
# print(resultado)

asd = euler_path_directed_multigraph(g)
# print(asd)

asd = isthere_euler_circuit_directed_multigraph(g)
# print(asd)

asd = euler_circuit_directed_multigraph(g)
# print(asd)

sequence = random_sequence(10)
# print(sequence)
spectrum = spectrum(sequence, 3)
#print(spectrum)
spectrum2 = spectrum_2(spectrum)
#print(spectrum2)

test = spectrum_2_graph(spectrum2)
#print(test)


print(spectrum2)
print(test)
test = spectrum_2_sequence(spectrum2)
print(test)

