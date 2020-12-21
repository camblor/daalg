import numpy as np
import time
import math

def get_minimum_power(t):
    """
    Calculation of minimum power of 2 greater than n.

    Parameters
    ----------
        n: input number to calculate with.
    Returns Minimum power of 2 greater than n.
    """
    length = len(t)
    minpow = 2 ** math.ceil(math.log(length, 2))

    # Power of 2 case
    fft_input = t

    # Not power of 2 case: Appending 0's
    if minpow > length:
        for _ in range(minpow-length):
            fft_input = np.append(fft_input, [0])

    return minpow, fft_input

def fft(t):
    """
    Calculation of FFT from a NumPy array.

    Parameters
    ----------
        t: Table.
    Returns Fast Fourier Transform of input.
    """
    # Obtaining table length
    length = len(t)

    # FFT STORAGE VARIABLES
    pares = []
    impares = []

    # FFT RETURN VARIABLE
    result = np.array([])    

    # Minimum power of 2 + 0's appending
    minpow, fft_input = get_minimum_power(t)

    # Not power of 2 case: Appending 0's
    if minpow > length:
        for _ in range(minpow-length):
            fft_input = np.append(fft_input, [0])
        
    # Base case with unique value
    if length <= 1: 
        return fft_input
    
    # Get odd and even elements
    for element, i in zip(fft_input, range(minpow)):
        if i % 2 == 0:
            pares.append(element)
        else:
            impares.append(element)


    # Recursive call
    f_e = fft(pares)
    f_o = fft(impares)
    

    # 2 ^ K-1
    previous2power = int(minpow/2)
    
    for i in range(minpow):
        # First Transformate Operand
        first = f_e[i % previous2power]
        
        #Second Transformate Operand
        tmp1 = (2 * np.pi * 1j) * (i / minpow)
        tmp2 = f_o[i % previous2power]
        second = np.exp(tmp1) * tmp2

        # Transformate result storage
        result = np.append(result, np.round(first, 2) + np.round(second,2 ))

    # Fast Fourier transform
    return result

def invert_fft(t, fft_func=fft):
    """
    Application of inversion algorithm of DFT.

    Parameters
    ----------
        t: Table.
        fft_func: FFT implementation function.
    Returns random protein sequence.
    """
    conjugate = np.conj(t)
    transformate = fft_func(conjugate)
    length = len(transformate)
    output = np.conj(t)
    return output

def rand_polinomio(long=2**10, base=10):
    """
    Random polynomial generation algorithm.

    Parameters
    ----------
        long: Polynomial length indication.
    Returns random polynomial of length long.
    """

    # Base argument restriction
    if base < 2 or base > 10:
        return None

    # Direct return of int generation
    return list(np.random.randint(base, size=long))

def poli_2_num(l_pol, base=10):
    """
    Calculation of resulting int from polynomial evaluation in base base.
    Usage of Horner's rule.

    Parameters
    ----------
        l_pol: Polynomial input.
        base: Base.
    Returns polynomial evaluation in base base.
    """
    tmp=0
    for coef in reversed(l_pol):
        res=int(coef)+tmp
        tmp=res*base
    return res

def rand_numero(num_digits, base=10):
    """
    Random integer generation with num_digits digits.

    Parameters
    ----------
        num_digits: Number of digits of the generated number.
        base: Base.
    Returns random integer with num_digits digits in base base.
    """

    # Random polynomial generation
    polynomial = rand_polinomio(num_digits, base)

    # Transformation from random polynomial to number (random aswell)
    return poli_2_num(polynomial)

def num_2_poli(num, base=10):
    """
    Calculation of integer list with polynomial coefficients equivalent to
    integer num representation in base base. Increasing order.

    Parameters
    ----------
        num: Integer input.
        base: Base.
    Returns integer list with polynomial coefficients.
    """
    poli = []
    while num > 0:
        num, poli_part = divmod(num, base)

        poli.append(poli_part)

    return poli

# TODO listas de enteros de python
def mult_polinomios_fft(l_pol_1, l_pol_2, fft_func=fft):
    """
    Calculation of the product of two polynomials with FFT.

    Parameters
    ----------
        l_pol_1: First polynomial input.
        l_pol_2: Second polynomial input.
        fft_func: FFT function.
    Returns product of input polynomials through fft_funct.
    """

    # Adjust of FFT input
    _, poli1 = get_minimum_power(l_pol_1)
    _, poli2 = get_minimum_power(l_pol_2)

    # Fast Fourier transformates
    coefficient1=fft_func(poli1)
    coefficient2=fft_func(poli2)

    # Coefficient multiplication (SECOND STEP)
    for i in range(len(coefficient1)):
        coefficient1[i] *= coefficient2[i]

    # Inversion algorithm with FFT (FINAL STEP)
    output = invert_fft(coefficient1, fft_func=fft_func)

    # Round and return NumPy Array
    return np.rint(output)

def mult_numeros_fft(num1, num2, fft_func=fft):
    """
    Calculation of the product of two polynomials with FFT.

    Parameters
    ----------
        n_enteros: Number of integer pairs.
        num_digits_ini: First value of the range.
        num_digits_fin: Final value of the range.
        step: Step between selected range.
        fft_func: FFT function.
    Returns product of input polynomials through fft_funct.
    """
    poli_num1 = num_2_poli(num1)
    poli_num2 = num_2_poli(num2)
    return poli_2_num(mult_polinomios_fft(poli_num1, poli_num2))

def time_fft(n_tablas, 
                        num_coefs_ini, 
                        num_coefs_fin, 
                        step, 
                        fft_func=fft):
    """
    Generation of n_tables and mean calculation.

    Parameters
    ----------
        n_tablas: Number of tables to operate with.
        num_coefs_ini: First value of the range.
        num_coefs_fin: Final value of the range.
        step: Step between selected range.
        fft_func: FFT function.
    Returns list with mean times of FFT function.
    """
    times = []
    for coeficientes in np.arange(num_coefs_ini, num_coefs_fin, step):
        count = 0
        for _ in range(n_tablas):
            poly = rand_polinomio(coeficientes)
            initial_time = time.time()
            fft(poly)
            final_time = time.time()
            count += final_time - initial_time
        
        times.append(count / n_tablas)

    return times

def time_mult_polinomios_fft(n_pairs, 
                        num_coefs_ini, 
                        num_coefs_fin, 
                        step, 
                        fft_func=fft):
    """
    Generation of n_pairs integer tables and mean calculation.

    Parameters
    ----------
        n_pairs: Number of integer pairs to operate with.
        num_coefs_ini: First value of the range.
        num_coefs_fin: Final value of the range.
        step: Step between selected range.
        fft_func: FFT function.
    Returns list with mean times of FFT function.
    """
    times = []
    for coeficientes in np.arange(num_coefs_ini, num_coefs_fin, step):
        count = 0
        for _ in range(n_pairs):
            poly1 = rand_polinomio(coeficientes)
            poly2 = rand_polinomio(coeficientes)
            fft_poly1 = fft(poly1)
            fft_poly2 = fft(poly2)
            initial_time = time.time()
            mult_polinomios_fft(fft_poly1, fft_poly2)
            final_time = time.time()
            count += final_time - initial_time
        
        times.append(count / n_pairs)

    return times

def time_mult_numeros_fft(n_enteros, 
                        num_coefs_ini, 
                        num_coefs_fin, 
                        step, 
                        fft_func=fft):
    """
    Generation of n_enteros integer tables and mean calculation.

    Parameters
    ----------
        n_enteros: Number of tables to operate with.
        num_coefs_ini: First value of the range.
        num_coefs_fin: Final value of the range.
        step: Step between selected range.
        fft_func: FFT function.
    Returns list with mean times of FFT function.
    """
    times = []
    for coeficientes in np.arange(num_coefs_ini, num_coefs_fin, step):
        count = 0
        for _ in range(n_enteros):
            poly1 = rand_polinomio(coeficientes)
            poly2 = rand_polinomio(coeficientes)
            num1 = poli_2_num(poly1)
            num2 = poli_2_num(poly2)
            initial_time = time.time()
            mult_numeros_fft(num1, num2)
            final_time = time.time()
            count += final_time - initial_time
        
        times.append(count / n_enteros)

    return times

"""
AUXILIARY
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


def dg_2_ma(g):
    """
    Computes the adjacency matrix of a given graph.

    Parameters
    ----------
        g: Given graph
    Returns Adjacency matrix of the graph
    """

    n = len(g)  # Get |g|

    adjmat = np.full((n, n), np.inf)  # Prepare matrix with NumPy

    for node in np.arange(n):  # Iterates node index
        adjmat[node][node] = 0  # Sets the value to same node to 0.
        # Iterates through destinations
        for destination, cost in g[node].items():
            # Updates the matrix with destination cost.
            adjmat[node][destination] = cost[0]

    return adjmat


def floyd_warshall(ma_g):
    """
    Computes Floyd-Warshall algorithm for a given adjacency matrix of a graph.

    Parameters
    ----------
        ma_g: Adjacency Matrix
    Returns Matrix with minimum costs.
    """

    result = np.copy(ma_g)  # Save matrix variable
    n = len(ma_g)  # Store number of rows (square matrix)
    prev = np.zeros(shape=(n,n)) # Previous path

    # Previous matrix initialization
    for i in range(n):
        prev[i][i] = i
        for j in range(n):
            if ma_g[i][j] != np.inf:
                prev[i][j] = i

    # Floyd - Warshall Computation
    for k in range(n):  # For every 1D
        for i in range(n):  # For every 2D
            for j in range(n):  # For every 3D
                # Compare the cost with other paths
                if result[i][j] > result[i][k] + result[k][j]:
                    result[i][j] = result[i][k] + result[k][j]
                    prev[i][j] = k
                result[i][j] = min(result[i][j], result[i][k] + result[k][j])

    return result, prev


def bellman_ford(u, ma_g):
    """
    Bellman-Ford algorithm implementation .

    Parameters
    ----------
        ma_g: Numpy Adjacency matrix.
    Returns lists with distances from u to other and list with previous nodes.
    """
    n = len(ma_g)
    adjlist = np.copy(ma_g)
    dist = np.full(n, np.inf)
    prev = np.full(n, None)
    dist[u] = 0

    c = 0

    while c <= n:
        for i in range(n):
            for j in range(n):
                weight = adjlist[i][j]
                if weight != np.inf:
                    if dist[i] + weight < dist[j]:
                        dist[j] = dist[i] + weight
                        prev[j] = i
        c += 1
    return dist, prev


def max_length_common_subsequence(str_1, str_2):
    """
    Calculation of maximum common partial sequences Matrix.

    Parameters
    ----------
        str_1: First string as input.
        str_2: Second string as input.
    Returns matrix with lengths of maximum common partial sequences.
    """
    m = len(str_1)
    n = len(str_2)
    L = np.zeros(shape=(m+1, n+1))
  
    # Following steps build L[m+1][n+1] in bottom up fashion. Note 
    # that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1]  
    for i in range(m+1): 
        for j in range(n+1): 
            if i == 0 or j == 0: 
                L[i][j] = 0
            elif X[i-1] == Y[j-1]: 
                L[i][j] = L[i-1][j-1] + 1
            else: 
                L[i][j] = max(L[i-1][j], L[i][j-1])

    return L

def find_max_common_subsequence(str_1, str_2):
    """
    Search of possible maximum length common subsequence in given input.

    Parameters
    ----------
        str_1: First string as input.
        str_2: Second string as input.
    Returns possible maximum length common subsequence.
    """
    L = max_length_common_subsequence(str_1, str_2)
    m = len(str_1)
    n = len(str_2)
    # Following code is used to print LCS 
    index = int(L[m][n]) 
  
    # Create a character array to store the lcs string 
    lcs = [""] * index
    lcs[index-1] = "" 

  
    # Start from the right-most-bottom-most corner and 
    # one by one store characters in lcs[] 
    i = m 
    j = n 
    while i > 0 and j > 0: 
  
        # If current character in X[] and Y are same, then 
        # current character is part of LCS 
        if X[i-1] == Y[j-1]: 
            lcs[index-1] = X[i-1] 
            i-=1
            j-=1
            index-=1
  
        # If not same, then find the larger of two and 
        # go in the direction of larger value 
        elif L[i-1][j] > L[i][j-1]: 
            i-=1
        else: 
            j-=1
  
    print ("Longest Common Subsequence:", "".join(lcs))
"""

prueba = fft(np.array([1, 2, 1, 0]))

salida = invert_fft(np.array([1, 2, 1, 0]))
print("prueba:", prueba)
print("salida:", salida)
randompoly = rand_polinomio(4)
print("random:", randompoly)

hola = poli_2_num(randompoly)
print(hola)

print(time_mult_polinomios_fft(10, 5, 15, 5))

# Longest subsequence
X = "biscuit"
Y = "suitcase"
print("Lengths (i,j) of LCS are ")
print (max_length_common_subsequence(X, Y)) 

# Driver program 
X = "AGGTAB"
Y = "GXTXAYB"
find_max_common_subsequence(X, Y) 
"""
"""
var = rand_weighted_multigraph(10)
mg = dg_2_ma(var)
mg2, mg3 = floyd_warshall(mg)
print(mg)
print(mg2)
print(mg3)
"""
var = rand_weighted_multigraph(5)
mg = dg_2_ma(var)
dist, prev = bellman_ford(0,mg)
print(dist)
print(prev)