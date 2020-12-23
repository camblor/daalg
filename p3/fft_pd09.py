import numpy as np
import time
import math


def get_minimum_power(n):
    """
    Calculation of minimum power of 2 greater than n.

    Parameters
    ----------
        n: Length to get minimum power of 2.
    Returns Minimum power of 2 greater than n.
    """
    return 2 ** math.ceil(math.log(n, 2))


def append_zeros(t):
    """
    Appending of 0 values to table.

    Parameters
    ----------
        t: Table to append 0's to.
    Returns table with 0's appended until len is minpow.
    """
    n = len(t)
    minpow = get_minimum_power(n)

    fft_input = np.array(t)
    zeros = np.array([0] * (minpow-n))

    return np.append(fft_input, zeros)


def recursive_fft(t):
    """
    Calculation of FFT from a NumPy array. Recursive implementation.

    Parameters
    ----------
        t: Table.
    Returns Fast Fourier Transform of input.
    """
    # FFT RETURN VARIABLE
    result = np.array([])

    # Get table length (Is power of 2)
    length = len(t)
    prev2k = math.ceil(length / 2)

    # Base case with unique value
    if length <= 1:
        return t

    # Separate into index even and odds
    even, odds = t[::2], t[1::2]

    # Recursive call
    f_e = recursive_fft(even)
    f_o = recursive_fft(odds)

    for i in range(length):
        # First Transformate Operand
        first = f_e[i % prev2k]

        # Second Transformate Operand
        tmp1 = (2 * np.pi * 1j) * (i / length)
        tmp2 = f_o[i % prev2k]
        second = np.exp(tmp1) * tmp2

        # Transformate result storage
        result = np.append(result, first + second)

    return result


def fft(t):
    """
    Calculation of FFT from a NumPy array.

    Parameters
    ----------
        t: Table.
    Returns Fast Fourier Transform of input.
    """

    # Minimum power of 2 + 0's appending
    fft_input = append_zeros(t)

    # Call to recursive function∫
    fft_output = recursive_fft(fft_input)

    # Output returnal
    return fft_output


def invert_fft(t, fft_func=fft):
    """
    Application of inversion algorithm of DFT.

    Parameters
    ----------
        t: Table.
        fft_func: FFT implementation function.
    Returns random protein sequence.
    """
    k = len(t)

    # Conjugate input
    conjugate = np.conj(t)

    # Compute DFT
    transformate = fft_func(conjugate)

    # Invert DFT
    output = [np.conj(n) / k for n in transformate]

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
    return np.random.randint(base, size=long).astype(type(int()))


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
    ret = int(0)
    for i in range(len(l_pol)):
        ret += int(np.real(l_pol[i]*pow(base, i)))
    return ret


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
    poli = np.array([])
    while num > 0:
        num, poli_part = divmod(num, base)

        poli = np.append(poli, poli_part)

    return poli.astype(type(int()))


def padding_mult(l_pol_1, l_pol_2):
    """
    0s padding to generate the fft multiplication array.

    Parameters
    ----------
        l_pol_1: Array 1.
        l_pol_2: Array 2.
    Returns both arrays with appended 0's.
    """
    # Get final length
    length = len(l_pol_1)+len(l_pol_2)-1

    # Get power of 2
    minpow = get_minimum_power(length)

    # Append to first array
    nzeros = minpow - len(l_pol_1)
    zerosarr = np.array([0]*nzeros)
    p1 = np.concatenate((l_pol_1, zerosarr), axis=0)

    # Append to first array
    nzeros = minpow - len(l_pol_2)
    zerosarr = np.array([0]*nzeros)
    p2 = np.concatenate((l_pol_2, zerosarr), axis=0)

    return p1, p2


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
    poli1, poli2 = padding_mult(l_pol_1, l_pol_2)

    # Fast Fourier transformates
    coefficient1 = fft_func(poli1)
    coefficient2 = fft_func(poli2)

    # Coefficient multiplication (SECOND STEP)
    res = [i*j for i, j in zip(coefficient1, coefficient2)]

    # Inversion algorithm with FFT (FINAL STEP)
    output = invert_fft(res, fft_func=fft_func)
    output = [int(np.real(i)) for i in np.rint(output)]

    # Round and return NumPy Array
    return output


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
    # Transform to polynomials
    poli_num1 = num_2_poli(num1)
    poli_num2 = num_2_poli(num2)

    # Reset to numbers after multiplication
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
    prev = np.full((n, n), None)  # Previous path
    

    # Previous array initialization
    for i in np.arange(n):
        prev[i][i] = i
        for j in np.arange(n):
            if result[i][j] != np.inf:
                prev[i][j] = i

    for k in range(n):  # For every 1D
        for i in range(n):  # For every 2D
            for j in range(n):  # For every 3D
                # Compare the cost with other paths
                if result[i][j] > result[i][k] + result[k][j]:
                    result[i][j] = result[i][k] + result[k][j]
                    # Previous update
                    prev[i][j] = prev[k][j]

    return result, prev


def bellman_ford(u, ma_g):
    """
    Bellman-Ford algorithm implementation .

    Parameters
    ----------
        ma_g: Numpy Adjacency matrix.
    Returns lists with distances from u to other and list with previous nodes.
    """
    # Get length
    length = len(ma_g)

    # Initialize arrays
    dist = np.array([np.inf for _ in range(length)])
    dist[u] = 0
    prev = np.array([None for i in range(length)])
    prev[u] = u

    # Algorithm implementation (k-1 times)
    for k in range(length-1):
        for i in range(length):
            for j in range(length):
                # There is an edge u -> i
                if ma_g[i][j] != np.inf:
                    found = dist[i] + ma_g[i][j]
                    # Update procedure
                    if dist[j] > found:
                        dist[j] = found
                        prev[j] = i

    # Output returnal
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

    # Counting maximum sequence length
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif str_1[i-1] == str_2[j-1]:
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
    index = int(L[m][n])

    # Create a character array to store the lcs string
    lcs = [""] * index
    lcs[index-1] = ""

    # Start from the right-bottom-most corner and
    # store characters in lcs[]
    i = m
    j = n
    while i > 0 and j > 0:
        # If current character str_1[] and str_2 are same
        # character is part of LCS
        if str_1[i-1] == str_2[j-1]:
            lcs[index-1] = str_1[i-1]
            i -= 1
            j -= 1
            index -= 1

        # If not same, find the larger and
        # go in the direction of larger
        elif L[i-1][j] > L[i][j-1]:
            i -= 1
        else:
            j -= 1

    return "".join(lcs)
