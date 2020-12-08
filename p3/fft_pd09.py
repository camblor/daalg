import numpy as np

#TODO devolver tabla con 0s añadidos y calcular mediante logaritmo
def get_minimum_power(n):
    """
    Calculation of minimum power of 2 greater than n.

    Parameters
    ----------
        n: input number to calculate with.
    Returns Minimum power of 2 greater than n.
    """
    minpow = 1
    while minpow < n:
        minpow *= 2
    return minpow

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
    minpow = get_minimum_power(length)

    # Power of 2 case
    fft_input = t

    # Not power of 2 case: Appending 0's
    if minpow > length:
        for i in range(minpow-length):
            fft_input = np.append(fft_input, [0])
        
    # Base case with unique value
    if length <= 1: 
        return fft_input
    
    # Get odd and even elements
    for element, i in zip(t, range(length)):
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
    # Degree calculation
    degree = len(l_pol_1) + len(l_pol_2) - 1

    # Adjust of FFT input TODO
    poli1 = _len2k(l_pol_1, degree)
    poli2 = _len2k(l_pol_2, degree)

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
    print(n_tablas)

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
    print(n_pairs)

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
    print(n_enteros)

def floyd_warshall(ma_g):
    """
    Floyd-Warshall algorithm implementation .

    Parameters
    ----------
        ma_g: Numpy Adjacency matrix.
    Returns lists with distances and list with previous nodes.
    """
    print(ma_g)

def bellman_ford(ma_g):
    """
    Bellman-Ford algorithm implementation .

    Parameters
    ----------
        ma_g: Numpy Adjacency matrix.
    Returns lists with distances from u to other and list with previous nodes.
    """
    print(ma_g)

def max_length_common_subsequence(str_1, str_2):
    """
    Calculation of maximum common partial sequences Matrix.

    Parameters
    ----------
        str_1: First string as input.
        str_2: Second string as input.
    Returns matrix with lengths of maximum common partial sequences.
    """
    print(str_1, str_2)

def find_max_common_subsequence(str_1, str_2):
    """
    Search of possible maximum length common subsequence in given input.

    Parameters
    ----------
        str_1: First string as input.
        str_2: Second string as input.
    Returns possible maximum length common subsequence.
    """
    print(str_1, str_2)


prueba = fft(np.array([1, 2, 1, 0]))

salida = invert_fft(np.array([1, 2, 1, 0]))
print("prueba:", prueba)
print("salida:", salida)
randompoly = rand_polinomio(4)
print("random:", randompoly)

hola = poli_2_num(randompoly)
print(hola)

