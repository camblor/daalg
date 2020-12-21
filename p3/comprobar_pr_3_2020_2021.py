#!/usr/bin/env python
# -*- coding: utf-8 -*-
import string, random
import numpy as np
import copy
import sys
import os

import argparse

eps = 1.e-5

############################################################################## funciones auxiliares
def dg_2_ma(mg):    
    """Returns the adj matrix of a multigraph.
    It only considers the firat edge.
    """
    n_v = len(mg.keys())
    ma_g = np.inf * np.ones( (n_v, n_v) )
    
    for u in sorted( mg.keys() ): #las claves del dicc pueden estar NO ordenadas
        ma_g[u, u] = 0.
    
        for v in sorted(mg[u].keys()):
            if v != u:
                ma_g[u, v] = mg[u][v][0]
    
    return ma_g

    
############################################################################## main
#def check_pr_2(args):
def check_pr_3(num_iters):
    """Prueba las funciones de secuencias.py desarrolladas en la práctica 2
    """
    
    np.set_printoptions(precision=3)
            
    # #### 1. Checking fft, invff and np.fft
    
    for i in range(num_iters):
        # * Generar aleatoriamente tabla de np.complex
        pot = np.random.randint(5, 9)
        len_t = int(2.**pot)
        t_re = np.random.randn(len_t)
        t_im = np.random.randn(len_t)
        t = t_re + t_im * np.complex(0, 1)
        
        # * Comprobar que fp.fft(t) - np.fp.fft(t) = 0
        #our fft
        fft_t = fp.fft(t)
        
        #Numpy's fft
        np_fft_t = np.conj(np.fft.fft(np.conj(t)))
        
        if np.linalg.norm(fft_t- np_fft_t) > eps:
            print("iter", i)
            print("error en fft:\n", fft_t, "\n", np_fft_t)
            break

        # * Comprobar que invert_fft(fp.fft(t) = t
        inv_fft_t = fp.invert_fft( fp.fft(t) )
        
        if np.linalg.norm(t - inv_fft_t[: len(t)]) > eps:
            print("iter", i)
            print("error en fft:\n", t, "\n", inv_fft_t)
            break

    print("\nvalues of the first five elements of last generated t, fft and inv_fft")
    print("\tt and inv_fft must coincide up to numerical precision")
    print("t      ", t[: 5])
    print("fft    ", fft_t[: 5])
    print("inv_fft", inv_fft_t[: 5])

    print("ok")

    _ = input("pulsar Intro para continuar ....................\n")

    # #### 2. Checking mult_polis estandar y fft, invff and np.fft
    print("\nchecking mult de polinomios ..........")

    for i in range(num_iters):
        # * Generar aleatoriamente polinomios
        pot  = np.random.randint(5, 9)
        base = np.random.randint(2, 11)

        p1 = fp.rand_polinomio(int(2**pot), base=base)
        p2 = fp.rand_polinomio(int(2**pot), base=base)

        # * Comprobar que mult fft y la mult de numpy coinciden
        prod_f = fp.mult_polinomios_fft(p1, p2)
        prod_np = np.polymul(p1[ : : -1], p2[ : : -1])
        
        if np.linalg.norm(prod_np[ : : -1] - np.array(prod_f[ : len(prod_np)])) > 0.:
            print("iter", i)
            print("error en prod:\n", prod_np, "\n", prod_ff[ : len(prod_np)])
            break

    print("first coefficients of the last generated p1, p2, prod_fft, prod_numpy")
    print("p1      ", p1[: 5])
    print("p2      ", p2[: 5])
    print("prod_fft", prod_f[: 5])
    print("prod_np ", list(prod_np[ : : -1])[ : 5])

    print("ok")
    _ = input("pulsar Intro para continuar ....................\n")

    # #### 3. Checking mult_numeros estandar, fft y python
    #
    # * Generar aleatoriamente polinomios
    # * Comprobar que mult fft, estandar y python coinciden
    print("\nchecking mult de números ..........")

    for i in range(num_iters):
        num_d = np.random.randint(50, 101)
        base = np.random.randint(4, 11)

        num1 = fp.rand_numero(num_d, base=base)
        num2 = fp.rand_numero(num_d, base=base)
        
        prod_f = fp.mult_numeros_fft(num1, num2)
        prod_p = num1 * num2

        if prod_p != prod_f:
            print("iter", i)
            print("error en prod s o f:\n", prod_p, "\n", prod_f)
            break

    print("values of the last generated num1, num2, prod_python y prod_fft")
    print("base", base, "num_digits", num_d)
    print("num1\t", num1)
    print("num2\t", num2)
    print("prod_python\t", prod_p)
    print("prod_fft\t", prod_f)

    print("ok")    
    _ = input("pulsar Intro para continuar ....................\n")
    
    # #### 4. Floyd-Warshall y Bellman-Ford
    print("checking floyd-warshall ...")
    g1 = {
          0: {1: {0: 10}, 2: {0: 1}}, 
          1: {2: {0: 1}}, 
          2: {3: {0: 1}},
          3: {1: {0: 1}}
         }
        
    g2 = {
         0: {1: {0:2}, 3: {0:1}},
         1: {2: {0:3}},
         2: {3: {0:-5}},
         3: {}
        }
  
    g3 = {
         0: {1: {0:2}, 4: {0:3}},
         1: {2: {0:2}},
		 2: {},
         3: {2: {0:3}},
         4: {3: {0:-3}},
        }

    for i, g in enumerate([g1, g2, g3]):
        print("checking in graph ", i)
        ma_g = dg_2_ma(g)
        d, p = fp.floyd_warshall(ma_g)
        print(d)
        print(p)
        
        _ = input("pulsar Intro para continuar ....................\n")
        print("checking bellman_ford ...")
        
        ma_g = dg_2_ma(g)
        for u in range(ma_g.shape[0]):
            d, p = fp.bellman_ford(u, ma_g)
            print("distances from ", u)
            print(d)
            print(p)
    
        _ = input("pulsar Intro para continuar ....................\n")

    # #### 5. Secuencias
    print("checking sequences ...")
    
    str_1 = "algoritmos"
    str_2 = "logaritmos"
    
    print(str_1, str_2)
    mlcs_matrix = fp.max_length_common_subsequence(str_1, str_2)
    print("lcs_matrix")
    print(mlcs_matrix)
    
    print("max_common_subseq\n", fp.find_max_common_subsequence(str_1, str_2))
    
##########################################################################################################        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="comprobador de la práctica 2")
    
    parser.add_argument("-p", "--pareja", type=str, default=None)
    parser.add_argument("-ni", "--num_iters", type=int, default=5)
    
    args = parser.parse_args()
    
    if args.pareja is not None:
        f_path = "." + "/fft_pd" + args.pareja + ".py"
        if os.path.isfile(f_path):
            #str_comm = "cp ./p3" + args.pareja + "/fft_pd" + args.pareja + ".py  ./fft_pd.py"
            #str_comm = "cp euler" + args.pareja + ".py  ./euler.py"
            #print(str_comm)
            #os.system(str_comm)
            #importar aquí vuestro fft_pdXX.py
            import fft_pd09 as fp
            
            _ = input("pulsar Intro para continuar ....................\n")

            check_pr_3(args.num_iters)
        
        else:
            print("no file " + f_path)
    else:
        print("falta num pareja")