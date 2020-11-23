#!/usr/bin/env python
# -*- coding: utf-8 -*-
import string, random
import numpy as np
import copy
import sys
import os

import argparse

from sklearn.linear_model import LinearRegression

import euler_secuencias_2020 as eu20

############################################################################## funciones auxiliares
    
############################################################################## main
#def check_pr_2(args):
def check_pr_2():
    """Prueba las funciones de secuencias.py desarrolladas en la práctica 2
    """
    
    np.set_printoptions(precision=3)
    
    ######################### guardar grafos como TFG
    print("\ncomprobamos las funciones de guardar grafos como TFG ..........")
    
    d_g = { 0: {0: {0:1}, 3: {0:1}}, 
            1: {2: {0:1}}, 
            2: {0: {0:1}, 3: {0:1}}, 
            3: {1: {0:1}, 2: {0:1}}
           }
    
    

    print("\nguardamos como TFG y mostramos archivo..........\n")
    f_name = 'my_graph.tfg'
    eu20.d_g_2_TGF(d_g, f_name)
    os.system("cat %s" % f_name)
    
    print("\n")
    f_name = 'my_graph_19.tfg'
    eu20.d_g_2_TGF(d_g, f_name)
    os.system("cat %s" % f_name)
    
    _ = input("\npulsar Intro para continuar ....................\n")
    print("\nreleemos TFG y comprobamos ..........")
    d_g2 = eu20.TGF_2_d_g(f_name)
    eu20.print_adj_list_mg(d_g2)
    
    print("\n")
    d_g2_19 = eu20.TGF_2_d_g(f_name)
    eu20.print_adj_list_mg(d_g2_19)
    
    _ = input("\npulsar Intro para continuar ....................\n")
    

    #print("....................................................................................................")
    print(".................... checking euler circuit and sequence reconstruction ....................")
    
    ##### varios ejemplos de grafos en principio dirigidos
    # Añadir ejemplos que se hayan usado
    
    
    print("\ncomprobamos existencia de circuitos eulerianos ..........")
    #d_g = { 0: {1: {0: 1}}, 1: {3: {0: 1}}, 2: {0: {0: 1}}, 3: {4: {0:1}}, 4: {9: {0:1}}, 5: {8: {0:1}}, 6: {5: {0: 1}}, 7: {6: {0:1}}, 8: {2: {0:1}}, 9: {7: {0:1}} }
    d_g = { 0: {0: {0:1}, 3: {0:1}}, 
            1: {2: {0:1}}, 
            2: {0: {0:1}, 3: {0:1}}, 
            3: {1: {0:1}, 2: {0:1}}
           }
    
    d_g_test_circ = eu20.TGF_2_d_g("tgf_test_circ.txt")
    eu20.print_adj_list_mg(d_g_test_circ)
    print("\ncalculamos tablas adj, inc ..........")
    print("adj_inc\n", eu.adj_inc_directed_multigraph(d_g_test_circ) )
    
    if eu.isthere_euler_circuit_directed_multigraph(d_g_test_circ):
        print("\nejemplo de paseo ..........")
        d_g_copy = copy.deepcopy(d_g_test_circ)
        print( eu.euler_walk_directed_multigraph(0, d_g_copy) )
        
        d_g_copy = copy.deepcopy(d_g_test_circ)
        print( eu20.euler_walk_directed_multigraph(0, d_g_copy) )
        
        print("\nencontramos circuitos eulerianos ..........")
        d_g_copy = copy.deepcopy(d_g_test_circ)
        print( eu.euler_circuit_directed_multigraph(d_g_copy, 0) )
        
        d_g_copy = copy.deepcopy(d_g_test_circ)
        print( eu20.euler_circuit_directed_multigraph(d_g_copy, 0) )
        #print([0, 0, 3, 2, 3, 1, 2, 0])
        
    else:
        print("\nparece que no hay circuito euleriano ..........")
    
    _ = input("pulsar Intro para continuar ....................\n")
    
    
    print("\ncomprobamos existencia de caminos eulerianos ..........")
    #d_g = { 0: {1: {0: 1}}, 1: {3: {0: 1}}, 2: {0: {0: 1}}, 3: {4: {0:1}}, 4: {9: {0:1}}, 5: {8: {0:1}}, 6: {}, 7: {6: {0:1}}, 8: {2: {0:1}}, 9: {7: {0:1}} }
    d_g = { 0: {0: {0:1}, 3: {0:1}}, 
            1: {2: {0:1}}, 
            2: {0: {0:1}, 3: {0:1}}, 
            3: {1: {0:1}}
           }
    d_g_test_path = eu20.TGF_2_d_g("tgf_test_path.txt")
    
    print("\ncalculamos tablas adj, inc ..........")
    print("adj_inc\n", eu.adj_inc_directed_multigraph(d_g_test_path) )
    
    if eu.isthere_euler_path_directed_multigraph(d_g_test_path):
        print("\nvértices inicial y final ..........")
        u, v = eu.first_last_euler_path_directed_multigraph(d_g_test_path)
        print(u, v)
    
        print("\nejemplo de paseo ..........")
        d_g_copy = copy.deepcopy(d_g_test_path)
        print( eu.euler_walk_directed_multigraph(u, d_g_copy) )
        
        print("\nencontramos circuitos eulerianos ..........")
        d_g_copy = copy.deepcopy(d_g_test_path)
        print( eu.euler_path_directed_multigraph(d_g_copy) )
        
        d_g_copy = copy.deepcopy(d_g_test_path)
        print( eu20.euler_path_directed_multigraph(d_g_copy) )
        #print([2, 0, 0, 3, 1, 2, 3])
    
    else:
        print("\nparece que no hay camino ..........")
    
    _ = input("pulsar Intro para continuar ....................\n")
    
        
    print("\ncomprobamos cálculos de espectros y grafos sobre secuencias sencillas ..........")
    #seq = "AACCGGTT"
    seq = "AAACCC"
    #seq = eu.random_sequence(4)
    print(seq)
    spec = eu.spectrum(seq, len_read=3)
    print("spec", spec)
    spec_2 = eu.spectrum_2(spec)
    print("spec_2", spec_2)
    d_g = eu.spectrum_2_graph(spec)
    print("dict", d_g)
    seq_rec = eu.spectrum_2_sequence(spec)
    print("seq_rec", seq_rec)
    
    _ = input("pulsar Intro para continuar ....................\n")

    print("comprobamos reconstrucción ..........")
    # generamos secuencia aleatoria y calculamos espectro y grafo
    len_seq = 100 + np.random.randint(50)
    len_read = 3
    print("len_seq", len_seq,  "len_read", len_read, '\n')
    
    seq = eu.random_sequence(len_seq)
    spec = eu.spectrum(seq, len_read)
    d_g  = eu.spectrum_2_graph(spec)
    
    print(seq)
    
    # buscamos circuito euleriano y reconstruimos cadenas 
    e_path = eu.euler_path_directed_multigraph(d_g)
    spec_2 = eu.spectrum_2(spec)
    seq_rec = eu.path_2_sequence(e_path, spec_2)
    print('\t'+seq+'\n', '\t'+seq_rec+'\n')
    
    # comprobamos coherencia de la reconstrucción
    spec_rec = eu.spectrum(seq_rec, len_read)
    print("seq  iguales?", seq == seq_rec)
    print("spec iguales?", sorted(list(spec)) == sorted(list(spec_rec)))
    
    _ = input("pulsar Intro para continuar ....................\n")
    
    for _ in range(10):
        len_seq = 200 + np.random.randint(100)
        len_read = np.random.randint(5, 10)
        print("len_seq", len_seq,  "len_read", len_read, '\n')
    
        seq = eu.random_sequence(len_seq)    
        spec = eu.spectrum(seq, len_read)
        d_g  = eu.spectrum_2_graph(spec)
        
        # buscamos circuito euleriano y reconstruimos cadenas 
        e_path = eu.euler_path_directed_multigraph(d_g)
        spec_2 = eu.spectrum_2(spec)
        
        seq_rec = eu.path_2_sequence(e_path, spec_2)
        seq_rec_2 = eu20.path_2_sequence(e_path, spec_2)
        
        print('\t'+seq[ : 50]+'\n', '\t'+seq_rec[ : 50]+'\n', '\t'+seq_rec_2[ : 50]+'\n')
        
        # comprobamos coherencia de la reconstrucción
        spec_rec   = eu.spectrum(seq_rec, len_read)
        spec_rec_2 = eu20.spectrum(seq_rec_2, len_read)
        
        print("seq  iguales?", seq == seq_rec)
        print("spec iguales?", sorted(list(spec)) == sorted(list(spec_rec)))
        if sorted(list(spec)) != sorted(list(spec_rec)) or sorted(list(spec)) != sorted(list(spec_rec_2)) :
            print(len_seq, len_read)
            print(seq)
            break
    
    
##########################################################################################################        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="comprobador de la práctica 2")
    
    parser.add_argument("-p", "--pareja", type=str, default=None)
    
    args = parser.parse_args()
    
    if args.pareja is not None:
        #f_path = "./p2" + args.pareja + "/euler" + args.pareja + ".py"
        f_path = "euler" + args.pareja + ".py"
        if os.path.isfile(f_path):
            #str_comm = "cp ./p2" + args.pareja + "/euler" + args.pareja + ".py  ./euler.py"
            str_comm = "cp euler" + args.pareja + ".py  ./euler.py"
            print(str_comm)
            os.system(str_comm)
            import euler as eu
        
            check_pr_2()
        
        else:
            print("no file " + f_path)
    else:
        print("falta num pareja")