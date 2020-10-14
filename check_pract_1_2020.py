#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import argparse
import textwrap

from sklearn.linear_model import LinearRegression

os.system("cp grafos09.py grafos.py")
import grafos as gr
import grafos09 as my_gr


def fit_plot(l, func_2_fit, size_ini, size_fin, step):
    l_func_values = [func_2_fit(i) for i in range(size_ini, size_fin + 1, step)]

    lr_m = LinearRegression()
    X = np.array(l_func_values).reshape(len(l_func_values), -1)
    lr_m.fit(X, l)
    y_pred = lr_m.predict(X)

    # plt.plot(l, '*', y_pred, '-')
    return y_pred


def n2_log_n(n):
    return n ** 2. * np.log(n)


def n3_log_n(n):
    return n ** 3. * np.log(n)


def n3(n):
    return n ** 3.


# l_values =[i*n2_log_n(i) *(1 + 0.25* np.random.rand()) for i in range(10, 500+1, 10)]
# fit_plot(l_values, n2_log_n, 10, 500, 10)


####################################### main
def main(n_graphs,
         n_nodes_ini,
         n_nodes_fin,
         step,
         probability,
         num_max_multiple_edges,
         max_weight):
    """Prueba las funciones de qs.py.
    
    Args: n_graphs, n_nodes_ini, n_nodes_fin, step, sparse_factor
    """

    # check print y conversión de matriz de adyacencia
    print("\ncomprobamos funciones básicas en un grafo predefinido ..........")
    mg = {
        0: {1: {0: 10}, 2: {0: 1}},
        1: {2: {0: 1}},
        2: {3: {0: 1}},
        3: {1: {0: 1}}
    }

    gr.print_adj_list_mg(mg)

    g = {
        0: {1: {0: 10}},
        1: {2: {0: 1}},
        2: {3: {0: 1}},
        3: {1: {0: 1}}
    }

    ma_g = gr.dg_2_ma(g)
    print("\nmatriz de adyacencia del grafo\n")
    for i in range(ma_g.shape[0]):
        print(ma_g[i, :])

    _ = input("\npulsar Intro para continuar ....................\n")

    # check generar/imprimir grafo, multigrafo
    print("\ncomprobamos la generación de multigrafos dirigidos aleatorios ..........")
    r_mg = gr.rand_weighted_multigraph(n_nodes=5,
                                       probability=probability,
                                       num_max_multiple_edges=num_max_multiple_edges,
                                       max_weight=max_weight,
                                       decimals=0,
                                       fl_unweighted=False,
                                       fl_diag=True)

    print("\nlista de adyacencia del grafo generado\n")
    gr.print_adj_list_mg(r_mg)

    print(
        "\ncomprobamos la generación de grafos estándar dirigidos aleatorios y obtención de matriz de adyacencia ..........")
    r_mg = gr.rand_weighted_multigraph(n_nodes=5,
                                       probability=probability,
                                       num_max_multiple_edges=1,
                                       max_weight=max_weight,
                                       decimals=0,
                                       fl_unweighted=False,
                                       fl_diag=True)

    ma_g = gr.dg_2_ma(r_mg)
    print("\nmatriz de adyacencia del grafo generado\n")
    for i in range(ma_g.shape[0]):
        print(ma_g[i, :])

    _ = input("\npulsar Intro para continuar ....................\n")

    print("\ncomprobamos la generación de grafos no dirigidos aleatorios ..........")

    r_mg = gr.rand_weighted_undirected_multigraph(n_nodes=5,
                                                  probability=probability,
                                                  num_max_multiple_edges=1,
                                                  max_weight=max_weight,
                                                  decimals=0,
                                                  fl_unweighted=False,
                                                  fl_diag=True)

    print("\nlista de adyacencia del grafo generado\n")
    gr.print_adj_list_mg(r_mg)

    ma_g = gr.dg_2_ma(r_mg)
    print("\nmatriz de adyacencia del grafo generado\n")
    for i in range(ma_g.shape[0]):
        print(ma_g[i, :])

    _ = input("\npulsar Intro para continuar ....................\n")

    # check dijkstra y caminos óptimos
    print("\nsingle source Dijkstra ....................")
    d, p = gr.dijkstra_mg(mg, 0)
    my_d, my_p = my_gr.dijkstra_mg(mg, 0)

    print("distancias", d, my_d)
    print("previos", p, my_p)

    d_path = gr.min_paths(p)
    my_d_path = my_gr.min_paths(p)

    for k in d_path.keys():
        print("paths_from ", k, d_path[k], my_d_path[k])

    _ = input("\npulsar Intro para continuar ....................\n")

    ## timing dijkstra
    # print("\ntiming dijkstra ....................")
    #
    # l_t = gr.time_dijkstra_mg(n_graphs, n_nodes_ini, n_nodes_fin, step, num_max_multiple_edges=num_max_multiple_edges, probability=probability)
    # t_pred = gr.fit_plot(l_t, n2_log_n, size_ini=n_nodes_ini, size_fin=n_nodes_fin, step=step)
    #
    ##print((np.array(l_t) - t_pred) / np.array(l_t) * 100.)
    #
    # _ = input("\npulsar Intro para continuar ....................\n")
    #

    print("\nDijkstra all pairs minimum distances ....................")
    dist_dijkstra = gr.dijkstra_all_pairs(mg)
    my_dist_dijkstra = my_gr.dijkstra_all_pairs(mg)

    print("all_dist_dijkstra\n", dist_dijkstra, '\n', my_dist_dijkstra)

    _ = input("\npulsar Intro para continuar ....................\n")
    # check floyd warshall
    print("\nFloyd-Warshall all pairs minimum distances ....................")
    ma_g = gr.dg_2_ma(mg)
    dist_fw = gr.floyd_warshall(ma_g)
    my_dist_fw = my_gr.floyd_warshall(ma_g)

    print("all_dist_fw\n", dist_fw, '\n', my_dist_fw)

    _ = input("\npulsar Intro para continuar ....................\n")

    # check tiempos djikstra/fw
    print("\ntiming all pairs dijkstra ....................")
    l_t_d = gr.time_dijkstra_mg_all_pairs(n_graphs, n_nodes_ini, n_nodes_fin, step,
                                          num_max_multiple_edges=num_max_multiple_edges, probability=probability)
    t_pred_d = fit_plot(l_t_d, n3_log_n, size_ini=n_nodes_ini, size_fin=n_nodes_fin, step=step)

    print("\ntiming Floyd-Warshall ....................")
    l_t_f = gr.time_floyd_warshall(n_graphs, n_nodes_ini, n_nodes_fin, step,
                                   num_max_multiple_edges=num_max_multiple_edges, probability=probability)
    t_pred_f = fit_plot(l_t_f, n3, size_ini=n_nodes_ini, size_fin=n_nodes_fin, step=step)

    print("\ntiempos_dijkstra_reales   ", np.array(l_t_d).round(4))
    print("tiempos_fw_reales         ", np.array(l_t_f).round(4))
    print('\n')
    print("tiempos_dijsktra_ajustados", t_pred_d.round(4))
    print("tiempos_fw_ajustados      ", t_pred_f.round(4))


###############################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
        Script para comprobar la corrección de la práctica 1.
        Recibe num_grafos, num_nodos_inicial, num_nodos_final, step, probability, num max_edges, max_weight.
        
        Use example:
        
        cd /mnt/d/Google\ Drive/g_drive_bck/Cursos/DAA/practicas/2020_2021/python/
        ./check_pract_1_2020.py -ng 5 -ni 10 -nf 50 -s 10 -p 0.75 -me 3 -mw 10
        
        """)
    )

    parser.add_argument("-ng", "--num_graphs", type=int, default=10,
                        help="num grafos a generar en cada paso; default=10")
    parser.add_argument("-ni", "--num_nodos_inicial", type=int, default=10, help="num inicial de nodos; default=10")
    parser.add_argument("-nf", "--num_nodos_final", type=int, default=20, help="num final de nodos; default=20")
    parser.add_argument("-s", "--step", type=int, default=5, help="paso; default=5")
    parser.add_argument("-p", "--probability", type=float, default=0.5, help="probabilidad de conexión; default=0.5")
    parser.add_argument("-me", "--max_edges", type=int, default=1, help="num max edges; default=1")
    parser.add_argument("-mw", "--max_weight", type=float, default=10., help="max weight; default=10")

    args = parser.parse_args()

    main(args.num_graphs,
         args.num_nodos_inicial,
         args.num_nodos_final,
         args.step,
         args.probability,
         args.max_edges,
         args.max_weight)
