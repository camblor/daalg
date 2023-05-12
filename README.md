# daalg
Asignatura Diseño y Análisis de Algoritmos - Universidad Autónoma de Madrid 2020/2021

## Practica 1

Graph theory and algorithms. Dijkstra and Floyd-Warshall. Functions to get real speed comparisons.

Graph Theory: Graph theory is a field of mathematics and computer science that studies graphs, which are mathematical structures used to model pairwise relations between objects. A graph is made up of vertices (also called nodes or points) which are connected by edges (also called arcs or lines).

Graphs can be used to represent many real-world systems, such as social networks (people are vertices, friendships are edges), web pages (web pages are vertices, hyperlinks are edges), biological networks (species are vertices, predator-prey relationships are edges), and many others.

Graph Algorithms: Graph algorithms are a set of procedures that can be used to solve common problems on graphs. These problems include finding the shortest path between two nodes, determining whether a graph is connected, finding a maximum matching, and many others. These algorithms are fundamental in many fields, including computer networking, transportation, and operations research.

Dijkstra's Algorithm: This is an algorithm developed by computer scientist Edsger Dijkstra in 1956. It is used to find the shortest path between nodes in a graph, which may represent, for example, road networks. It operates by building a set of nodes that have minimum distance from the starting node and visiting neighbors of these nodes to find the overall shortest path.

Floyd-Warshall Algorithm: The Floyd-Warshall algorithm is another important graph algorithm. It solves the all-pairs shortest path problem, which involves finding the shortest path between all pairs of nodes in a weighted graph (a graph where each edge has a "length" or "cost"). This is in contrast to Dijkstra's algorithm, which only solves the single-source shortest path problem.

The Floyd-Warshall algorithm works by repeatedly improving an estimate on the shortest path between two vertices, until ultimately, the estimate is the shortest path. The key idea is to consider each node and update all shortest paths that could be improved by going through that node.

The Floyd-Warshall algorithm is particularly useful when dealing with dense graphs, where most or all pairs of nodes are connected by edges, or when all shortest paths are required, not just a single shortest path from one node to another.

## Practica 2

Graph theory and Eulerian paths. Protein sequencing.

An Eulerian path in a graph is a path that traverses every edge exactly once. This concept comes from graph theory, a branch of mathematics. The path is named after the Swiss mathematician Leonhard Euler, who first proposed the problem of the Seven Bridges of Königsberg, a historical problem that led to the development of graph theory.

The problem of the Seven Bridges of Königsberg asked if it was possible to walk through the city of Königsberg (now Kaliningrad, Russia) crossing each of its seven bridges exactly once. Euler proved that such a path did not exist, thus laying the groundwork for the theory of graphs.

Eulerian paths have an interesting connection to the field of bioinformatics, particularly in the problem of protein sequencing and DNA sequencing. The process of reconstructing the original sequence of a DNA or protein from smaller fragments is often likened to an Eulerian path problem.

In the classic Eulerian path problem, each "edge" of the graph corresponds to a fragment of DNA or protein, and we want to find the sequence (path) that includes each fragment exactly once. This analogy isn't perfect because in biological sequencing problems, errors and repeats can occur, making the problem more complex than a pure Eulerian path problem.

For example, De Bruijn graphs, a specific type of directed graph, have been used in modern sequencing techniques to solve the problem of assembling short DNA or protein sequences into a complete genome or protein sequence. In these graphs, each vertex represents a sequence of length k (a k-mer), and directed edges connect vertices with overlapping sequences. A Eulerian path in this graph corresponds to a sequence that uses every k-mer exactly once, which could be a plausible reconstruction of the original sequence.

However, it's important to note that in real-world sequencing, various complications arise, such as errors in sequencing, repeats in the genome, and the fact that sequencing reads come from both strands of the DNA double helix. Therefore, sophisticated algorithms and software tools have been developed to deal with these complications in practice.

## p3

Fast Fourier Transform in order to perform fast addition. Floyd-Warshall algorithm, Bellman-Ford algorithm, maximum length common subsequence.

Fast Fourier Transform (FFT): The Fast Fourier Transform is an algorithm used to compute the Discrete Fourier Transform (DFT) and its inverse. The DFT is a mathematical technique used in signal processing to transform a sequence of complex numbers into another sequence of complex numbers, providing a frequency domain representation of the original sequence. FFT reduces the computational complexity of calculating DFT from O(n^2) to O(n log n), where n is the size of the data set.

In the context of fast addition or multiplication, FFT is typically used in multiplying large integers and polynomials. This is achieved by representing the numbers or polynomials as points on the complex plane, performing the FFT, multiplying each corresponding point together, and then performing the Inverse FFT to get back the result in the time domain. This approach is much faster than traditional long multiplication for very large numbers or polynomials.

Floyd-Warshall Algorithm: This is a dynamic programming algorithm used to find the shortest paths between all pairs of vertices in a weighted graph with positive or negative edge weights but no negative cycles. The algorithm works by incrementally improving an estimate of the shortest path between every pair of vertices until the estimate is optimal.

Bellman-Ford Algorithm: The Bellman-Ford algorithm is another dynamic programming algorithm used for finding the shortest paths from a single source vertex to all other vertices in a weighted digraph. It handles both positive and negative edge weights and is capable of detecting negative cycles, which is something the Dijkstra’s algorithm can't do.

Maximum Length Common Subsequence (MLCS): The MLCS problem is a classic computer science problem that involves finding a maximum-length sequence of characters that appear left-to-right (but not necessarily in a contiguous block) in each of a given set of sequences. This problem arises in many contexts, such as bioinformatics, where it can be used to align DNA or protein sequences to identify similar regions. The MLCS problem can be solved using dynamic programming, where we construct a table that contains the lengths of longest common subsequences for prefixes of the input sequences, and then use this table to construct the MLCS itself.
