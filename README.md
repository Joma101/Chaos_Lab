# Chaos_Lab

README for Quotient_Network_Expansion_Algorithm.m

INPUT:

When inputting the quotient network tuple into the algorithm, only the removal of separators such as commas and parentheses in the set notation is necessary, averting the need to re-organize any values. See the example input below.

Given the following quotient network tuple example, R:

V = {1,2,3}
E = {(2,1),(3,2)}
A = [0 0 0; 1 2 0; 0 1 0]

The input of R to the algorithm can be written as the following (Qadj is the adjacency matrix of the quotient network, A):

V = 1:3;
E = [2 1 3 2];
Qadj = [0 0 0; 1 2 0; 0 1 0];

VARIABLE N_USER:

The variable "n_user" exists to enable user input for the array "n" for Case II network expansions. When desiring a Case I network expansion, you must write "n_user = 0;". For Case II network expansions, n_user must be the same number of elements as V above, and will instruct the script to set n equal to n_user. For example, "n_user = [10 10 10];" will instruct the script to generate a 30-node expanded network from the quotient, in which there are 10 nodes per cluster. The first element in n_user corresponds to the first element in the input V, and so on.

FEATURES:

Although the inputs follow standard graph/network notation (V,E,A), our proposed algorithm does not run on this language. It runs on arrays n, s, c, and f translated seamlessly from the tuple, R, then translates back into a node and edge list output along with the expanded network adjacency matrix A.

The script will save the node and edge lists of the expanded network produced toward the end in the working directory, overwriting any previously saved node/edge lists with the same filename.

(Within the script workspace viewed using the MATLAB application) The in-degree of each expanded network node is reported in the variable "in_degree". The index of each element in in_degree same as the index of the corresponding node in the expanded network.

The matrix "Gephi_Colors" is an RGB color matrix that can be used to distinguish nodes in different clusters by color, particularly in the Gephi network graphing software. Each row in the matrix corresponds to each element in the input V.

If a cluster is isolated from its neighbors, the script will ensure there still exists at least one node in the expanded network for that isolated cluster.