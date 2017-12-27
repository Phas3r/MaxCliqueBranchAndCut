# MaxCliqueBranchAndCut
# Description

Realisation of branch and cut algorithm for finding exact solution for the Maximum clique problem using greedy graph heuristics.<br>
Project also contains Cliquer solver:  https://users.aalto.fi/~pat/cliquer.html (can be used for debugging purposes)
This is C++ realization for branch-and-cut (for example, see http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.40.9473&rep=rep1&type=pdf) algorithm for the maximum clique problem: 
In every branch we solve the corresponding LP relaxation for the maximum clique problem using Cplex
## Instructions

To run the executable you should write in the command line:
- [executable_name] [path_to_the_file] [time_limit_in_seconds]

Default timeout is 0

## Test results (for those graphs, where clique was found before the program was terminated by 3600 sec timeout)

Graph name|#Nodes|#Edges|Found clique size|Time (ms)
---|---|---|---|---
brock200_2.clq|200|9876|12|2839031
brock200_4.clq|200|13089|21|1897371
c-fat200-1.clq|200|1534|12|1120
c-fat200-2.clq|200|3235|24|552
c-fat200-5.clq|200|8473|58|1252
c-fat500-1.clq|500|4459|14|26139
c-fat500-2.clq|500|9139|26|54228
c-fat500-5.clq|500|23191|64|55274
c-fat500-10.clq|500|46627|126|2223
hamming6-2.clq|64|1824|32|5
hamming6-4.clq|64|704|4|20
hamming8-4.clq|256|20864|16|981
johnson8-2-4.clq|28|210|4|0
johnson8-4-4.clq|70|1855|14|18
johnson16-2-4.clq|120|5460|8|33
MANN_a9.clq|45|918|16|20
MANN_a27.clq|378|70551|126|207761
p_hat300_1.clq|300|10933|8|703782
p_hat300_2.clq|300|21928|25|2350891
p_hat500_1.clq|500|31569|9|1637892
san200_0.7_1.clq|200|13930|30|200
san200_0.7_2.clq|200|13930|18|8526
san200_0.9_1.clq|200|17910|70|120
san200_0.9_2.clq|200|17910|60|65
san200_0.9_3.clq|200|17910|44|870755
san400_0.5_1.clq|400|39900|13|161552
C125.9.clq|125|6963|34|733421
gen200_p0.9_55.clq.txt|200|17910|55|2480
keller4.clq|171|9435|11|607662
