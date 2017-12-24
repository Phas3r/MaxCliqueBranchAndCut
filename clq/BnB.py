
# coding: utf-8

# In[15]:


import cplex
import numpy as np
import sys
import math
import networkx as nx
import os


max_clique = []

def parser(filename):
    G = nx.Graph()   
    asd = os.getcwd()    
    with open(filename) as f:
        for line in f:
            if line.startswith("e"):
                splitted_line = line.split(" ")
                G.add_edge(int(splitted_line[1]), int(splitted_line[2]))   
    return(G)

def independent_sets(G):    
    CG = nx.coloring.greedy_color(G, strategy=nx.coloring.strategy_independent_set)
    print(type(CG))
    independent_sets = [[] for i in range(max(CG.values())+1)]
    color = []
    for i in range(1,G.number_of_nodes()+1):
        print(type(CG[i]))
        if CG[i] not in color:
            color.append(CG[i])        
        independent_sets[CG[i]].append(i)
    return independent_sets

def independent_sets2(G):
    ind_sets = []
    strategies = [nx.coloring.strategy_independent_set]
    for strategy in strategies:
        d = nx.coloring.greedy_color(G, strategy=strategy)
        
        for color in set(color for node, color in d.items()):
            ind_sets.append(
                [key for key, value in d.items() if value == color])
    return ind_sets


def candidates(node, G):   
    return list(nx.complement(G)[node].keys())

def cand_in_set(candidates, ind_set):
    for i in range(len(candidates)):
        if i in ind_set:
            candidates.remove(i)
    return candidates

def base_model_solution(independent_sets, vertex_number, G):
    problem = cplex.Cplex()
    problem.objective.set_sense(problem.objective.sense.maximize)
    types = 'C' * vertex_number    
    problem.variables.add(obj=[1.0] * vertex_number,
                          ub=[1.0] * vertex_number,
                          names=['x{0}'.format(x) for x in range(1,vertex_number+1)], types = types)
    constraints = []
    not_connected = nx.complement(G).edges()
    for ind_set in independent_sets: 
        constraints.append([['x{0}'.format(i) for i in ind_set], [1.0] * len(ind_set)])
    for i, j in not_connected: 
        constraints.append([['x{0}'.format(i), 'x{0}'.format(j)], [1.0] * 2])

    problem.linear_constraints.add(lin_expr = constraints,
                                   senses = [ 'L' ] * (len(independent_sets) + len(not_connected)),
                                   rhs = [1.0] * (len(independent_sets) + len(not_connected)),
                                   names = ['c{0}'.format(x) for x in range(len(independent_sets)+len(not_connected))])
    return problem

#def max_on_inclusion(independent_sets, G):
    #for ind_set in independent_sets:  #в множестве независимых множеств
        #for i in range(len(ind_set)):  #для каждой вершины из независимого множества 
            #cand_in_set(candidates(node, G), ind_set)
            
def get_non_zero_variables(solution): #ненулевые переменные 
    return {i: value for i, value in enumerate(solution) if abs(value) > sys.float_info.epsilon} 

def get_branch_var(non_zero_variables): # вернёт номер вершины для ветвления 
    not_integer_vars_dict = {vertex: value for vertex, value in non_zero_variables.items() 
                                    if abs(value - 1) > sys.float_info.epsilon} 
    if len(not_integer_vars_dict) == 0:
        return None
    return max(not_integer_vars_dict, key=not_integer_vars_dict.get)



def branching(problem):
    
    def branching_constraint(problem, branching_var, rhs):   
        status = problem.linear_constraints.add(lin_expr=[[[branching_var], [1.0]]],
                                                senses=['E'],
                                                rhs=[rhs],
                                                names=['branch_{0}_{1}'.format(branching_var, rhs)]) 
        return problem

    global max_clique
    problem.solve()
    solution = problem.solution.get_values()
    s = sum(solution)
    print(s)
    non_zero = get_non_zero_variables(solution)
    if sum(solution)>len(max_clique):
        var = get_branch_var(non_zero)
        if var is not None:
            print(var)
            branching(branching_constraint(problem, var, 0.0))
            branching(branching_constraint(problem, var, 1.0))
            
        else:
            max_clique = list(non_zero.keys())
            return 0 
    return max_clique
  

def main():
    G = parser('MANN_a9.clq.txt')
    branching(base_model_solution(independent_sets2(G),G.number_of_nodes(),G))
    print("Length of max clique:", len(max_clique))
    print("Nodes of max clique:", max_clique)
    
main()

