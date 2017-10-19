# -*- coding: utf-8 -*-
"""
March 2017
@author: Tim Pearce
Estimating network properties using random walks
"""

import networkx as nx
import matplotlib.pyplot as plt
#import pylab
import numpy as np
#from scipy.sparse import csr_matrix
import scipy.sparse as sparse # eigenvalue stuff
import random
import datetime			# for working with dates/times

def fn_transition_m(G, walk_type):
    """ return transtition matrix from a graph as sci sparse matrix
    in preparation for a random walk
    walk_type may be "random", "lazy"
    where lazy, 1/2 chance of staying at same state
    adapted from https://networkx.github.io/documentation/development/_modules/networkx/linalg/
    laplacianmatrix.html#directed_laplacian_matrix
    """
    from scipy.sparse import identity, spdiags
    import scipy as sp
    M = nx.to_scipy_sparse_matrix(G, dtype=float) # get adjacency matrix
    n, m = M.shape
    DI = spdiags(1.0/sp.array(M.sum(axis=1).flat), [0], n, n)
    if walk_type == "random":
        P =  DI * M
    elif walk_type == "lazy":
        I = identity(n)
        P = (I + DI * M) / 2.0
    else:
        raise "walk_type not valid"
    return P


def fn_load_graph(G_type):
    """ load a graph - can be google, gen_large, gen_new """
    if G_type == "google":
        # google graph prep
        G_path = “…\\web-Google.txt"
        G = nx.read_edgelist(G_path,comments='#',create_using=nx.DiGraph(), nodetype=int) # load the graph
        # get largest weakly connected component
        G = max(nx.weakly_connected_component_subgraphs(G),key=len)
        G = G.to_undirected() # convert to undirected
        
    elif G_type == "gen_large":
        # originally used
        #n = 600000 # nodes
        #m = 3 # edges per node
        #G = nx.barabasi_albert_graph(n,m,seed=1) # preferential attachment
        G_path = "random walks\\practical\\keeper results\\graph_01.txt"
        nx.write_edgelist(G, G_path)
        G = nx.read_edgelist(G_path,nodetype=int)
        
    elif G_type == "gen_new":
        n = 50 # nodes
        m = 3 # edges per node
        G = nx.barabasi_albert_graph(n,m,seed=1)
        if draw_graph:
            pos = nx.fruchterman_reingold_layout(G)
            nx.draw(G,pos,with_labels=True)
    
    return G
        
def fn_get_graph_properties():
    """ return actual properties of a graph """
    total_nodes = G.number_of_nodes()
    total_edges = G.number_of_edges()
    total_triangles = int(sum(list(nx.triangles(G).values()))/3)
    return total_nodes, total_edges, total_triangles

def fn_compile_dict():
    """ sometimes a loaded graph will have nodes missing            """
    """ in this case we need to compile a lookup dictionary from    """
    """ the node number to its index                                """
    """ e.g. node 99 might be in index position 87                  """
    """ we need to know this later when working directly with matrices """
    count = 0
    deletions = 0
    list_nodes = G.nodes() # dictionary G.node | index in G.node
    dict_nodes = {} # quicker to compile this stuff first (still takes ages)
    for u in G.nodes():
        count += 1
        if count % 10000 == 0:
            print("\tlookup dictionary creation :: node",count,"of",total_nodes, "\t" ,datetime.datetime.now().strftime('%H:%M:%S'))
            del list_nodes[:9999]
            deletions += 1 # do deletions to increase speed
        dict_nodes[u] = list_nodes.index(u) + deletions*9999 # need manual adjustment for deletions
        
    return dict_nodes


def fn_set_edge_weight():
    """ for weighted random walk we must set edge weights according to the """
    """ property we are estimating """
    if what_est == "edges" or what_est == "cycle_all":
        return # nothing needed -simple random walk
    elif what_est == "nodes":
        for u,v,d in G.edges(data=True):
            d['weight'] = 1/G.degree(u) + 1/G.degree(v) 
    elif what_est == "triangles":   
        adj_m = nx.adjacency_matrix(G)
        count = 0
        for u,v,d in G.edges(data=True):
            count += 1
            if count % 100000 == 0:
                print("\tcalculating edge weights :: edge",count,"of",total_edges, "\t" ,datetime.datetime.now().strftime('%H:%M:%S'))
            u_index = dict_nodes[u]
            v_index = dict_nodes[v]
            u_array = adj_m.getrow(u_index).nonzero()[1] # get lists of adjacent nodes for both nodes
            v_array = adj_m.getrow(v_index).nonzero()[1]
            # find the intersection of the two sets - the size of these is no. triangles
            d['weight'] = 1 + c*len(np.intersect1d(v_array,u_array))

            
def fn_get_trans_m():
    """ get transition matrix and 2 largest eigenvalues """
    """ compile transition matrix into cumulative """
    trans_m = fn_transition_m(G, walk_type) # get transition matrix
    if eig_return == True: # work out eigen value if required
        eig_vals = sparse.linalg.eigs(trans_m, k=2, which="LR",return_eigenvectors=False) 
    else:
        eig_vals = 1 # just dummy to return
    # for each row, scroll through and make cumulative transition matrix
    total_rows = trans_m.shape[1]
    for row in range(0, total_rows): # for each row
        temp_array = trans_m.getrow(row).nonzero()[1] # indices of non zeroes in row
        run_total = 0 # reset
        if row % 100000 == 0: # every x times
            print("\taccumulating trans_m row",row,"of",total_rows, "\t" ,datetime.datetime.now().strftime('%H:%M:%S'))
        for elm in range(0, len(temp_array)): # for each non-zero element (starts from largest)
            temp_val = trans_m[row,temp_array[elm]]
            trans_m[row,temp_array[elm]] = 1 - run_total
            run_total = run_total + temp_val
    
    return trans_m, eig_vals

    
def fn_start_weight():
    """ find the node with the largest starting weight according to """
    """ the property we are searching for """
    if what_est == "edges" or what_est == "nodes" or what_est == "cycle_all":
        # for edges and nodes we choose highest degree
        total_rows = G.number_of_nodes()
        highest_weight = 0
        for i in G.nodes(): # for each row !! google is missing nodes (as trimmed)
            curr_weight = G.degree(i)
            if i % 100000 == 0: # every x returns
                print("\tsearching for largest start weight :: node",i,"of",total_rows, "\t" ,datetime.datetime.now().strftime('%H:%M:%S'))
            if highest_weight < curr_weight:
                start_vertex = i
                start_weight = curr_weight
                highest_weight = curr_weight
        start_vertex_index = G.nodes().index(start_vertex)
        if what_est == "nodes":
            start_weight = 1 # from p229, 3.7
            for val in trans_m.getrow(start_vertex_index).nonzero()[1]:
                start_weight += 1/len(trans_m.getrow(val).nonzero()[1]) 
          
    elif what_est == "triangles":
        total_rows = G.number_of_nodes()
        max_weight = 0
        count = 0
        for u in G.nodes():
            count += 1
            if count % 10000 == 0:
                print("\tsearching for largest start weight :: node",count,"of",total_rows, "\t" ,datetime.datetime.now().strftime('%H:%M:%S'))
            curr_weight = G.degree(u) + 2*nx.triangles(G,u)
            if curr_weight > max_weight:
                start_vertex = u
                start_weight = curr_weight
                max_weight = curr_weight
        start_vertex_index = G.nodes().index(start_vertex)
    
    return start_vertex, start_vertex_index, start_weight

    
def fn_random_walk():
    """ do the damn random walk """
    cycle_fn_running_n = 0
    cycle_fn_running_t = 0
    result_list = []
    no_steps = 0
    return_count = 0
    node_list = G.nodes()
    #start_vertex_index = G.nodes().index(start_vertex)
    curr_vertex = start_vertex_index
    while return_count < return_count_limit:
        no_steps = no_steps + 1
        temp_array = trans_m.getrow(curr_vertex).nonzero()[1] # indices of non zeroes in row
        temp_array = temp_array[::-1] # reverse order
        rnd_no = random.uniform(0, 1)
        for elm in range(0, len(temp_array)): # for each element (starts from smallest)
            if rnd_no <= trans_m[curr_vertex, temp_array[elm]]:
                curr_vertex = temp_array[elm] # select new vertex
                if what_est == "cycle_all":
                    # cycle formula fn----------------
                    curr_degree = len(trans_m.getrow(curr_vertex).nonzero()[1])
                    cycle_fn_running_n += 1/curr_degree
                    curr_no_t = nx.triangles(G,node_list[curr_vertex])
                    cycle_fn_running_t += curr_no_t/(3*curr_degree)
                    # --------------------------------
                if curr_vertex == start_vertex_index: # check if returned
                    return_count = return_count + 1
                    if what_est != "cycle_all":
                        result_list.append([return_count,no_steps])
                    else:
                        result_list.append([return_count,no_steps,cycle_fn_running_n,cycle_fn_running_t])
                    if return_count % 100 == 0: # every x returns
                        print("\trandom walk returns:", return_count,"/", return_count_limit,"\t", datetime.datetime.now().strftime('%H:%M:%S'))
                break
    if what_est == "edges":
        no_edges_est = (no_steps*start_weight)/(2*return_count)
        print("\nactual edges:", total_edges)
        print("no_edges_est:",no_edges_est)
        
    elif what_est == "nodes":
        no_verts_est = (no_steps*start_weight)/(2*return_count) # p229, 3.8
        print("\nactual vertices:", total_nodes)
        print("no_verts_est:",no_verts_est)
        
    elif what_est == "triangles":
        no_triang_est = max((no_steps*start_weight/(6*return_count)) -(G.number_of_edges()/3),0) # p228, 3.6
        print("\nactual triangles:", total_triangles)
        print("no_triang_est:",no_triang_est)
        
    elif what_est == "cycle_all":
        no_verts_est = cycle_fn_running_n * G.degree(start_vertex) / return_count # p227, 3.2
        print("\nactual vertices:", total_nodes)
        print("no_verts_est:",no_verts_est)

        no_trian_est = cycle_fn_running_t * G.degree(start_vertex) / return_count # p227, 3.2
        print("\nactual triangles:", total_triangles)
        print("no_trian_est:", no_trian_est)
            
    return result_list

    
def fn_save_results():
    """ save results to wd - 2 cols for regular, 4 cols for cycle """
    if save_results:
        out_filename = "temp_results" + datetime.datetime.now().strftime('%Y_%m_%d_%H_%M') + "_vandT_est_10k_cycle.txt"
        if what_est == "cycle_all":
            file_output = open(out_filename,"w")
            file_output.write("\n".join(str(row[0])+";"+str(row[1])+";"+str(row[2])+";"+str(row[3]) for row in result_list))
            file_output.close
        else:
            file_output = open(out_filename,"w")
            file_output.write("\n".join(str(row[0])+";"+str(row[1]) for row in result_list))
            file_output.close


def fn_compute_results(result_list):
    """ figure out what the results actually mean """
    #np_results = np.array(result_list)
    
    if what_est != "cycle_all":
        #                      0              1               2              3            4            5
        # results_input = return_count | no_steps_cum | no_steps_diff | estimate_cum | M1_stddev | exp_stddev
        #        6                 7        ---       23             24
        # bootstrap_cum_1 | bootstrap_est_1 --- bootstrap_est_9 | average
        result_list_safe = result_list.copy()
        
        # no_steps_diff - calculate no steps for that return only
        result_list = [x + [0] for x in result_list] # add extra col
        for i in range(0,len(result_list)):
            if i == 0:
                result_list[i][2] = result_list[i][1]
            else:
                result_list[i][2] = result_list[i][1]-result_list[i-1][1]
        
        # estimate_cum - calculate cumulative estimate based on no_steps_cum 
        if what_est == "edges":
            result_list = [x + [(x[1]*start_weight)/(2*x[0])] for x in result_list]
        elif what_est == "nodes":
            result_list = [x + [(x[1]*start_weight)/(2*x[0])] for x in result_list]
        elif what_est == "triangles":
            result_list = [x + [max((x[1]*(start_weight)/(6*x[0]*c)) - (total_edges/(3*c)),0)] for x in result_list]
        
        # M1_stddev - variance using second eigenvalue
        eig_val_2 = np.real(eig_vals[0])
        Zvv = 1/(1-eig_val_2) # p231 3.15
        print("eig_val_2",eig_val_2)
        
        if what_est == "edges":
            stat_var = start_weight/(2*total_edges) # p230 =w(v)/wG
            var_t = (2*Zvv + stat_var - 1) / (stat_var**2) # p231 3.13
            var = var_t * ((start_weight/2)**2)
        elif what_est == "nodes":
            stat_var = start_weight/(2*total_nodes) # p230 =w(v)/wG p229 3.7
            var_t = (2*Zvv + stat_var - 1) / (stat_var**2) # p231 3.13
            var = var_t * ((start_weight/2)**2)
        elif what_est == "triangles":
            stat_var = start_weight/((2*total_edges) + (6*total_triangles)) # p230 =w(v)/wG p228
            var_t = (2*Zvv + stat_var - 1) / (stat_var**2) # p231 3.13
            var = var_t * ((start_weight/(6*c))**2)
            
        std_dev = var**0.5 # stddev = sqrt(var)
        result_list = [x + [(std_dev)/(x[0]**0.5)] for x in result_list] # stddev/sqrt(n)
                                
        # exp_stddev - empirical experimental deviation
        no_steps_diff = [row[2] for row in result_list]
        exp_var_t = np.var(no_steps_diff)
        if what_est == "edges":
            exp_var = exp_var_t * ((start_weight/2)**2)
        elif what_est == "nodes":
            exp_var = exp_var_t * ((start_weight/2)**2)
        elif what_est == "triangles":
            exp_var = exp_var_t * ((start_weight/(6*c))**2)
            
        print("var:",var)
        print("exp_var:",exp_var)
        exp_std_dev = exp_var**0.5 # stddev = sqrt(var)
        result_list = [x + [(exp_std_dev)/(x[0]**0.5)] for x in result_list] # stddev/sqrt(n)
        
        # bootstraps
        random.seed(2)
        for bootstrap_no in range(0,9):
            result_list = [x + [0] for x in result_list] # add extra col
            for i in range(0,len(result_list)): # for each row choose a sample at random, with replacement
                if i == 0:
                    result_list[i][6+(2*bootstrap_no)] = random.sample(no_steps_diff, 1)[0]
                else:
                    result_list[i][6+(2*bootstrap_no)] = random.sample(no_steps_diff, 1)[0]+result_list[i-1][6+(2*bootstrap_no)]
            # estimate_cum - calculate cumulative estimate based on bootstraps_cum 
            if what_est == "edges":
                result_list = [x + [(x[6+(2*bootstrap_no)]*start_weight)/(2*x[0])] for x in result_list]
            elif what_est == "nodes":
                result_list = [x + [(x[6+(2*bootstrap_no)]*start_weight)/(2*x[0])] for x in result_list]
            elif what_est == "triangles":
                result_list = [x + [max((x[6+(2*bootstrap_no)]*(start_weight)/(6*x[0]*c)) - (total_edges/(3*c)),0)] for x in result_list]
        # calc the mean of all
        result_list = [x + [np.mean([ x[3],x[7],x[9],x[11],x[13],x[15],x[17],x[19],x[21],x[23] ])] for x in result_list]

    else: # if cycle method
        #                       0            1              2                      3
        # results_input = return_count | no_steps | cycle_fn_running_n_cum | cycle_fn_running_t_cum
        #      4                5             6                7            8      9
        # running_n_dif | running_t_dif | exp_std_dev_n | exp_std_dev_t | est_n | est_t
        result_list_safe = result_list.copy()
        
        # no_steps_diff - calculate no steps for that return only
        result_list = [x + [0] + [0] for x in result_list] # add 2 extra col
        for i in range(0,len(result_list)):
            if i == 0:
                result_list[i][4] = result_list[i][2]
                result_list[i][5] = result_list[i][3]
            else:
                result_list[i][4] = result_list[i][2]-result_list[i-1][2]
                result_list[i][5] = result_list[i][3]-result_list[i-1][3]
         
        # exp_stddev - empirical experimental deviation
        running_n_dif = [row[4] for row in result_list]
        running_t_dif = [row[5] for row in result_list]
        exp_var_n = np.var(running_n_dif)
        exp_var_t = np.var(running_t_dif)
        exp_var_n = exp_var_n * (start_weight**2)
        exp_var_t = exp_var_t * (start_weight**2)
        print("exp_var_n:",exp_var_n)
        print("exp_var_t:",exp_var_t)
        exp_std_dev_n = exp_var_n**0.5 # stddev = sqrt(var)
        exp_std_dev_t = exp_var_t**0.5 # stddev = sqrt(var)
        result_list = [x + [(exp_std_dev_n)/(x[0]**0.5)] for x in result_list] # stddev/sqrt(n)
        result_list = [x + [(exp_std_dev_t)/(x[0]**0.5)] for x in result_list] # stddev/sqrt(n)
        
        # estimates
        result_list = [x + [(x[2]*start_weight/x[0])] for x in result_list]
        result_list = [x + [(x[3]*start_weight/x[0])] for x in result_list]

    return result_list
    
    
def fn_visualise_results(result_list):
    """ plot the results """
    if what_est != "cycle_all":
        if what_est == "edges":
            actual_property = total_edges
        elif what_est == "nodes":
            actual_property = total_nodes
        elif what_est == "triangles":
            actual_property = total_triangles
            
        no_returns = [row[0] for row in result_list]
        estimate_property = [row[3] for row in result_list]
        m1_var_u = [actual_property+row[4] for row in result_list]
        m1_var_l = [actual_property-row[4] for row in result_list]
        exp_var_u = [actual_property+row[5] for row in result_list]
        exp_var_l = [actual_property-row[5] for row in result_list]
        boot_estimate_1 = [row[7] for row in result_list]
        boot_estimate_2 = [row[9] for row in result_list]
        boot_estimate_3 = [row[11] for row in result_list]
        boot_estimate_4 = [row[13] for row in result_list]
        boot_estimate_5 = [row[15] for row in result_list]
        boot_estimate_6 = [row[17] for row in result_list]
        boot_estimate_7 = [row[19] for row in result_list]
        boot_estimate_8 = [row[21] for row in result_list]
        boot_estimate_9 = [row[23] for row in result_list]
        avg_estimate = [row[24] for row in result_list]
        
        plt.plot(no_returns,estimate_property,color="k",linewidth=1) # original estimate
        plt.plot(no_returns,avg_estimate, color="b",linewidth=2) # averaged estimate
        plt.plot(no_returns,m1_var_u, linestyle="dashdot", color = "r",linewidth=3) # std devs
        plt.plot(no_returns,exp_var_u, linestyle="--", color = "g",linewidth=2)
        plt.plot((min(no_returns), max(no_returns)), (actual_property, actual_property), 'k-',linewidth=2) # actual
        
        # bootstraps
        plt.plot(no_returns,boot_estimate_1,color="k",linestyle="",marker=".", markersize=2,alpha=0.1) 
        plt.plot(no_returns,boot_estimate_2,color="k",linestyle="",marker=".", markersize=2,alpha=0.1)
        plt.plot(no_returns,boot_estimate_3,color="k",linestyle="",marker=".", markersize=2,alpha=0.1)
        plt.plot(no_returns,boot_estimate_4,color="k",linestyle="",marker=".", markersize=2,alpha=0.1)
        plt.plot(no_returns,boot_estimate_5,color="k",linestyle="",marker=".", markersize=2,alpha=0.1)
        plt.plot(no_returns,boot_estimate_6,color="k",linestyle="",marker=".", markersize=2,alpha=0.1)
        plt.plot(no_returns,boot_estimate_7,color="k",linestyle="",marker=".", markersize=2,alpha=0.1)
        plt.plot(no_returns,boot_estimate_8,color="k",linestyle="",marker=".", markersize=2,alpha=0.1)
        plt.plot(no_returns,boot_estimate_9,color="k",linestyle="",marker=".", markersize=2,alpha=0.1)    
        plt.plot(no_returns,m1_var_l, linestyle="dashdot", color = "r",linewidth=3)
        plt.plot(no_returns,avg_estimate, color="b",linewidth=2) # averaged estimate # repeat
        plt.plot(no_returns,exp_var_u, linestyle="--", color = "g",linewidth=2) # repeat
        plt.plot(no_returns,exp_var_l, linestyle="--", color = "g",linewidth=2)
        plt.ylim([min(avg_estimate),max(avg_estimate)])
        plt.xlim([min(no_returns),max(no_returns)])
        plt.xscale('log')
        #plt.legend(['orig est','avg est', 'std dev : eigen', 'std dev : exp', 'act. value'], bbox_to_anchor=(1.1, 1.1))
        fig = plt.gcf()
        plt.show()
        fig.set_size_inches(6, 4)
        fig.savefig('google_rw_.png', dpi=300)
        
    else: # if is cycle_all
        # first plot nodes
        actual_property = total_nodes
        
        no_returns = [row[0] for row in result_list]
        estimate_property = [row[8] for row in result_list]
        exp_var_u = [actual_property+row[6] for row in result_list]
        exp_var_l = [actual_property-row[6] for row in result_list]
        
        plt.plot(no_returns,estimate_property,color="k",linewidth=1) # original estimate
        #plt.plot(no_returns,avg_estimate, color="b",linewidth=2) # averaged estimate
        plt.plot(no_returns,exp_var_u, linestyle="--", color = "g",linewidth=2)
        plt.plot((min(no_returns), max(no_returns)), (actual_property, actual_property), 'k-',linewidth=2) # actual
        plt.plot(no_returns,exp_var_l, linestyle="--", color = "g",linewidth=2)

        plt.ylim([min(estimate_property),max(estimate_property)])
        plt.xlim([min(no_returns),max(no_returns)])
        plt.xscale('log')
        fig = plt.gcf()
        plt.show()
        fig.set_size_inches(6, 4)
        fig.savefig('cycle_nodes.png', dpi=300)
        
        # now plot triangles
        actual_property = total_triangles
        
        estimate_property = [row[9] for row in result_list]
        exp_var_u = [actual_property+row[7] for row in result_list]
        exp_var_l = [actual_property-row[7] for row in result_list]
        
        plt.plot(no_returns,estimate_property,color="k",linewidth=1) # original estimate
        #plt.plot(no_returns,avg_estimate, color="b",linewidth=2) # averaged estimate
        plt.plot(no_returns,exp_var_u, linestyle="--", color = "g",linewidth=2)
        plt.plot((min(no_returns), max(no_returns)), (actual_property, actual_property), 'k-',linewidth=2) # actual
        plt.plot(no_returns,exp_var_l, linestyle="--", color = "g",linewidth=2)

        plt.ylim([min(estimate_property),max(estimate_property)])
        plt.xlim([min(no_returns),max(no_returns)])
        plt.xscale('log')
        fig = plt.gcf()
        plt.show()
        fig.set_size_inches(6, 4)
        fig.savefig('cycle_triang.png', dpi=300)
    
    
# =============================================================================
# inputs

G_type = "gen_new" # "google" "gen_large" "gen_new"
what_est = "cycle_all" # "nodes" "edges" "triangles" "cycle_all"
draw_graph = False # True False
eig_return = True # True False
walk_type = "random" # "random" "lazy"
return_count_limit = 500
save_results = False # True False
c = 1 # for triangles

# =============================================================================
# program

print("\n-- 00 started: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
print("-- 01 loading graph: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
G = fn_load_graph(G_type) # load graph
print("-- 02 calculating properties: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
total_nodes, total_edges, total_triangles = fn_get_graph_properties() # properties
print("-- 03 creating lookup dict: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
dict_nodes = fn_compile_dict() # lookup dict from node number to index number
print("-- 04 setting edge weights: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
fn_set_edge_weight() # set edge weights
print("-- 05 getting transition matrix and eigenvalues: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
trans_m, eig_vals = fn_get_trans_m() # get transition matrix and eigenvalues
print("-- 06 finding start vertex: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
start_vertex, start_vertex_index, start_weight = fn_start_weight()
print("-- 07 doing random walk: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
result_list = fn_random_walk()
print("-- 08 saving results: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
fn_save_results()
print("-- 09 computing results: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
result_list = fn_compute_results(result_list)
print("-- 10 visualising results: ",datetime.datetime.now().strftime('%H:%M:%S'),"--")
fn_visualise_results(result_list)
