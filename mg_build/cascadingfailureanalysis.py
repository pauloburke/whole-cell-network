import numpy as np

def cascade_failure(g,node,target_reaction='',verbose=False):
    if type(node) != list:
        node = [node]
    for n in node:
        if n not in g:
            print("WARNING: node not found in network:",n)
            return False,False,False
    vertices_counter = []
    cascade_counter = 0
    reactions = []
    deleted = {n:False for n in g.nodes()}
    deleted[''] = False
    deleted_edges = {n:0 for n in g.nodes()}
    for n in node:
        deleted[n] = True
        if verbose:
            print("Deleted node:",n,"type:",g.nodes[n]['type'])
        if g.nodes[n]['type'] == 'm':
            reacts = [key for key in g.successors(n)]
            for r in reacts:
                reactions.append(r)
        else:
            reactions.append(n)
    while(len(reactions)):
        reactions = list(set(reactions))
        if verbose:
            print("cascade:",cascade_counter,"\tReactions:",len(reactions))        
        molecules = []
        for r in reactions:
            products = [key for key in g.successors(r)]
            for m in products:
                if not deleted[m]:
                    if 'biosynthesisreaction' in g.nodes[m].keys():
                        if r in g.nodes[m]['biosynthesisreaction']:
                            molecules.append(m)
                            continue
                    deleted_edges[m] += 1
                    if (g.nodes[m]['indegree']-deleted_edges[m]) == 0:
                        molecules.append(m)
            deleted[r] = True
        molecules = list(set(molecules))
        reactions = []
        for m in molecules:
            mol_reactions = [key for key in g.successors(m)]
            for r in mol_reactions:
                if not deleted[r]:
                    reactions.append(r)
            deleted[m] = True
        cascade_counter += 1
    final_g = g.subgraph([key for key,val in deleted.items() if not val])
    return len(g)-len(final_g),deleted[target_reaction],final_g

def randomize_bipartite_network(g,p,mol_nodes,react_nodes):    
    to_remove = [e for e in g.edges() if np.random.rand()<p and (e[0] in react_nodes or e[1] in react_nodes)]
    g.remove_edges_from(to_remove)
    for e in to_remove:
        source = 0
        target = 0
        if e[0] in mol_nodes:
            source = mol_nodes[np.random.randint(len(mol_nodes))]
            target = react_nodes[np.random.randint(len(react_nodes))]
        else:
            target = mol_nodes[np.random.randint(len(mol_nodes))]
            source = react_nodes[np.random.randint(len(react_nodes))]
        g.add_edge(source,target)
