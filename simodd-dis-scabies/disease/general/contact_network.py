import numpy as np
import matplotlib as plt

try:
    import networkx as nx
except ImportError:
    pass
from population.utils import sample


class ContactNetwork(object):
    """
    A community contact network for the current population
    
    ind_ages: a map of (ID, age) pairs
    """

    def __init__(self):
        self.G = nx.Graph()

    def create_network(self, pop, cmatrix, rng):
        """
        Populate contact network nodes and edges.

        So far as possible, edges are created such that:
    
        (a) the degree of node x equals the expected number
            of contacts given x's age; and
    
        (b) the propensity for node x and y to share an edge
            is relative to the level of contact between x and
            y, given their respective ages.
        """
        print("creating contact network...")

        # create a disconnected graph to hold the new network (each node = an individual in the population)
        self.G.add_nodes_from(list(pop.I.keys()))

        # create a dictionary mapping (lower bound of) age group to a list of individuals in that age group
        I_by_age = {}
        for cur_min_age, cur_max_age in zip(cmatrix.age_classes, cmatrix.age_classes[1:] + [102]):
            I_by_age[cur_min_age] = pop.individuals_by_age(cur_min_age, cur_max_age - 1)

        # create a dictionary mapping node ID to the number of edges remaining to be allocated for that node 
        remaining_edges = dict(
            [(x[0], int(cmatrix.comm_a_levels[cmatrix.age_map[x[1].age]]))
             for x in pop.I.items()])

        # print remaining edges structure (for validation)
        #for i in remaining_edges:
        #    print i, pop.I[i].age, remaining_edges[i]

        # get a matrix of age-specific contact rates rescaled such that each row sums to 1
        n_matrix = cmatrix.get_normalised_matrix()

        # create a shuffled list of node IDs to iterate over
        node_ids = list(pop.I.keys())
        rng.shuffle(node_ids)

        # for each source node (current indivdiual)
        for source in node_ids:
            #print "allocating edges for node", cur_ind, "--", pop.I[cur_ind].age, "years"

            # skip iteration if cur_ind already has all edges allocated
            if remaining_edges[source] <= 0:
                continue

            # allocate requisite number of remaining edges
            # keep track of how many attempts we have made to allocate this edge;
            # if maximum is exceeded 
            remaining_edge_attempts = 100
            while remaining_edges[source] > 0:
                remaining_edge_attempts -= 1
                if remaining_edge_attempts < 0:
                    break
                #print "remaining edges:", remaining_edges[cur_ind]

                # get target age_class index
                sample_empty = True
                fail = False
                remaining_attempts = 100
                while sample_empty:
                    remaining_attempts -= 1
                    # sample age of target node from matrix of normalised age x age mixing rates
                    target_age_class = cmatrix.age_classes[sample(n_matrix[pop.I[source].age], rng)]
                    # check that there is at least one potential target node in the resulting sample
                    sample_empty = (len(I_by_age[target_age_class]) == 0)
                    if remaining_attempts < 0:
                        fail = True
                        break
                if fail:
                    continue

                #print "target age class:", target_age_class

                # randomly choose an individual from target sample (who is not self and not already a neighbour)
                target = source
                #print pop.I[cur_ind].age
                remaining_attempts = 100
                while target == source or target in nx.all_neighbors(self.G, source):
                    remaining_attempts -= 1
                    # check that sample contains more than just the source node!
                    if [pop.I[source]] == I_by_age[target_age_class] or remaining_attempts < 0:
                        target = -1
                        break
                    target = rng.sample(I_by_age[target_age_class], 1)[0].ID
                #print cur_ind, "considering target:", target, "(age_class:", \
                #    target_age_class, ") : age :", pop.I[target].age, \
                #    ": remaining edges :", remaining_edges[target]
                #print "target's current edges", self.G.edges(target)

                # if successfully identified a target
                if target > 0:
                    # create the new edge
                    self.G.add_edge(source, target)
                    #print "added edge: (", cur_ind, ",", target, ")"

                    # update remaining edges of source and target
                    remaining_edges[source] -= 1
                    remaining_edges[target] -= 1

                    #print "remaining edges for target", target, ":", remaining_edges[target]
                    # if target has all their edges filled, remove from available nodes
                    if remaining_edges[target] <= 0:
                        #print "target has no remaining edges: index", cmatrix.age_ids_map[pop.I[target].age]
                        I_by_age[cmatrix.age_ids_map[pop.I[target].age]].remove(pop.I[target])

                        # remove cur_ind from available nodes

                        #print cur_ind, "finished, removing from I_by_age, index:", \
                        #    cmatrix.age_ids_map[pop.I[cur_ind].age]
            I_by_age[cmatrix.age_ids_map[pop.I[source].age]].remove(pop.I[source])

            #print "TOTAL: after index", node_ids.index(cur_ind), "remaining edges:", sum(remaining_edges.values())
            #the last few edges are typically very difficult to pair;
            #get quicker results (and acceptable networks) by ignoring them...
            if sum(remaining_edges.values()) <= 10:
                break

        print("contact network created!")

    def get_contacts(self, pop, ind):
        """
        Returns contacts of specified individual.
        """
        return [pop.I[x] for x in self.G.neighbors(ind.ID)]

    def write_age_matrix(self, pop):
        """
        Dump the age matrix resulting from the contact network to a file.
        """
        age_matrix = np.zeros((101, 101))

        for cur_edge in self.G.edges_iter():
            age_1 = pop.I[cur_edge[0]].age
            age_2 = pop.I[cur_edge[1]].age
            age_matrix[age_1][age_2] += 1
            age_matrix[age_2][age_1] += 1

        age_matrix /= 2.0

        np.savetxt('age_matrix.csv', age_matrix, delimiter=',')

    def dump_stats(self):
        """
        Dump various statistics on the contact network structure.
        """
        print("Is connected?", nx.is_connected(self.G))
        print("Number connected components:", nx.number_connected_components(self.G))
        print("Size of connected components:", [len(x) for x in nx.connected_components(self.G)])

        degrees = list(self.G.degree().values())
        degree_dist = np.histogram(degrees, bins=25, range=(0, 25), normed=False)
        print(degree_dist[0])
        print(degree_dist[1])

        pos = nx.spring_layout(self.G)
        nx.draw(self.G, pos)
        plt.savefig('cnetwork.pdf')
        nx.draw_graphviz(self.G)
        nx.write_dot(self.G, 'cnetwork.dot')
