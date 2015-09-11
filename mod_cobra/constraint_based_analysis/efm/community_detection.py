from igraph import Graph
from louvain import find_partition

__author__ = 'anna'

g = Graph()
g.add_vertices([0, 1, 2, 3, 4, 5])
for (i, j) in [(0, 1), (0, 2), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5)]:
    g.add_edge(i, j)
    g.add_edge(j, i)
partition = find_partition(graph=g, method='Significance')
print partition