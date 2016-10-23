import unittest
from contigGenerator import *


class ContigGeneratorTest(unittest.TestCase):

    def test_build_weighted_de_bruijn_graph(self):

        # Should successfully build de bruijn graph with edge weights
        kmers = ['ACGT', 'CGTA', 'GTAC', 'TACG', 'ACGT', 'ACGA']
        weighted_edge_graph = build_weighted_de_bruijn_graph(kmers)

        # Verify 4 nodes: ACG, CGT, GTA, TAC
        self.assertEqual(4, len(weighted_edge_graph))

        # Verify exit edges
        self.assertEqual(2, len(weighted_edge_graph['ACG']))
        self.assertEqual(1, len(weighted_edge_graph['CGT']))
        self.assertEqual(1, len(weighted_edge_graph['GTA']))
        self.assertEqual(1, len(weighted_edge_graph['TAC']))

        # Verify exit edge weights
        self.assertEqual(2, weighted_edge_graph['ACG']['CGT'])
        self.assertEqual(1, weighted_edge_graph['ACG']['CGA'])
        self.assertEqual(1, weighted_edge_graph['CGT']['GTA'])
        self.assertEqual(1, weighted_edge_graph['GTA']['TAC'])
        self.assertEqual(1, weighted_edge_graph['TAC']['ACG'])

    def test_filter_weighted_de_bruijn_graph(self):

        # TEST 1

        # Should not filter when edge_filter is 0
        edge_filter = 0
        weighted_edge_graph = {'A': {'A': 1, 'C': 2, 'G': 3}, 'C': {'A': 1, 'C': 2, 'G': 3, 'T':5}}
        filtered_graph = filter_weighted_de_bruijn_graph(edge_filter, weighted_edge_graph)

        # Should correctly convert to non-weighted graph
        self.assertEqual(2, len(filtered_graph))
        self.assertEqual(3, len(filtered_graph['A']))
        self.assertEqual(4, len(filtered_graph['C']))

        filtered_graph['A'].sort()
        filtered_graph['C'].sort()
        self.assertEqual(['A', 'C', 'G'], filtered_graph['A'])
        self.assertEqual(['A', 'C', 'G', 'T'], filtered_graph['C'])

        # TEST 2

        # Should filter all edges with 1
        edge_filter = 1
        weighted_edge_graph = {'A': {'A': 1, 'C': 2, 'G': 3}, 'C': {'A': 1, 'C': 2, 'G': 3, 'T': 5}}
        filtered_graph = filter_weighted_de_bruijn_graph(edge_filter, weighted_edge_graph)

        self.assertEqual(2, len(filtered_graph))
        self.assertEqual(2, len(filtered_graph['A']))
        self.assertEqual(3, len(filtered_graph['C']))

        filtered_graph['A'].sort()
        filtered_graph['C'].sort()
        self.assertEqual(['C', 'G'], filtered_graph['A'])
        self.assertEqual(['C', 'G', 'T'], filtered_graph['C'])

        # TEST 3

        # Should filter all edges <= 2
        edge_filter = 2
        weighted_edge_graph = {'A': {'A': 1, 'C': 2, 'G': 3}, 'C': {'A': 1, 'C': 2, 'G': 3, 'T': 5}}
        filtered_graph = filter_weighted_de_bruijn_graph(edge_filter, weighted_edge_graph)

        self.assertEqual(2, len(filtered_graph))
        self.assertEqual(1, len(filtered_graph['A']))
        self.assertEqual(2, len(filtered_graph['C']))

        filtered_graph['A'].sort()
        filtered_graph['C'].sort()
        self.assertEqual(['G'], filtered_graph['A'])
        self.assertEqual(['G', 'T'], filtered_graph['C'])

        # TEST 4

        # Should filter all edges <= 3
        edge_filter = 3
        weighted_edge_graph = {'A': {'A': 1, 'C': 2, 'G': 3}, 'C': {'A': 1, 'C': 2, 'G': 3, 'T': 5}}
        filtered_graph = filter_weighted_de_bruijn_graph(edge_filter, weighted_edge_graph)

        self.assertEqual(2, len(filtered_graph))
        self.assertEqual(0, len(filtered_graph['A']))
        self.assertEqual(1, len(filtered_graph['C']))

        filtered_graph['A'].sort()
        filtered_graph['C'].sort()
        self.assertEqual([], filtered_graph['A'])
        self.assertEqual(['T'], filtered_graph['C'])

        # TEST 5

        # Should filter all edges <= 4
        edge_filter = 4
        weighted_edge_graph = {'A': {'A': 1, 'C': 2, 'G': 3}, 'C': {'A': 1, 'C': 2, 'G': 3, 'T': 5}}
        filtered_graph = filter_weighted_de_bruijn_graph(edge_filter, weighted_edge_graph)

        self.assertEqual(2, len(filtered_graph))
        self.assertEqual(0, len(filtered_graph['A']))
        self.assertEqual(1, len(filtered_graph['C']))

        filtered_graph['A'].sort()
        filtered_graph['C'].sort()
        self.assertEqual([], filtered_graph['A'])
        self.assertEqual(['T'], filtered_graph['C'])

        # TEST 6

        # Should filter all edges <= 5
        edge_filter = 5
        weighted_edge_graph = {'A': {'A': 1, 'C': 2, 'G': 3}, 'C': {'A': 1, 'C': 2, 'G': 3, 'T': 5}}
        filtered_graph = filter_weighted_de_bruijn_graph(edge_filter, weighted_edge_graph)

        self.assertEqual(2, len(filtered_graph))
        self.assertEqual(0, len(filtered_graph['A']))
        self.assertEqual(0, len(filtered_graph['C']))

        filtered_graph['A'].sort()
        filtered_graph['C'].sort()
        self.assertEqual([], filtered_graph['A'])
        self.assertEqual([], filtered_graph['C'])


if __name__ == "__main__":
    unittest.main()