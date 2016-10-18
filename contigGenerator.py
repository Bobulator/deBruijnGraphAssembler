import sys


def read_fasta_file(input_file):
    return ""


def get_kmers(sequence, k):
    return [sequence[i:k] for i in xrange(len(sequence) - k + 1)]


def build_de_bruijn_graph(kmers):
    graph = {}

    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]

        if prefix not in graph:
            graph[prefix] = []

        graph[prefix].append(suffix)

    return graph


def find_branching_nodes(graph):
    edge_counts = {}

    for node, edges in graph.iteritems():

        if node not in edge_counts:
            edge_counts[node] = [0, 0]

        edge_counts[node][1] += len(edges)

        for edge in edges:
            if edge not in edge_counts:
                edge_counts[edge] = [0, 0]
            edge_counts[edge][0] += 1

    return [node for node, edges in edge_counts.iteritems() if not (edges[0] == 1 and edges[1] == 1)]


def generate_contigs(input_file, k):
    sequence = read_fasta_file(input_file)
    kmers = get_kmers(sequence, k)
    graph = build_de_bruijn_graph(kmers)
    branching_nodes = find_branching_nodes(graph)

    contigs = []
    for node in branching_nodes:
        if node not in graph:
            continue

        for edge in graph[node]:
            contig = list(node)

            current_node = edge
            while current_node not in branching_nodes:
                contig.append(current_node[-1])
                current_node = graph[current_node][0]

            contig.append(current_node[-1])
            contigs.append(''.join(contig))

    contigs.sort()
    return contigs


if __name__ == "__main__":
    generate_contigs(sys.argv[1], int(sys.argv[2]))