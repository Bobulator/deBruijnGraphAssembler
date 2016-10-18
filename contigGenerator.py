import sys

def write_results(result_file, results):
    if type(results) is list:
        results = "\n".join(results)

    with open(result_file, 'w') as f:
        f.write(results)

def read_fasta_file(input_file):
    kmers = []
    with open(input_file, 'r') as f:
        for line in f:
            if line[0] != '>':
                kmers.append(line.strip())

    return kmers


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

    for node, edges in graph.items():

        if node not in edge_counts:
            edge_counts[node] = [0, 0]

        edge_counts[node][1] += len(edges)

        for edge in edges:
            if edge not in edge_counts:
                edge_counts[edge] = [0, 0]
            edge_counts[edge][0] += 1

    return [node for node, edges in edge_counts.items() if not (edges[0] == 1 and edges[1] == 1)]


def generate_contigs(input_file):
    kmers = read_fasta_file(input_file)
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
    c = generate_contigs(sys.argv[1])
    write_results("tests\\result-output.txt", c)