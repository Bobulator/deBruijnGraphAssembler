import sys


def write_results(result_file, results):
    if type(results) is list:
        results = "->".join(results)

    with open(result_file, 'w') as f:
        f.write(results)

def read_fasta_file(input_file):
    kmers = []
    with open(input_file, 'r') as f:
        for line in f:
            if line[0] != '>':
                kmers.append(line.strip())

    return kmers


def kmer_counts(reads, k, threshold=0):
    kmer_counts = {}

    for read in reads:
        for kmer in [read[i:i + k] for i in xrange(len(read) - k + 1)]:
            if kmer not in kmer_counts:
                kmer_counts[kmer] = 0

            kmer_counts[kmer] += 1

    results = []
    for kmer, count in kmer_counts.iteritems():
        if count >= threshold:
            for i in xrange(count):
                results.append(kmer)

    return results


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


def eulerian_cycle(edge_dict):
    current_node = edge_dict.keys()[0]
    path = [current_node]

    # Get the initial cycle.
    while True:
        path.append(edge_dict[current_node][0])

        if len(edge_dict[current_node]) == 1:
            del edge_dict[current_node]
        else:
            edge_dict[current_node] = edge_dict[current_node][1:]

        if path[-1] in edge_dict:
            current_node = path[-1]
        else:
            break

    # Continually expand the initial cycle until we're out of edge_dict.
    while len(edge_dict) > 0:
        for i in xrange(len(path)):
            if path[i] in edge_dict:
                current_node = path[i]
                cycle = [current_node]
                while True:
                    cycle.append(edge_dict[current_node][0])

                    if len(edge_dict[current_node]) == 1:
                        del edge_dict[current_node]
                    else:
                        edge_dict[current_node] = edge_dict[current_node][1:]

                    if cycle[-1] in edge_dict:
                        current_node = cycle[-1]
                    else:
                        break

                path = path[:i] + cycle + path[i+1:]
                break
    return path


def eulerian_path(edge_dict):
    # Determine the unbalanced edges.
    out_values = reduce(lambda a, b: a + b, edge_dict.values())
    for node in set(out_values+edge_dict.keys()):
        out_value = out_values.count(node)
        if node in edge_dict:
            in_value = len(edge_dict[node])
        else:
            in_value = 0

        if in_value < out_value:
            unbalanced_from = node
        elif out_value < in_value:
            unbalanced_to = node

    # Add an edge connecting the unbalanced edges.
    if unbalanced_from in edge_dict:
        edge_dict[unbalanced_from].append(unbalanced_to)
    else:
        edge_dict[unbalanced_from] = [unbalanced_to]

    # Get the Eulerian Cycle from the edges, including the unbalanced edge.
    cycle = eulerian_cycle(edge_dict)

    # Find the location of the unbalanced edge in the eulerian cycle.
    divide_point = filter(lambda i: cycle[i:i+2] == [unbalanced_from, unbalanced_to], xrange(len(cycle)-1))[0]

    # Remove the unbalanced edge, and shift appropriately, overlapping the head and tail.
    return cycle[divide_point+1:] + cycle[1:divide_point+1]


def generate_eulerian_path(input_file):

    reads = read_fasta_file(input_file)
    kmers = kmer_counts(reads=reads, k=6)
    graph = build_de_bruijn_graph(kmers)
    #branching_nodes = find_branching_nodes(graph)

    result = eulerian_path(graph)
    return result[0] + ''.join([result[i][-1] for i in xrange(1, len(result))])


if __name__ == "__main__":
    c = generate_eulerian_path(sys.argv[1])
    write_results("tests\\result-output.txt", c)