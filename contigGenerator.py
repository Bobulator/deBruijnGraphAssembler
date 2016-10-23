import numpy
import sys
import os
import errno


def convert_number_to_string(num):
    s = ""
    if num < 10:
        s += "0"
    s += str(num)
    return s


def calculate_median_contig_size(contigs):
    contig_lengths = []

    for contig in contigs:
        contig_lengths.append(len(contig))

    return numpy.median(numpy.array(contig_lengths))


def calculateN50(contigs):
        total_len_of_contigs = 0
        for contig in contigs:
            total_len_of_contigs += len(contig)

        halfway = total_len_of_contigs / 2

        for i in range(0, len(contigs)):
            contig = contigs[i]
            if halfway - len(contig) <= 0:
                return len(contig)

            halfway -= len(contig)

        return -1


def write_results(result_file, longest_contig_length, median_contig_size, n50, results):
    if type(results) is list:
        results = "\n".join(results)

    if not os.path.exists(os.path.dirname(result_file)):
        try:
            os.makedirs(os.path.dirname(result_file))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    with open(result_file, 'w') as f:
        f.write("***** STATISTICS *****\n")
        f.write("\tLongest contig length: " + str(longest_contig_length) + "\n")
        f.write("\tMedian contig size: " + str(median_contig_size) + "\n")
        f.write("\tN50: " + str(n50) + "\n")
        f.write("\n\n\n")
        f.write("***** CONTIGS *****\n")
        f.write(results)

        
def read_fasta_file(input_file):
    kmers = []
    with open(input_file, 'r') as f:
        for line in f:
            if line[0] != '>':
                kmers.append(line.strip())

    return kmers


def kmer_counts(reads, k, threshold=0):
    kmer_occurence_per_read = {}

    for index, read in enumerate(reads):
        for kmer in [read[i:i + k] for i in xrange(len(read) - k + 1)]:
            if kmer not in kmer_occurence_per_read:
                kmer_occurence_per_read[kmer] = {}

            if index not in kmer_occurence_per_read[kmer]:
                kmer_occurence_per_read[kmer][index] = 0

            kmer_occurence_per_read[kmer][index] += 1

    results = []
    for kmer, occurence_per_read in kmer_occurence_per_read.iteritems():
        if len(occurence_per_read) >= threshold:
            for _, count in kmer_occurence_per_read[kmer].iteritems():
                results.append(kmer)

    return results


def build_weighted_de_bruijn_graph(kmers):
    graph = {}

    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]

        if prefix not in graph:
            graph[prefix] = {}

        if suffix not in graph[prefix]:
            graph[prefix][suffix] = 0

        graph[prefix][suffix] += 1

    return graph


def filter_weighted_de_bruijn_graph(filter, graph):
    for node, weighted_edges in graph.items():
        for edge, weight in weighted_edges.items():
            if weight <= filter:
                weighted_edges.pop(edge)

        graph[node] = list(weighted_edges.keys())

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

    return [node for node, degrees in edge_counts.items() if not (degrees[0] == 1 and degrees[1] == 1)]


def generate_contigs(input_file, k, coverage_filter, weight_filter):
    reads = read_fasta_file(input_file)
    kmers = kmer_counts(reads, k, coverage_filter)
    print kmers
    graph = build_weighted_de_bruijn_graph(kmers)
    print graph
    graph = filter_weighted_de_bruijn_graph(weight_filter, graph)
    print graph
    branching_nodes = find_branching_nodes(graph)

    contigs = []
    for node in branching_nodes:
        if node not in graph:
            continue

        for exit_edge in graph[node]:
            contig = list(node)

            current_node = exit_edge
            while current_node not in branching_nodes:
                contig.append(current_node[-1])
                current_node = graph[current_node][0]

            contig.append(current_node[-1])
            contigs.append(''.join(contig))

    contigs.sort()
    return list(set(contigs))


def generate_output(contigs, k, coverage_filter, weighted_edge_filter):
    # Sort contigs first
    contigs.sort(key=len, reverse=True)

    if len(contigs) == 0:
        raise Exception("No contigs were found")

    # Calculate length of longest contig
    longest_contig_length = len(contigs[0])

    # Calculate median contig size
    median_contig_size = calculate_median_contig_size(contigs)

    # Calculate N50
    n50 = calculateN50(contigs)

    output_file = "output\\" + filename_only + "-" + convert_number_to_string(
        kmer_len) + "-" + convert_number_to_string(coverage_filter) + "-" + convert_number_to_string(
        weighted_edge_filter) + ".txt"
    write_results(output_file, longest_contig_length, median_contig_size, n50, contigs)

'''
def assembleContigs(contigs, k):

    for kmer_size in xrange(k, 0, -1):

        len_changed = True
        while len_changed:
            temp = len(contigs)

            for i in xrange(len(contigs)):
                if not contigs[i]:
                    continue

                for j in [x for x in xrange(len(contigs)) if x != i]:
                    if contigs[i][-kmer_size + 1:] == contigs[j][:kmer_size - 1]:
                        contigs[i] += contigs[j][kmer_size - 1:]
                        contigs[j] = []

            contigs = [x for x in contigs if x]
            print len(contigs)
            if temp == len(contigs):
                len_changed = False

    return [''.join(x) for x in contigs]
'''


if __name__ == "__main__":
    filename = sys.argv[1]
    filename_only = filename.split(os.sep)[-1]
    kmer_len = int(sys.argv[2])
    coverage_filter = int(sys.argv[3])
    weighted_edge_filter = int(sys.argv[4])

    # Calculate the contigs
    contigs = generate_contigs(filename, kmer_len, coverage_filter, weighted_edge_filter)
    generate_output(contigs, kmer_len, coverage_filter, weighted_edge_filter)

