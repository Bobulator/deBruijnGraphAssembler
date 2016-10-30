import sys
import os
import errno


def convert_number_to_string(num):
    s = ""
    if num < 10:
        s += "0"
    s += str(num)
    return s


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


def write_results(results_file_name, contigs_file_name, longest_contig_length, num_contigs, n50, results):

    if not os.path.exists(os.path.dirname(contigs_file_name)):
        try:
            os.makedirs(os.path.dirname(contigs_file_name))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    # This appends a line to the results file, which
    #  is the same as a row in a spreadsheet
    # The file columns are as follows:
    # run name,read length,kmer length,coverage filter, weighted edge filter,longest contig length,number of contigs,n50
    #  TODO: Use the actual read length
    with open(results_file_name, 'a') as f:
        run_name = contigs_file_name.replace('contigs-', '')
        run_name = run_name.replace('.txt', '')
        f.write(run_name)
        f.write(',')
        f.write('read_length')
        f.write(',')
        f.write(str(kmer_len))
        f.write(',')
        f.write(str(coverage_filter))
        f.write(',')
        f.write(str(weighted_edge_filter))
        f.write(',')
        f.write(str(longest_contig_length))
        f.write(',')
        f.write(str(num_contigs))
        f.write(',')
        f.write(str(n50))
        f.write('\n')
        f.close()

    # Write contigs to file
    with open(contigs_file_name, 'w') as f:
        if type(results) is list:
            results = "\n".join(results)

        f.write(results)
        f.close()

        
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
                for c in xrange(count):
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
    graph = build_weighted_de_bruijn_graph(kmers)
    graph = filter_weighted_de_bruijn_graph(weight_filter, graph)
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


def generate_output(filename_only, contigs, coverage_filter, weighted_edge_filter):
    # Sort contigs first
    contigs.sort(key=len, reverse=True)

    if len(contigs) == 0:
        raise Exception("No contigs were found")

    # Calculate length of longest contig
    longest_contig_length = len(contigs[0])

    # Calculate number of contigs
    num_contigs = len(contigs)

    # Calculate N50
    n50 = calculateN50(contigs)


    # This file contains all of the accumulated results
    # for all of the runs of this gene
    results_file_name = "output\\" + filename_only + "-results.csv"
    # This file contains all of the contigs created by this run
    contigs_output_file_name = "output\\" + filename_only + "-contigs-" + convert_number_to_string(
        kmer_len) + "-" + convert_number_to_string(coverage_filter) + "-" + convert_number_to_string(
        weighted_edge_filter) + ".txt"
    write_results(results_file_name, contigs_output_file_name, longest_contig_length, num_contigs, n50, contigs)

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
    # filename_only = filename.split('/')[-1]
    filename_only = filename.split(os.sep)[-1]
    kmer_len = int(sys.argv[2])
    coverage_filter = int(sys.argv[3])
    weighted_edge_filter = int(sys.argv[4])

    # Calculate the contigs
    contigs = generate_contigs(filename, kmer_len, coverage_filter, weighted_edge_filter)
    generate_output(filename_only, contigs, coverage_filter, weighted_edge_filter)

