# The program is given a fasta file with DNA sequences, a TRANSFAC matrix table, an 'A' prior and a score
# threshold for the motif. The motif in the TRANSFAC table can be of varying length, as can the DNA sequences.
# The output should be the sites that score above the provided threshold (where on which sequence).
# Site scores are calculated as the sum of PSSM weights selected at each position.

import sys
import math
# get TRANSFAC matrix
# define a function that open the transfac file and extract the matrix from the file
def transfac_matrix(filename):
    # catch file opening error
    try:
        infile = open(filename, 'r')
    except IOError as error:
        print('Can’t open file, reason:', str(error))
        sys.exit(1)
    # define a flag to indicate the position of transfac matrix
    flag_matrix = False
    # initialise the transfac_matrix as an empty list
    transfac_matrix = []
    for line in infile:
        while flag_matrix:
            # end the loop when meet XX
            if line.startswith('XX'):
                flag_matrix = False
                break
            else:
                # extract the positional base count in the matrix
                transfac_line = line.split()[1:len(transfac_matrix[0])+1]
                # convert all the str into int type
                transfac_line = list(map(int, transfac_line))
                transfac_matrix.append(transfac_line)
                break
        if line.startswith('PO'):
            transfac_matrix.append(line.split()[1:])
            flag_matrix = True
    infile.close()
    return transfac_matrix


# define a function that convert from a positional base count matrix to a log-likelihood matrix
def ll_matrix(transfac_matrix, p_a):
    # calculate the priors using given the prior for 'A'
    p_t = p_a
    p_c = p_g = (1-p_a*2)/2
    # set a dictionary to save the priors
    prior_dic = {'A': p_a, 'T': p_t, 'C': p_c, 'G': p_g}
    # introduce pseudo-counts to the matrix
    for i in range(1, len(transfac_matrix)):
        for j in range(len(transfac_matrix[0])):
            transfac_matrix[i][j] += prior_dic[transfac_matrix[0][j]]
    # initialise the log-likelihood matrix
    ll_matrix = [['A', 'C', 'G', 'T']]
    # convert the positional base count matrix into probability matrix and then into log-likelihood matrix
    for i in range(1, len(transfac_matrix)):
        ll_matrix.append([])
        for j in range(len(transfac_matrix[0])):
            prob = transfac_matrix[i][j] / sum(transfac_matrix[i])
            ll_prob = format(math.log2(prob / prior_dic[ll_matrix[0][j]]), '.4f')
            ll_matrix[i].append(ll_prob)
    print(ll_matrix)
    return ll_matrix


# define a function that compute the score of the given sequence using the log-likelihood matrix
def score_compute(sequence, ll_matrix):
    score_sum = 0
    for i in range(len(sequence)):
        score = ll_matrix[i+1][(ll_matrix[0].index(sequence[i]))]
        score = float(score)
        score_sum += score
    return score_sum


# use the log-likelihood matrix on a sequence
if len(sys.argv) <= 4:
    print('Usage:', sys.argv[0], '<fasta filename>', '<TRANSFAC filename>', '<A prior>', '<score threshold>')
    sys.exit(2)
else:
    ll_matrix = ll_matrix(transfac_matrix(sys.argv[2]), float(sys.argv[3]))
    thold = float(sys.argv[4])
    # catch file opening error
    try:
        infile = open(sys.argv[1], 'r')
    except IOError as error:
        print('Can’t open file, reason:', str(error))
        sys.exit(1)
    # save the accession number into acc_num_list
    # save the sequence into seq_list
    motif_len = len(ll_matrix) - 1
    acc_num_list = []
    seq_list = []
    seq_idx = -1
    for line in infile:
        if line.startswith('>'):
            acc_num_list.append(line.split()[0])
            seq_idx = seq_idx + 1
            seq_list.append("")
        else:
            seq_list[seq_idx] = seq_list[seq_idx] + line.strip()
    # iterate all the sequences and calculate the score
    # output the sites that scores above the given threshold
    for idx in range(len(acc_num_list)):
        motif_pos = 0
        while motif_pos <= len(seq_list[idx]) - motif_len:
            score_sum = score_compute(seq_list[idx][motif_pos:motif_pos + motif_len], ll_matrix)
            if score_sum >= thold:
                print(acc_num_list[idx], ':', motif_pos)
            motif_pos += 1
    infile.close()
