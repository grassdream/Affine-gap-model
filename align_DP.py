# coding=gbk
# from io_functions import *
# import os
import sys
import time
import math
import numpy as np


# function to write output.txt
def write_output(outfile, max_score, seq1, seq2):
    with open(outfile, 'w') as out_f:
        out_f.write("score = %d\n" % max_score)
        out_f.write(">seq1\n%s\n\n" % seq1)
        out_f.write(">seq2\n%s" % seq2)


# function to read input.txt
def read_input(ifile):
    seq = {}
    try:
        f = open(ifile, 'r')
    except OSError:
        print("Error open ", ifile)
    else:
        tmp_line = ""
        for line in f.readlines():
            line = line.strip()
            if line.startswith(";"):
                continue
            if ";" in line:
                line = line.split(";")[0]
                line = line.strip()
            if line.startswith(">"):
                tmp_line = line[1:]
                seq[tmp_line] = ""
                continue
            if len(tmp_line) != 0:
                seq[tmp_line] = seq[tmp_line] + line
        f.close()
    return seq


# function to read parameter
def read_parameter(parameter_file):
    para = {"Init_gap": 0, "Base_Indel": 0, "alphabet": [], "matrix": {}}
    tmp_lines = []
    with open(parameter_file) as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith(";"):
                pass
            if ";" in line:
                line = line.split(";")[0]
                line = line.strip()
            if line:
                tmp_lines.append(line)

    para["Init_gap"] = int(tmp_lines[0])
    para["Base_Indel"] = int(tmp_lines[1])
    #    print("Gap penalty:", para["Init_gap"], para["Base_Indel"])

    para["alphabet"] = tmp_lines[2].split(" ")
    #    print("Alpha",para["alphabet"])

    alphabet_size = len(para["alphabet"])
    d = np.empty((alphabet_size, alphabet_size))

    matrix_lines = tmp_lines[3:]
    tmp_dir = {}
    for i in range(alphabet_size):
        t = matrix_lines[i].split(" ")
        while "" in t:
            t.remove("")
        for j in range(alphabet_size):
            tmp_dir[para["alphabet"][i] + para["alphabet"][j]] = int(t[j])
    para["matrix"] = tmp_dir

    return para


#### Write your functions below
def affine_gap_model(seq1, seq2, alphabets, score_matrix, gap_open_penalty, gap_extend_penalty):
    """
       Calculates the optimal global alignment between two sequences using the affine gap penalty model.

       Args:
           seq1 (str): The first input sequence.
           seq2 (str): The second input sequence.
           alphabets (str): A string of characters representing the alphabet used in the input sequences.
           score_matrix (ndarray): A ndarray of score matrix
           gap_open_penalty (int): The penalty for opening a gap in one of the sequences.
           gap_extend_penalty (int): The penalty for extending an existing gap in one of the sequences.

       Returns:
           aligned sequence 1, aligned sequence 2, alignment score

       """
    n = len(seq1)
    m = len(seq2)
    V = np.zeros((n + 1, m + 1))
    # E矩阵，用来存储gap在seq1中的值 F矩阵，用来存储gap在seq2中的值
    E = np.zeros((n + 1, m + 1))
    F = np.zeros((n + 1, m + 1))
    # 这样可以方便索引，可以直接根据字母判断是否匹配
    replace = {}
    for i, alphabet in enumerate(alphabets):
        replace[alphabet] = i
    # 初始化矩阵
    for i in range(n + 1):
        E[i][0] = -np.inf
        F[i][0] = gap_open_penalty + i * gap_extend_penalty
        V[i][0] = gap_open_penalty + i * gap_extend_penalty
    for j in range(m + 1):
        E[0][j] = gap_open_penalty + j * gap_extend_penalty
        F[0][j] = -np.inf
        V[0][j] = gap_open_penalty + j * gap_extend_penalty
    E[0][0] = -np.inf
    F[0][0] = -np.inf
    V[0][0] = 0
    # 填矩阵
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            E[i][j] = max(E[i][j - 1] + gap_extend_penalty, V[i][j - 1] + gap_open_penalty + gap_extend_penalty)
            F[i][j] = max(F[i - 1][j] + gap_extend_penalty, V[i - 1][j] + gap_open_penalty + gap_extend_penalty)
            match = V[i - 1][j - 1] + score_matrix[replace[seq1[i - 1]], [replace[seq2[j - 1]]]]
            delete = E[i][j]
            insert = F[i][j]
            V[i][j] = max(match, delete, insert)
    score = V[n][m]
    # 回溯
    aligned_seq_1 = ""
    aligned_seq_2 = ""
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and V[i][j] == V[i - 1][j - 1] + score_matrix[replace[seq1[i - 1]], [replace[seq2[j - 1]]]]:
            aligned_seq_1 = seq1[i - 1] + aligned_seq_1
            aligned_seq_2 = seq2[j - 1] + aligned_seq_2
            i -= 1
            j -= 1
        elif i > 0 and V[i][j] == F[i][j]:
            aligned_seq_1 = seq1[i - 1] + aligned_seq_1
            aligned_seq_2 = "_" + aligned_seq_2
            i -= 1
        elif j > 0 and V[i][j] == E[i][j]:
            aligned_seq_1 = "_" + aligned_seq_1
            aligned_seq_2 = seq2[j - 1] + aligned_seq_2
            j -= 1
    return aligned_seq_1, aligned_seq_2, score

def process_input(in_seq, para):
    '''
    The function process_input takes two arguments: in_seq, a dictionary containing two protein sequences, and para, a dictionary containing parameters such as alphabet, score matrix, gap penalties, etc. The function returns a tuple containing seq1, seq2, alphabets, score_matrix, gap_open_penalty, and gap_extend_penalty.

    Args:
        in_seq (dict): A dictionary containing two protein sequences with keys "seq1" and "seq2".
        para (dict): A dictionary containing the following keys:
            - "alphabet" (list): A list of characters representing the alphabet of the protein sequences.
            - "matrix" (dict): A dictionary containing scores for all possible pairs of characters in the alphabet.
            - "Init_gap" (float): Gap opening penalty.
            - "Base_Indel" (float): Gap extension penalty.

    Returns:
        seq1, seq2, alphabets, score_matrix, gap_open_penalty, gap_extend_penalty
    '''
    alphabets = para["alphabet"]
    score_dict = para["matrix"]
    score_matrix = np.zeros((len(alphabets), len(alphabets)))
    for i, a1 in enumerate(alphabets):
        for j, a2 in enumerate(alphabets):
            key = a1 + a2
            score_matrix[i][j] = score_matrix[j][i] = score_dict[key]
    alphabets = ''.join([str(elem) for elem in alphabets])
    seq1 = in_seq['seq1']
    seq2 = in_seq['seq2']
    gap_open_penalty = para["Init_gap"]
    gap_extend_penalty = para["Base_Indel"]
    return seq1, seq2, alphabets, score_matrix, gap_open_penalty, gap_extend_penalty


#### End of your functions


##### Main program

if len(sys.argv) != 4:
    print("python " + sys.argv[0] + " parameter.txt input.txt output.txt")
    sys.exit(0)
start_time = time.time()

parameter_file = sys.argv[1]
input_file = sys.argv[2]
out_file = sys.argv[3]
in_seq = read_input(input_file)
para = read_parameter(parameter_file)

### write your main script below

# The input is in_seq, para["Init_gap"], para["Base_Indel"], para["alphabet"], para["matrix"].
# You need to compute the global alignment
# You can comment the following five print statements when you start programming.
# print("Input sequences: ", in_seq)
# print("Init gap penalty: ", para["Init_gap"])
# print("Base extension penalty: ", para["Base_Indel"])
# print("Alphabet: ", para["alphabet"])
# print("Score matrix: ", para["matrix"])

seq1, seq2, alphabets, score_matrix, gap_open_penalty, gap_extend_penalty = process_input(in_seq, para)
align1, align2, max_score = affine_gap_model(seq1, seq2, alphabets, score_matrix, gap_open_penalty, gap_extend_penalty)
write_output(out_file, max_score, align1, align2)
### End of your main script

end_time = time.time()
print("time: ", end_time - start_time)

sys.exit()
