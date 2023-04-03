# _*_ coding: utf-8 _*_
# @Time : 2023-04-03 15:59 
# @Author : YingHao Zhang(池塘春草梦)
# @Version：v0.1
# @File : Affine_gap_model_global_alignment_NW.py
# @desc : Global alignment with affine gap model
'''
这个只考虑了匹配不匹配的情况，没有考虑更复杂的得分矩阵，而且也只有一条回溯路径。
'''

import numpy as np

# 要比对的两条序列、匹配得分、不匹配得分、空格起始罚分、空格延申罚分、序列包含的字母表
seq1 = 'ACCGA'
seq2 = 'AGTTA'
match_score = 1
mismatch_score = -3
gap_open_penalty = -1
gap_extend_penalty = -1
alphabets = "ACGT"


def affine_gap_model(seq1, seq2, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty, alphabets):
    """
       Calculates the optimal global alignment between two sequences using the affine gap penalty model.

       Args:
           seq1 (str): The first input sequence.
           seq2 (str): The second input sequence.
           match_score (int): The score for a match between two characters in the same position.
           mismatch_score (int): The score for a mismatch between two characters in the same position.
           gap_open_penalty (int): The penalty for opening a gap in one of the sequences.
           gap_extend_penalty (int): The penalty for extending an existing gap in one of the sequences.
           alphabets (str): A string of characters representing the alphabet used in the input sequences.

       Returns:
           aligned sequence 1, aligned sequence 2, alignment score

       Examples:
           # >>> affine_gap_model('ACCGA', 'AGTTA', 1, -3, -1, -1, "ACGT")
           ACCG--A A--GTTA -3.0
       """
    n = len(seq1)
    m = len(seq2)
    V = np.zeros((n + 1, m + 1))
    # E矩阵，用来存储gap在seq1中的值 F矩阵，用来存储gap在seq2中的值
    E = np.zeros((n + 1, m + 1))
    F = np.zeros((n + 1, m + 1))
    # 创建得分矩阵
    score_matrix = np.full((len(alphabets), len(alphabets)), mismatch_score)
    np.fill_diagonal(score_matrix, match_score)
    score_matrix = np.array(score_matrix)
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
            aligned_seq_2 = "-" + aligned_seq_2
            i -= 1
        elif j > 0 and V[i][j] == E[i][j]:
            aligned_seq_1 = "-" + aligned_seq_1
            aligned_seq_2 = seq2[j - 1] + aligned_seq_2
            j -= 1
    return aligned_seq_1, aligned_seq_2, score


align1, align2, score = affine_gap_model(seq1, seq2, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty,
                                         alphabets)
print("Global Alignmet using affine gap model between {} and {}:\n{}\n{}\nAlignment score is {}.".format(seq1, seq2,
                                                                                                         align1, align2,
                                                                                                         score))
