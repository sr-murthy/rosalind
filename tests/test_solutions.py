# -- IMPORTS --

# -- Standard libraries --
import collections
import os
import sys
import typing

from collections import Counter
from itertools import chain, permutations, product

# -- 3rd party libraries --
import Bio
import pytest

from Bio import SeqIO

# -- Internal libraries --import os, sys

sys.path.insert(0, 'src')
from solutions import *

# NOTE: Test cases are based mainly on the examples provided
#       in the ROSALIND problem pages, although there are
#       some which differ.


def test_basecount():
    res = basecount("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC")
    assert isinstance(res, collections.Counter)
    assert res['A'] == 20
    assert res['C'] == 12
    assert res['G'] == 17
    assert res['T'] == 21


def test_transcribe_dna_to_rna():
    assert transcribe_dna_to_rna("GATGGAACTTGACTACGTAAATT") == 'GAUGGAACUUGACUACGUAAAUU'


def test_reverse_complement():
    assert reverse_complement("AAAACCCGGT") ==  'ACCGGGTTTT'


def test_fibo_rabbits():
    assert fibo_rabbits(n=5, k=1) == 5
    assert fibo_rabbits(n=5, k=3) == 19


def test_max_gc_content():
    assert max_gc_content(SeqIO.parse("tests/rosalind_gc.txt", "fasta")) == ('Rosalind_6127', 26.282722513089006)


def test_point_mutations():
    assert point_mutations('ACCGGGTTTT', 'ACCGGGTTTT') == 0
    assert point_mutations("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT") == 7


def test_translate_rna_into_protein():
    assert translate_rna_to_protein("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA") == 'MAMAPRTEINSTRING'


def count_dna_motif():
    assert count_dna_motif("GATATATGCATATACTT", "TCAT") == ()
    assert count_dna_motif("GATATATGCATATACTT", "ATAT") == (2, 4, 10)


def test_protein_mass():
    assert protein_mass("SKADYEK") == 821.39192


def test_longest_common_substring():
    assert longest_common_substring(["ACCC", "GGTT"]) == ""
    assert longest_common_substring(["GATTACA", "TAGACCA", "ATACA"]) == "TA"


def test_sequence_distance_matrix():
    assert sequence_distance_matrix(["TTTCCATTTA", "GATTCATTTC", "TTTCCATTTT", "GTTCCATTTA"]) == (
        [[0.0, 0.4, 0.1, 0.1],
         [0.4, 0.0, 0.4, 0.3],
         [0.1, 0.4, 0.0, 0.2],
         [0.1, 0.3, 0.2, 0.0]]
    )


def test_profile_matrix():
    assert profile_matrix(['ATCCAGCT', 'GGGCAACT', 'ATGGATCT', 'AAGCAACC', 'TTGGAACT', 'ATGCCATT', 'ATGGCACT']) == (
        'ATGCAACT',
        [[5, 1, 0, 0, 5, 5, 0, 0],
         [0, 0, 1, 4, 2, 0, 6, 1],
         [1, 1, 6, 3, 0, 1, 0, 0],
         [1, 5, 0, 0, 0, 1, 1, 6]]
    )


def test_signed_permutations():
    assert list(signed_permutations(2)) == [
        (1, 2), (2, 1), (1, -2), (-2, 1), (-1, 2), (2, -1), (-1, -2), (-2, -1)
    ]


def test_lexicographic_kmers():
    assert list(lexicographic_kmers('ACGT', k=2)) == [
        'AA',
        'AC',
        'AG',
        'AT',
        'CA',
        'CC',
        'CG',
        'CT',
        'GA',
        'GC',
        'GG',
        'GT',
        'TA',
        'TC',
        'TG',
        'TT'
    ]


def test_kmer_composition():
    assert tuple(kmer_composition('CTTCGAAAG', 'ACGT', 2)) == (
        2, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1
    )


def test_variable_length_lexicographic_ordering():
    assert list(variable_length_lexicographic_ordering('DNA', k=2)) == [
        'D',
        'DD',
        'DN',
        'DA',
        'N',
        'ND',
        'NN',
        'NA',
        'A',
        'AD',
        'AN',
        'AA'
    ]


def test_linguistic_sequence_complexity():
    assert linguistic_sequence_complexity("ATTTGGATT", {"A", "C", "G", "T"}) == 0.875
