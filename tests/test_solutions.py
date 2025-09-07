# -- IMPORTS --

# -- Standard libraries --
import collections
import decimal
import os
import sys
import typing

from collections import Counter
from decimal import Decimal
decimal.getcontext().prec = 10
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
    assert fibo_rabbits(5, k=1) == 5
    assert fibo_rabbits(5, k=3) == 19


def test_max_gc_content():
    max_gc_seq_id, max_gc_content_ = max_gc_content(SeqIO.parse("tests/rosalind_gc.txt", "fasta"))
    assert max_gc_seq_id == 'Rosalind_6344'
    assert round(max_gc_content_, 5) == Decimal('51.68884')


def test_point_mutations():
    assert point_mutations('ACCGGGTTTT', 'ACCGGGTTTT') == 0
    assert point_mutations("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT") == 7


def test_transition_transversion_ratio():
    ttr = transition_transversion_ratio("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT")
    assert round(ttr, 5) == Decimal('0.16667')

    ttr = transition_transversion_ratio(
        "GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT",
        "TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT"
    )
    assert round(ttr, 5) == Decimal('1.21429')


def test_edit_distance():
    assert edit_distance("ACGT", "AGCT") == Decimal('2')
    assert edit_distance("AAGACTCTGG", "CGTTTAACTT") == Decimal('8')
    assert edit_distance("ACGT", "ACGT") == Decimal('0')
    assert edit_distance("ACGT", "") == Decimal('4')
    assert edit_distance("", "ACGT") == Decimal('4')


def test_translate_rna_into_protein():
    assert translate_rna_to_protein("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA") == 'MAMAPRTEINSTRING'


def test_splice_rna():
    assert splice_rna(
        "ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG",
        ("ATCGGTCGAA", "ATCGGTCGAGCGTGT")
    ) == "MVYIADKQHVASREAYGHMFKVCA"


def test_count_dna_motif():
    assert count_dna_motif("GATATATGCATATACTT", "TCAT") == ()
    assert count_dna_motif("GATATATGCATATACTT", "ATAT") == (2, 4, 10)
    assert count_dna_motif("", "AB") == ()
    assert count_dna_motif("ABC", "") == ()


def test_find_spliced_motif():
    assert find_spliced_motif("ACGTACGTGACG", "GTA") == (3, 4, 5)
    assert find_spliced_motif("ACGTACGTGACG", "AGTCC") == (1, 3, 4, 6, 11)
    assert find_spliced_motif("ACGTACGTGACG", "GTX") == ()
    assert find_spliced_motif("", "AB") == ()
    assert find_spliced_motif("ABC", "") == ()


def test_protein_mass():
    assert round(protein_mass("SKADYEK"), 5) == Decimal('821.39192')


def test_longest_common_shared_motif():
    assert longest_common_shared_motif(("ACCC", "GGTT")) == ""
    assert longest_common_shared_motif(("GATTACA", "TAGACCA", "ATACA")) == "TA"


def test_sequence_distance_matrix():
    assert sequence_distance_matrix(("TTTCCATTTA", "GATTCATTTC", "TTTCCATTTT", "GTTCCATTTA")) == (
        [[0.0, Decimal('0.4'), Decimal('0.1'), Decimal('0.1')],
         [Decimal('0.4'), 0.0, Decimal('0.4'), Decimal('0.3')],
         [Decimal('0.1'), Decimal('0.4'), 0.0, Decimal('0.2')],
         [Decimal('0.1'), Decimal('0.3'), Decimal('0.2'), 0.0]]
    )


def test_consensus_string():
    assert consensus_string(('ATCCAGCT', 'GGGCAACT', 'ATGGATCT', 'AAGCAACC', 'TTGGAACT', 'ATGCCATT', 'ATGGCACT')) == (
        'ATGCAACT',
        ([5, 1, 0, 0, 5, 5, 0, 0],
         [0, 0, 1, 4, 2, 0, 6, 1],
         [1, 1, 6, 3, 0, 1, 0, 0],
         [1, 5, 0, 0, 0, 1, 1, 6])
    )


def test_overlap_graph():
    seqs = ["AAATAAA", "AAATTTT", "TTTTCCC", "AAATCCC", "GGGTGGG"]
    sids = ["Rosalind_0498", "Rosalind_2391", "Rosalind_2323", "Rosalind_0442", "Rosalind_5013"]
    srecs = [Bio.SeqRecord.SeqRecord(seq, id=sid) for seq, sid in zip(seqs, sids)]
    assert overlap_graph(srecs, 3) == (
        ('Rosalind_0498', 'Rosalind_2391'),
        ('Rosalind_0498', 'Rosalind_0442'),
        ('Rosalind_2391', 'Rosalind_2323')
    )


def test_oriented_gene_orderings():
    assert list(oriented_gene_orderings(2)) == [
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
    assert linguistic_sequence_complexity("ATTTGGATT", "ACGT") == 0.875


def test_random_dna_strings():
    assert random_dna_strings(
        "ACGATACAA",
        (0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783),
        roundto=3
    ) == (Decimal('-5.737'), Decimal('-5.217'), Decimal('-5.263'), Decimal('-5.360'), Decimal('-5.958'), Decimal('-6.628'), Decimal('-7.009'))
