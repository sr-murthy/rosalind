__all__ = [
    'basecount',
    'consensus_string',
    'count_dna_motif',
    'edit_distance',
    'fibo_rabbits',
    'find_spliced_motif',
    'kmer_composition',
    'lexicographic_kmers',
    'linguistic_sequence_complexity',
    'longest_common_shared_motif',
    'max_gc_content',
    'MONOISOTOPIC_MASS_TABLE',
    'oriented_gene_orderings',
    'overlap_graph',
    'point_mutations',
    'protein_mass',
    'random_dna_strings',
    'reverse_complement',
    'RNA_CODON_TABLE',
    'sequence_distance_matrix',
    'splice_rna',
    'transcribe_dna_to_rna',
    'transition_transversion_ratio',
    'translate_rna_to_protein',
    'variable_length_lexicographic_ordering',
]


# -- IMPORTS --

# -- Standard libraries --
import collections
import decimal
import functools
import math
import typing

from collections import Counter
from decimal import Decimal
from itertools import permutations

# -- 3rd party libraries --
import Bio

from Bio import Seq, SeqIO


# -- Internal libraries --
from utils import (
    find_subsequence,
    find_substring,
    hamming_difference,
    hamming_distance,
    levenshtein_distance,
    longest_common_substring,
    signed_permutations,
    word_grams,
    word_k_grams,
)


#: From ROSALIND (protein mass problem) see the page below:
#:
#:  https://rosalind.info/glossary/monoisotopic-mass-table/
MONOISOTOPIC_MASS_TABLE = {
    "A":   71.03711,
    "C":   103.00919,
    "D":   115.02694,
    "E":   129.04259,
    "F":   147.06841,
    "G":   57.02146,
    "H":   137.05891,
    "I":   113.08406,
    "K":   128.09496,
    "L":   113.08406,
    "M":   131.04049,
    "N":   114.04293,
    "P":   97.05276,
    "Q":   128.05858,
    "R":   156.10111,
    "S":   87.03203,
    "T":   101.04768,
    "V":   99.06841,
    "W":   186.07931,
    "Y":   163.06333 
}


#: From ROSALIND (protein translation problem): see the page below:
#:
#:  https://rosalind.info/glossary/rna-codon-table/
RNA_CODON_TABLE = {
    'UUU': 'F',
    'CUU': 'L',
    'AUU': 'I',
    'GUU': 'V',
    'UUC': 'F',
    'CUC': 'L',
    'AUC': 'I',
    'GUC': 'V',
    'UUA': 'L',
    'CUA': 'L',
    'AUA': 'I',
    'GUA': 'V',
    'UUG': 'L',
    'CUG': 'L',
    'AUG': 'M',
    'GUG': 'V',
    'UCU': 'S',
    'CCU': 'P',
    'ACU': 'T',
    'GCU': 'A',
    'UCC': 'S',
    'CCC': 'P',
    'ACC': 'T',
    'GCC': 'A',
    'UCA': 'S',
    'CCA': 'P',
    'ACA': 'T',
    'GCA': 'A',
    'UCG': 'S',
    'CCG': 'P',
    'ACG': 'T',
    'GCG': 'A',
    'UAU': 'Y',
    'CAU': 'H',
    'AAU': 'N',
    'GAU': 'D',
    'UAC': 'Y',
    'CAC': 'H',
    'AAC': 'N',
    'GAC': 'D',
    'UAA': '',
    'CAA': 'Q',
    'AAA': 'K',
    'GAA': 'E',
    'UAG': '',
    'CAG': 'Q',
    'AAG': 'K',
    'GAG': 'E',
    'UGU': 'C',
    'CGU': 'R',
    'AGU': 'S',
    'GGU': 'G',
    'UGC': 'C',
    'CGC': 'R',
    'AGC': 'S',
    'GGC': 'G',
    'UGA': '',
    'CGA': 'R',
    'AGA': 'R',
    'GGA': 'G',
    'UGG': 'W',
    'CGG': 'R',
    'AGG': 'R',
    'GGG': 'G'
}


@functools.cache
def basecount(s: str | Bio.Seq.Seq, /) -> collections.Counter:
    """:py:class:`collections.Counter` : Returns a dict of DNA bases and their counts in a given DNA sequence.

    Solution to the Counting DNA Nucleotides problem (DNA):

    https://rosalind.info/problems/dna/

    Parameters
    ----------
    s : str, Bio.Seq.Seq
        The input DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    dict
        A dict of bases and their integer counts in the sequence.

    Examples
    --------
    >>> basecount("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC")
    Counter({'T': 21, 'A': 20, 'G': 17, 'C': 12})
    """
    return Counter(s)


@functools.cache
def transcribe_dna_to_rna(s: str | Bio.Seq.Seq, /) -> str:
    """:py:class:`str` : Returns an RNA transcription of a DNA string.

    Solution to the Transcribing DNA into RNA problem (RNA):

    https://rosalind.info/problems/rna/

    Parameters
    ----------
    s : str
        The input DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    A string of the RNA sequence transcribed from the given DNA sequence.

    Examples
    --------
    >>> transcribe_dna_to_rna("GATGGAACTTGACTACGTAAATT")
    'GAUGGAACUUGACUACGUAAAUU'
    """
    return str(s).translate(str.maketrans('T', 'U'))


@functools.cache
def reverse_complement(s: str | Bio.Seq.Seq, /) -> str:
    """:py:class:`str` : Returns the reverse complement of a DNA string, which is the reversed string with the bases complemented.

    Solution to the Complementing a Strand of DNA problem (REVC):

    https://rosalind.info/problems/revc/

    Parameters
    ----------
    s : str
        The input DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    The reverse complement of an input DNA sequence: this is the original DNA
    sequence reversed with complements taken of the nucleobases:
    ::
        A -> T
        T -> A
        C -> G
        G -> C

    Examples
    --------
    >>> reverse_complement("AAAACCCGGT")
    'ACCGGGTTTT'
    """
    return str(s)[::-1].translate(str.maketrans("ATCG", "TAGC"))


@functools.cache
def fibo_rabbits(n: int, k: int) -> int:
    """:py:class:`int` : The number of rabbit pairs alive after ``n`` months, starting with 1 pair and with each reproductive-age pair producing ``k`` pairs.

    Solution to the Rabbits and Recurrence Relations problem (FIB):

    https://rosalind.info/problems/fib/

    .. note::

       In month ``n`` the number of rabbit pairs present is equal to the
       number of the previous month's pairs plus any new pairs (offspring
       pairs). We are told that each reproductive-age rabbit pair produces
       ``k`` pairs. The number of new pairs produced in month ``n`` is equal
       to the number that were alive two months prior. If ``f`` is the
       function then this gives the formula:
       ::

            f(n, k) = f(n - 1, k) + 3 * f(n - 2, k)

    This solution is iterative (and non-recursive), to avoid recursion errors
    caused by large values of ``n``.

    Parameters
    ----------
    n  : int
        The number of months elapsed.

    k : int
        The number of pairs produced by each reproduction-age rabbit pair.

    Returns
    -------
    int
        The number of rabbit pairs alive after ``n`` months, starting with 1
        pair, and with each reproductive-age pair producing ``k`` pairs.

    Examples
    --------
    >>> fibo_rabbits(5, 3)
    19
    """
    # An inner Fibo generator depending on ``k`` only, that can
    # generate infinitely.
    def fibo_inf(k: int) -> typing.Generator[int, None, None]:
        a, b = 0, 1
        yield a
        yield b

        while True:
            a, b = b, k * a + b
            yield b

    for i, x in enumerate(fibo_inf(k)):
        if i == n:
            return x


def max_gc_content(fasta_records: Bio.SeqIO.FastaIO.FastaIterator | typing.Iterable[Bio.SeqRecord.SeqRecord], /) -> tuple[str, decimal.Decimal]:
    """:py:class:`str` : Returns a tuple containing a record ID and GC content percentage for the FASTA record with the highest GC content in a FASTA file.

    Solution to the Computing GC Content problem (GC):

    https://rosalind.info/problems/gc/

    Parameters
    ----------
    fasta_records : Bio.SeqIO.FastaIO.FastaIterator, typing.Iterable
        Either an iterator of (``biopython``) FASTA records
        (:py:class:`Bio.SeqIO.FastaIO.FastaIterator`), or an iterable
        of such records (:py:class:`Bio.SeqRecord.SeqRecord`).

    Returns
    -------
    tuple
        A tuple containing a record ID and GC content percentage for the FASTA
        record with the highest GC content in the given FASTA file.

    Examples
    --------
    >>> from Bio import SeqIO
    >>> max_gc_content(SeqIO.parse("./rosalind_gc.txt", "fasta"))
    ('Rosalind_6344', Decimal('51.68884339815762'))
    """
    def gc_content(s: str | Bio.Seq.Seq, /) -> decimal.Decimal:
        return Decimal(sum(1 for base in s if base in ['C', 'G'])) / Decimal(len(s))

    sorted_gc_contents: dict[str, Decimal] = sorted(
        {
            record.id: gc_content(record.seq) * 100
            for record in fasta_records
        }.items(),
        key=lambda t: t[1]
    )

    return sorted_gc_contents[-1]


@functools.cache
def point_mutations(s: str | Bio.Seq.Seq, t: str | Bio.Seq.Seq, /) -> int:
    """:py:class:`int` : Returns an integer count of the mutations in two equal-length DNA sequences.

    Solution to the Counting Point Mutations problem (HAMM):

    https://rosalind.info/problems/hamm/

    This is an application of the **Hamming distance** metric.

    Parameters
    ----------
    s : str, Bio.Seq.Seq
        The first DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    t : str, Bio.Seq.Seq
        The second DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    int
        A (non-negative) integer count of the mutations in the two given DNA
        sequences.

    Examples
    --------
    >>> point_mutations("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT")
    7
    """
    return hamming_distance(s, t)


@functools.cache
def transition_transversion_ratio(s: str | Bio.Seq.Seq, t: str | Bio.Seq.Seq, /) -> decimal.Decimal:
    """:py:class:`decimal.Decimal`` : The ratio of transitions to transversions in two equal-length DNA sequences.

    Solution to the Transitions and Traversions problem (TRAN):

    https://rosalind.info/problems/tran/

    In two equal-length DNA sequences transitions and transversions are defined as
    follows:

        transition:  A <-> G, C <-> T
        tranversion: A <-> C, A <-> T, G <-> C, G <-> T

    This means that transitions and transversions are mutually exclusive: a
    change of base is either a transition or a transversion.

    Parameters
    ----------
    s : str, Bio.Seq.Seq
        The first DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    t : str, Bio.Seq.Seq
        The second DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    decimal.Decimal
        The ratio of transitions to transversions in the two strings.

    Examples
    --------
    >>> transition_transversion_ratio("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT")
    """
    transitions: int = 0
    transversions: int = 0

    for _, diff in hamming_difference(s, t):
        if diff in [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]:
            transitions += 1
        else:
            transversions += 1

    return Decimal(transitions) / Decimal(transversions)


@functools.cache
def edit_distance(s: str, t: str, /) -> int:
    """:py:class:`int` : Returns the edit distance (more generally called Levenshtein distance) between two DNA/RNA/protein strings/sequences.

    This is a solution to the Edit Distance problem (EDIT):

    https://rosalind.info/problems/edit/

    Calls on :py:func:`utils.levenshtein_distance` with insertion, deletion and
    substitution costs set to 1.

    .. note::

       This cached recursive implementation is slower than equivalent iterative
       implementations, but is more readable and thus easier to understand.

    Parameters
    ----------
    s : str
        The first string.

    t : str
        The second string, equal in length to the first.

    Returns
    -------
    int
        The edit distance between the two strings.

    Examples
    --------
    >>> edit_distance("ACGT", "AGCT")
    2
    >>> edit_distance("AAGACTCTGG", "CGTTTAACTT")
    8
    >>> edit_distance("ACGT", "ACGT")
    0
    >>> edit_distance("ACGT", "")
    4
    >>> edit_distance("", "ACGT")
    4
    """
    return levenshtein_distance(s, t, insertion_cost=1, deletion_cost=1, substitution_cost=1)


@functools.cache
def translate_rna_to_protein(s: str | Bio.Seq.Seq, /) -> str:
    """:py:class:`str` : Returns a protein sequence encoded by an RNA seq.

    Solution to the Translating RNA into Protein problem (PROT):

    https://rosalind.info/problems/prot/

    .. note::

       An assumption is made that the length of the RNA sequence is a multiple
       of 3, as codons in the RNA codon table are all of length 3.

    Parameters
    ----------
    s : str, Bio.Seq.Seq
        The RNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    str
        The protein sequence encoded by the RNA sequence.

    Examples
    --------
    >>> translate_rna_to_protein("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")
    'MAMAPRTEINSTRING'
    """
    encoding: str = ''
    i: int = 0

    while i < len(s) - 2:
        codon: str = s[i: i + 3]
        encoding += RNA_CODON_TABLE[codon]
        i += 3

    return encoding


@functools.cache
def splice_rna(s: str | Bio.Seq.Seq, introns: tuple[str | Bio.Seq.Seq], /) -> str:
    """:py:class:`str` : Solution to the RNA Splicing problem.

    Solution to the RNA Splicing problem (SPLC):

    https://rosalind.info/problems/splc/

    Parameters
    ----------
    s : str
        The input DNA string/sequence.

    introns : tuple
        A tuple of segments of the input string which represent entities
        called introns, which must be removed.  The tuple requirement is due
        to the fact that the function uses a cache, which requires arguments
        to be hashable.

    Returns
    -------
    str
        The translated protein string.

    Examples
    --------
    >>> s = "ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG"
    >>> splice_rna(s, ["ATCGGTCGAA", "ATCGGTCGAGCGTGT"])
    'MVYIADKQHVASREAYGHMFKVCA'
    """
    # A local initial copy of ``s``
    s_copy: str = s

    # Remove the introns - note that using ``re.sub`` won't work if the inputs are
    # ``Bio.Seq.Seq`` objects, hence the use of ``str.replace``.
    for i in introns:
        s_copy = s_copy.replace(i, '')

    # Now perform the DNA -> RNA -> Protein transformation and return
    return translate_rna_to_protein(transcribe_dna_to_rna(s_copy))


@functools.cache
def count_dna_motif(s: str | Bio.Seq.Seq, t: str | Bio.Seq.Seq, /) -> tuple[int]:
    """:py:class:`tuple` : A tuple of all starting indices of occurrences of a DNA subsequence in the given DNA sequence.

    Solution to the Finding a Motif in DNA problem (SUBS);

    https://rosalind.info/problems/subs/

    .. note::

       The answer is given in terms of indices of 1-indexed arrays - to get
       zero indices just subtract 1.

    Parameters
    ----------
    s : str, Bio.Seq.Seq
        The DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    t : str, Bio.Seq.Seq
        The DNA subsequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    tuple
        Tuple of starting positions of all occurrences of the DNA subsequence
        in the given DNA sequence.

    Examples
    --------
    >>> count_dna_motif("GATATATGCATATACTT", "ATAT")
    (2, 4, 10)
    """
    base: tuple[int] | typing.Literal[()] = find_substring(s, t)

    return tuple(map(lambda i: i + 1, base)) if base else ()


@functools.cache
def find_spliced_motif(s: str | Bio.Seq.Seq, t: str | Bio.Seq.Seq, /) -> tuple[int] | typing.Literal[()]:
    """:py:class:`tuple` or None : Returns a tuple of 1-indexed array indices of a sequence ``t`` if it is a subsequence of ``s``, or null if not.

    Solution to the Finding a Spliced Motif problem (SSEQ):

    https://rosalind.info/problems/sseq/

    .. note::

       The indices are given in terms of 1-indexed arrays - just subtract 1
       from them all to convert them to 0-indexed array indices. Also, the
       indices will not necessarily be unique, as there may be multiple
       occurrences - the earliest occurring indices will be returned,
       in case of a subsequence match, and later occurrences, if they exist,
       are ignored.

    Parameters
    ----------
    s : str, Bio.Seq.Seq
        The DNA sequence in which to search for another sequence.

    t : str, Bio.Seq.Seq
        The second DNA sequence to search for in the first sequence.

    Returns
    -------
    tuple
        The 1-indexed ``s``-indices of the bases of ``t``, if ``t`` occurs as a
        subsequence of ``s``, or null if not. Not necessarily unique - the
        earliest occurring indices are returned in case of a subsequence match.

    Examples
    --------
    >>> find_spliced_motif("ACGTACGTGACG", "GTA")
    (3, 4, 5)
    """
    ixs: tuple[int] | typing.Literal[()] = find_subsequence(s, t)

    return tuple(map(lambda i: i + 1, ixs)) if ixs else ()


@functools.cache
def protein_mass(s: str | Bio.Seq.Seq, /) -> decimal.Decimal:
    """:py:class:`decimal.Decimal`` : Calculates the mass of a protein sequence using the monoisotopic mass table.

    Solution to the Calculating Protein Mass problem (PRTM):

    https://rosalind.info/problems/prtm/

    Parameters
    ----------
    s : str, Bio.Seq.Seq
        The protein sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    decimal.Decimal
        The protein string mass.

    Examples
    --------
    >>> round(protein_mass("SKADYEK"), 5)
    Decimal('821.39192')
    """
    return Decimal(sum(MONOISOTOPIC_MASS_TABLE[char] for char in s))


@functools.cache
def longest_common_shared_motif(seqs: tuple[str | Bio.Seq.Seq], /) -> str:
    """:py:class:`str` : Returns a longest common substring among an iterable of strings.

    Solution to the Finding a Shared Motif problem (LCSM):

    https://rosalind.info/problems/lcsm/

    Parameters
    ----------
    seqs : tuple
        A tuple of strings, which may be given as plain strings or
        Biopython genetic sequences (:py:class:`Bio.Seq.Seq`).  The tuple
        requirement is due to the fact that the function uses a cache,
        which requires arguments to be hashable.

    Returns
    -------
    str
        A longest common substring (not necessarily unique) - the first
        one encountered is returned.

    Examples
    --------
    >>> seqs = ["GATTACA", "TAGACCA", "ATACA"]
    >>> longest_common_shared_motif(seqs)
    'TA'
    """
    return longest_common_substring(seqs)


@functools.cache
def sequence_distance_matrix(seqs: tuple[str | Bio.Seq.Seq], /) -> list[list[decimal.Decimal]]:
    """:py:class:`list` : Returns a (symmetric) matrix of relative Hamming distances for all (ordered) pairs of sequences from an iterable of equal-length genetic sequences.
    
    Solution to the Creating a Distance Matrix problem (PDST):

    https://rosalind.info/problems/pdst/

    As the distance between a sequence and itself is zero, and the distance
    between two sequences is independent of order, the matrix will be symmetric
    with zeros on the diagonal.

    Parameters
    ----------
    seqs : tuple
        An iterable of equal-length, plain Python strings or genetic sequences
        given as Biopython sequences (:py:class:`Bio.Seq.Seq`).  The tuple
        requirement is due to the fact that the function uses a cache, which
        requires arguments to be hashable.

    Returns
    -------
    list
        An ``m`` by ``m`` symmetric matrix of relative Hamming distances of all
        pairs of sequences from the given iterable of sequences, where ``m`` is
        the common (equal) length of the sequences.

    Examples
    --------
    # Define a tuple input for the function caching to work
    >>> seqs = ("TTTCCATTTA", "GATTCATTTC", "TTTCCATTTT", "GTTCCATTTA")
    >>> sequence_distance_matrix(seqs)
    [[0.0, Decimal('0.4'), Decimal('0.1'), Decimal('0.1')],
     [Decimal('0.4'), 0.0, Decimal('0.4'), Decimal('0.3')],
     [Decimal('0.1'), Decimal('0.4'), 0.0, Decimal('0.2')],
     [Decimal('0.1'), Decimal('0.3'), Decimal('0.2'), 0.0]]
    """
    n: int = len(seqs)

    # Set up ``n x n`` matrix of zeros - note that this construction below
    # using a list comprehension is designed to ensure that all of the zero
    # arrays are different objects in memory.
    mat: list = []
    [mat.append(list([0.] * n)) for i in range(n)]

    # A variable to store the common sequence length from the length of the
    # 1st sequence
    m: int = len(seqs[0])

    # The outer loop on rows
    i: int = 0
    while i < n:
        # The column index always starts from ``i`` so that we always have
        # ``i <= j``
        j = i
        # The inner loop on columns
        while j < n:
            # if row and column indices are equal, there's nothing to do, so
            # increment column index and continue to next iteration
            if i == j:
                j += 1
                continue
            # otherwise it must be that ``i < j``, in which case set
            # the Hamming distance of the current sequence pair ``(i, j)``,
            # and use the the same value for the transposed entry ``(j, i)``,
            # avoiding another calculation.
            else:
                mat[i][j] = Decimal(point_mutations(seqs[i], seqs[j])) / Decimal(m)
                mat[j][i] = mat[i][j]
            j += 1
        i += 1

    return mat


@functools.cache
def consensus_string(seqs: tuple[str | Bio.Seq.Seq]) -> tuple[str, list[int]]:
    """:py:class:`tuple` : Returns the consensus string and profile matrix for a collection of equal-lenth DNA sequences/strings.

    Solution to the Consensus and Profile problem (CONS):

    https://rosalind.info/problems/cons/

    .. note::

       The iterable must be hashable, e.g. tuples not lists: it must contain
       strings or ``Bio.Seq.Seq`` objects, which are hashable, and not
       ``Bio.SeqRecord.SeqRecord`` objects, which are not hashable.

       The profile matrix is returned as a tuple of integer-valued lists.

    Parameters
    ----------
    seqs : tuple
        A tuple of equal length DNA sequences (or strings). The tuple
        requirement is due to the fact that the function uses a cache, which
        requires arguments to be hashable.

    Returns
    -------
    tuple
        The consensus string and profile matrix for the DNA sequences.

    Examples
    --------
    >>> seqs = ('ATCCAGCT', 'GGGCAACT', 'ATGGATCT', 'AAGCAACC', 'TTGGAACT', 'ATGCCATT', 'ATGGCACT')
    >>> cs, pm = consensus_string(seqs)
    >>> print(cs)
    ATGCAACT
    >>> for i, b in enumerate('ACGT'):
    ...     print(f'{b}: ' + ' '.join(map(str, pm[i])))
    A: 5 1 0 0 5 5 0 0
    C: 0 0 1 4 2 0 6 1
    G: 1 1 6 3 0 1 0 0
    T: 1 5 0 0 0 1 1 6
    """
    # Define the width of the matrix.
    n: int = len(seqs[0])

    # Create the initial ``4 x n`` profile matrix of zeros
    profile_matrix: list[list[int]] = [list([0] * n) for i in range(4)]

    # A list of blanks to store the characters of the consensus string
    consensus_str: list[str] = list([''] * n)

    # The bases string (in order, and with no repetitions)
    base_str: str = 'ACGT'

    # The outer loop on columns (letters of the consensus string)
    for j in range(n):
        # The inner loop on rows (bases)
        for i, base in enumerate(base_str):
            # Set the ``(i, j)``-th value as the number of sequences whose
            # ``j``-th letter is equal to the current base ``b``
            profile_matrix[i][j] = sum(1 for seq in seqs if seq[j] == base)
            # Having completed the last row of the profile matrix for
            # column ``j`` calculate the ``j``-th letter of the consensus
            # string as the base with highest frequency in the column.
            if i == 3:
                consensus_str[j] = max(base_str, key=lambda x: profile_matrix[base_str.index(x)][j])
    
    return ''.join(consensus_str), tuple(profile_matrix)


def overlap_graph(seqs: typing.Iterable[Bio.SeqRecord.SeqRecord], k: int, /) -> tuple[tuple[str, str]]:
    """:py:class:`tuple` : Returns a tuple of edges of the ``O_k`` overlap graph of a tuple of DNA sequences/strings.

    Solution to the Overlap Graphs problem (GRPH):

    https://rosalind.info/problems/grph/

    The ``O_k`` overlap graph of a set of DNA sequences is defined as a
    directed graph where nodes/vertices are sequences, and there is
    an edge for every ordered pair of sequences such that the ``k``-length
    suffix of the first matches the ``k``-length prefix of the second. Pairs
    of identical sequences are excluded.

    .. note::
    Under this definition it is possible that the graph may include edges
    that connect the same pair of sequences in both directions, e.g. if
    ``s`` is ``"AAACGGG"``, ``t``is ``"GGGTAAA"``, and ``k`` is ``3``. This
    means that we must use ordered pairs of sequences, not unordered pairs.

    In this implementation, as it is a solution for the ROSALIND problem which
    requires the output of sequence pairs to be given in terms of sequence IDs,
    the nodes of the graph are sequence IDs.

    Parameters
    ----------
    seqs : typing.Iterable
        An iterable of DNA sequences/strings given as
        :py:class:`Bio.SeqRecord.SeqRecord` objects.

    k : int
        The graph parameter ``k``.

    Returns
    -------
    tuple
        A tuple of edges of the ``O_k`` overlap graph of the given sequences,
        where each edge is a pair of sequence IDs which satisfy the overlap
        graph requirement described above.
    """
    G: list[tuple[str, str]] = []

    for s, t in permutations(seqs, 2):
        # Note that the use of ``itertools.permutations`` ensures that pairs
        # involving an identical sequence are excluded, which means we don't
        # have to check if they are different. 
        if s.seq[-k:] == t.seq[:k]:
            G.append((s.id, t.id))

    return tuple(G)


def oriented_gene_orderings(n: int, /) -> typing.Generator[tuple[int], None, None]:
    """:py:class:`typing.Generator` : Solution to the SIGN problem.

    Solution to the Enumerating Oriented Gene Orderings problem (SIGN):

    https://rosalind.info/problems/sign/

    Parameters
    ----------
    n : int
        The (common) length of the signed permutations of ``1..n``.

    Yields
    ------
    tuple
        Each generated value is a tuple of ints.

    Examples
    --------
    >>> list(oriented_gene_orderings(3))
    [(1, 2, 3),
     (1, 3, 2),
     (2, 1, 3),
     (2, 3, 1),
     (3, 1, 2),
     (3, 2, 1),
     (1, 2, -3),
     (1, -3, 2),
     (2, 1, -3),
     (2, -3, 1),
     (-3, 1, 2),
     (-3, 2, 1),
     (1, -2, 3),
     (1, 3, -2),
     (-2, 1, 3),
     (-2, 3, 1),
     (3, 1, -2),
     (3, -2, 1),
     (1, -2, -3),
     (1, -3, -2),
     (-2, 1, -3),
     (-2, -3, 1),
     (-3, 1, -2),
     (-3, -2, 1),
     (-1, 2, 3),
     (-1, 3, 2),
     (2, -1, 3),
     (2, 3, -1),
     (3, -1, 2),
     (3, 2, -1),
     (-1, 2, -3),
     (-1, -3, 2),
     (2, -1, -3),
     (2, -3, -1),
     (-3, -1, 2),
     (-3, 2, -1),
     (-1, -2, 3),
     (-1, 3, -2),
     (-2, -1, 3),
     (-2, 3, -1),
     (3, -1, -2),
     (3, -2, -1),
     (-1, -2, -3),
     (-1, -3, -2),
     (-2, -1, -3),
     (-2, -3, -1),
     (-3, -1, -2),
     (-3, -2, -1)]
    """
    yield from signed_permutations(n)


def lexicographic_kmers(s: str, k: int) -> typing.Generator[str, None, None]:
    """:py:class:`typing.Generator` : Returns a generator of ``k``-length substrings formed from subsequences of a given sequence of characters.
    
    Solution to the Enumerating k-mers Lexicographically problem (LEXF):

    https://rosalind.info/problems/lexf/

    Parameters
    ----------
    s : typing.Sequence
        A sequence of characters, which can be any single letter Unicode
        strings. The given ordering is used in building lexicographically
        ordered substrings.

    k : int
        The length of the substrings of the sequence to be generated.

    Yields
    ------
    str
        ``k``-length substrings formed from the sequence generted in
        lexicographic order (with respect to the original sequence).

    Examples
    --------
    >>> list(lexicographic_kmers('ACGT', k=2))
    ['AA',
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
     'TT']
    """
    yield from word_k_grams(s, k=k)


def kmer_composition(s: str | Bio.Seq.Seq, A: str, k: int) -> typing.Generator[int, None, None]:
    """:py:class:`typing.Generator` : Generates the ``k``-mer composition of a string with respect to an (ordered) alphabet.

    Solution to the k-Mer Composition problem (KMER):

    https://rosalind.info/problems/kmer/

    .. note::

       The function has been implemented as a generator, because ``k``-mers
       can grow very rapidly depending on the size of the alphabet and the
       value of ``k``.

    Parameters
    ----------
    s : str
        The string/sequence for which the ``k``-mer composition is sought.

    A : str, list, tuple
        The (ordered) alphabet from which the string/sequence is formed.

    k : int
        The (fixed) length of the kmer

    Yields
    -------
    int
        A generator of ``k``-mer frequencies in the given string/sequence, where
        the ``i``th value is the ``i``-th ``k``-mer of ``A`` in
        ``s``, given in the natural ordering of all ``k``-mers of ``A``.

    Examples
    --------
    >>> tuple(kmer_composition('AAAACCCGGT', 'ACGT', k=2))
    (3, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0)
    """
    yield from (len(find_substring(s, kmer)) for kmer in word_k_grams(A, k=k))


def variable_length_lexicographic_ordering(s: str, k: int) -> typing.Generator[str, None, None]:
    """:py:class:`typing.Generator` : Returns a generator of lexicographically ordered, variable-length substrings ,of size at most ``k`` , of a given string / sequence.
    
    Solution to the Ordering Strings of Varying Length Lexicographically problem (LEXV):

    https://rosalind.info/problems/lexv/

    Parameters
    ----------
    s : typing.Sequence
        A sequence of characters, which can be any single letter Unicode
        strings. The given ordering is used in building lexicographically
        ordered substrings.

    Yields
    ------
    str
        Variable-length substrings of length at most ``k``, formed from
        the sequence, generated in lexicographic order (with respect to
        the original sequence).

    Examples
    --------
    >>> list(variable_length_lexicographic_ordering('DNA', k=3))
    ['D',
     'DD',
     'DDD',
     'DDN',
     'DDA',
     'DN',
     'DND',
     'DNN',
     'DNA',
     'DA',
     'DAD',
     'DAN',
     'DAA',
     'N',
     'ND',
     'NDD',
     'NDN',
     'NDA',
     'NN',
     'NND',
     'NNN',
     'NNA',
     'NA',
     'NAD',
     'NAN',
     'NAA',
     'A',
     'AD',
     'ADD',
     'ADN',
     'ADA',
     'AN',
     'AND',
     'ANN',
     'ANA',
     'AA',
     'AAD',
     'AAN',
     'AAA']
    """
    yield from word_grams(s, k=k)


@functools.cache
def linguistic_sequence_complexity(s: str, A: str, /) -> decimal.Decimal:
    """:py:class:`float` : Returns the linguistic complexity (LC) of a string ``s`` formed over an alphabet ``A``.

    Solution to the Linguistic Complexity of a Genome problem (LING):

    https://rosalind.info/problems/ling/

    The LC measure is defined on the ROSALIND problem page:

        https://rosalind.info/problems/ling/

    but a clearer definition is also in the following paper (page W630):

        Orlov YL, Potapov VN. Complexity: an internet resource for analysis
        of DNA sequence complexity. Nucleic Acids Res. 2004 Jul 1;32(Web Server
        issue):W628-33. doi: 10.1093/nar/gkh466. PMID: 15215465; PMCID:
        PMC441604.

    This implementation uses the LC formula from this paper, which is called
    CL.

    Parameters
    ----------
    s : str
        The input string.

    A : str
        The underlying alphabet of the string, as a string (unique letters
        of the alphabet in the natural order, without repetitions).

    Returns
    -------
    decimal.Decimal
        The linguistic complexity of the string.

    Examples
    --------
    >>> linguistic_sequence_complexity("ATTTGGATT", "ACGT")
    Decimal('0.875')
    """
    n: int = len(s)
    m: int = len(A)

    # Calculate the number of observed substrings as the sum of the number of
    # ``k``-length substrings found in the string, for ``k`` in ``1..n``.
    num_observed_substrings: int = len(set(s)) + sum(
        len(set([s[j: j + k] for j in range(n - k + 1)]))
        for k in range(2, n + 1)
    )

    # Calculate the number of possible substrings as the sum of the maximum
    # number of ``k``-length substrings formed from the alphabet that can fit
    # into the given string, overlapping as necessary. For each ``k`` the
    # maximum number of possible ``k``-length substrings that can fit into
    # the string is equal to ``min(m^k, n - k + 1)``.
    num_possible_substrings: int = sum(min(m ** k, n - k + 1) for k in range(1, n + 1))

    return Decimal(num_observed_substrings) / Decimal(num_possible_substrings)


@functools.cache
def random_dna_strings(s: str | Bio.Seq.Seq, A: tuple[float | decimal.Decimal], /, *, roundto: int = 3) -> tuple[decimal.Decimal]:
    """:py:class:`tuple` : An array of common logarithms of probabilities determined by a given string and an array of possible GC content values.

    Solution to the Introduction to Random Strings problem (PROB):

    https://rosalind.info/problems/prob/

    .. note::

       The array of GC content values must be given as tuple, because the
       function uses caching (:py:func:`functools.cache`), which requires
       the arguments to be hashable: in Python tuples are immutable and
       therefore hashable, while mutable containers such as lists, sets and
       dicts are not hashable.

    Returns an array :math:`B` having the same length as :math:`A` in which
    :math:`B[i]` represents the common logarithm (:math:`log_10`) of the
    probability that a random string :math:`s`
    constructed with the GC-content found in :math:`A[i]` will match
    :math:`s` exactly.

    If :math:`s` is :math:`s_0s_1 \\cdots s_{n - 1}` then the probability that
    a random string with a GC content :math:`x` will match :math:`s` exactly
    is given by:

    .. math::

       \\begin{align}
       \\prod_{i=0}^{n - 1} P(s_i) &= \\sum_{i=0}^n P(s_i) \\
                             &= (n_A + n_T)\\frac{x}{2} + (n_C + n_G)\\frac{(1 - x)}{2}
       \\end{align}

    where :math:`n_A, n_C, n_G, n_T` are the base counts of :math:`s`, and
    :math:`P(s_0),\\ldots,P(s_{n - 1})` are the probabilities of the letters
    of :math:`s`.

    Parameters
    ----------
    s : str
        The input DNA string (or sequence).

    A : tuple
        A tuple of floats or :py:class``decimal.Decimal`` objects, representing
        possible GC content values of the given string.  The tuple requirement
        is due to the fact that the function uses a cache
        (:py:func:`functools.cache`) which requires arguments to be hashable.

    roundto : int, default=3
        An optional number of digits to round the values to, with a default
        of ``3``.

    Returns
    -------
    tuple
        A tuple of common logarithms of probabilities, where the ``i``-th
        value represents the base 10 log of the probability that a random
        string constructed with the GC-content found in ``A[i]`` will match
        ``s`` exactly.

    Examples
    --------
    >>> random_dna_strings("ACGATACAA", (0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783))
    (Decimal('-5.737'), Decimal('-5.217'), Decimal('-5.263'), Decimal('-5.360'), Decimal('-5.958'), Decimal('-6.628'), Decimal('-7.009'))
    """
    # A function to build a frequency table for the bases based on the given
    # GC content.
    def base_frequency_table(gc_content: decimal.Decimal) -> dict[str, decimal.Decimal]:
        x: decimal.Decimal = Decimal(gc_content)

        return {
            'A': (1 - x) / 2,
            'C': x / 2,
            'G': x / 2,
            'T': (1 - x) / 2
        }

    logs: list[decimal.Decimal] = []

    for gc_content in A:
        table: dict[str, decimal.Decimal] = base_frequency_table(Decimal(gc_content))
        logs.append(round(Decimal(sum(math.log10(table[base]) for base in s)), roundto))

    return tuple(logs)
