__all__ = [
    'basecount',
    'count_dna_motif',
    'fibo_rabbits',
    'gc_content',
    'kmer_composition',
    'lexicographic_kmers',
    'linguistic_sequence_complexity',
    'longest_common_substring',
    'max_gc_content',
    'MONOISOTOPIC_MASS_TABLE',
    'point_mutations',
    'profile_matrix',
    'protein_mass',
    'reverse_complement',
    'RNA_CODON_TABLE',
    'sequence_distance_matrix',
    'signed_permutations',
    'transcribe_dna_to_rna',
    'translate_rna_to_protein',
    'variable_length_lexicographic_ordering',
]


# -- IMPORTS --

# -- Standard libraries --
import collections
import typing

from collections import Counter
from itertools import chain, permutations, product

# -- 3rd party libraries --
import Bio

from Bio import SeqIO

# -- Internal libraries --

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
        yield 0
        yield 1

        a = 0
        b = 1
        while True:
            a, b = b, k * a + b
            yield b

    for i, x in enumerate(fibo_inf(k)):
        if i == n:
            return x


def gc_content(s: str | Bio.Seq.Seq, /) -> float:
    """:py:class:`float` : The proportion of bases in the DNA sequence which are either ``'C'`` or ``'G'``.

    Utility function for the max. GC content problem (GC):

    .. note::

       The returned value is a proportion (or ratio), not a percentage.

    Parameters
    ----------
    s : str, Bio.Seq.Seq
        The DNA sequence.

    Returns
    -------
    float
        The proportion of bases in the DNA sequence which are either ``'C'``
        or ``'G'``.

    Examples
    --------
    >>> gc_content("CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG")
    0.3125
    """
    return sum(1 for c in s if c in ['C' or 'G']) / len(s)


def max_gc_content(fasta_records: Bio.SeqIO.FastaIO.FastaIterator | typing.Iterable[Bio.SeqRecord.SeqRecord], /) -> tuple[str, float]:
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
    ('Rosalind_6127', 26.282722513089006)
    """
    sorted_gc_contents = sorted(
        {
            record.id: gc_content(record.seq) * 100
            for record in fasta_records
        }.items(),
        key=lambda t: t[1]
    )

    return sorted_gc_contents[-1]


def point_mutations(s: str | Bio.Seq.Seq, t: str | Bio.Seq.Seq, /) -> int:
    """:py:class:`int` : Returns a (non-negative) integer count of the mutations (or differences) in two DNA sequences of equal length.

    Solution to the Counting Point Mutations problem (HAMM):

    https://rosalind.info/problems/hamm/

    This is an application of the **Hamming distance** metric from computer
    science, as applied to DNA sequences.

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
    return sum(1 for i in range(len(s)) if s[i] != t[i])


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
    encoding = ''
    i = 0

    while i < len(s) - 2:
        codon = s[i: i + 3]
        encoding += RNA_CODON_TABLE[codon]
        i += 3

    return encoding


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
    n = len(s)
    m = len(t)

    return tuple([
         i + 1 for i in range(n - m + 1)
         if s[i: i + m] == t
    ])


def protein_mass(s: str | Bio.Seq.Seq, /) -> float:
    """:py:class:`float` : Calculates the mass of a protein sequence using the monoisotopic mass table.

    Solution to the Calculating Protein Mass problem (PRTM):

    https://rosalind.info/problems/prtm/

    Parameters
    ----------
    s : str, Bio.Seq.Seq
        The protein sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    float
        The protein string mass.

    Examples
    --------
    >>> protein_mass("SKADYEK")
    821.39192
    """
    return sum(MONOISOTOPIC_MASS_TABLE[c] for c in s)


def longest_common_substring(strs: typing.Iterable[str | Bio.Seq.Seq], /) -> str:
    """:py:class:`str` : Returns a longest common substring among an iterable of strings.

    Solution to the Finding a Shared Motif problem (LCSM):

    https://rosalind.info/problems/lcsm/

    Parameters
    ----------
    str : typing.Iterable
        An iterable of strings, which may be given as plain strings or
        Biopython genetic sequences (:py:class:`Bio.Seq.Seq`).

    Returns
    -------
    str
        A longest common substring (not necessarily unique) - the first
        one encountered is returned.

    Examples
    --------
    >>> seqs = ["GATTACA", "TAGACCA", "ATACA"]
    >>> longest_common_substring(seqs)
    'TA'
    """
    # Some key initial steps, starting with a sort of the incoming iterable
    # by length, in order to get a maximum length string.
    sorted_strs = sorted(strs, key=lambda s: len(s))
    max_len_str = sorted_strs[-1]
    max_str_len = len(max_len_str)
    str_len = len(max_len_str)

    # The outermost ``while`` loop on substring length, which descends from the
    # length of a largest substring (not necessarily unique) to minimal
    # substrings (of length 1).
    while str_len > 0:
        i = 0
        # The innermost ``while`` loop on substrings of length ``substr_len``
        while i < max_str_len - str_len + 1:
            cur_substr = max_len_str[i: i + str_len]
            # If the current substring isn't common then skip to the next
            # substring
            if any(cur_substr not in str for str in sorted_strs):
                i += 1
                continue

            # The current substring must be the common largest, as it hasn't
            # failed the checks, so return it.
            return cur_substr

        # Otherwise decrement the substring length, and start the next
        # iteration of the outer loop
        str_len -= 1

    # At this point the LCS must be an empty string, as there was no non-empty
    # substring returned in the inner loop, so just return it.
    return ''


def sequence_distance_matrix(seqs: typing.Iterable[str | Bio.Seq.Seq], /) -> list[list[float]]:
    """:py:class:`list` : Returns a (symmetric) matrix of relative Hamming distances for all (ordered) pairs of sequences from an iterable of equal-length genetic sequences.
    
    Solution to the Creating a Distance Matrix problem (PDST):

    https://rosalind.info/problems/pdst/

    As the distance between a sequence and itself is zero, and the distance
    between two sequences is independent of order, the matrix will be symmetric
    with zeros on the diagonal.

    Parameters
    ----------
    seqs : typing.Iterable
        An iterable of equal-length, plain Python strings or genetic sequences
        given as Biopython sequences (:py:class:`Bio.Seq.Seq`).

    Returns
    -------
    list
        An ``m`` by ``m`` symmetric matrix of relative Hamming distances of all
        pairs of sequences from the given iterable of sequences, where ``m`` is
        the common (equal) length of the sequences.

    Examples
    --------
    >>> seqs = ["TTTCCATTTA", "GATTCATTTC", "TTTCCATTTT", "GTTCCATTTA"]
    >>> sequence_distance_matrix(seqs)
    [[0.0, 0.4, 0.1, 0.1],
     [0.4, 0.0, 0.4, 0.3],
     [0.1, 0.4, 0.0, 0.2],
     [0.1, 0.3, 0.2, 0.0]]
    """
    n = len(seqs)

    # Set up ``n x n`` matrix of zeros - note that this construction below
    # using a list comprehension is designed to ensure that all of the zero
    # arrays are different objects in memory.
    mat = []
    [mat.append(list([0.] * n)) for i in range(n)]

    # A variable to store the common sequence length from the length of the
    # 1st sequence
    m = len(seqs[0])

    # The outer loop on rows
    i = 0
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
                mat[i][j] = point_mutations(seqs[i], seqs[j]) / m
                mat[j][i] = mat[i][j]
            j += 1
        i += 1

    return mat


def profile_matrix(seqs: typing.Iterable[str | Bio.Seq.Seq]) -> tuple:
    """:py:class:`tuple` : Returns the consensus string and profile matrix for a collection of DNA sequences/strings.

    Solution to the Consensus and Profile problem (CONS):

    https://rosalind.info/problems/cons/

    Parameters
    ----------
    seqs : typing.Iterable
        A collection of DNA sequences (or strings).

    Returns
    -------
    tuple
        The consensus string and profile matrix for the DNA sequences.

    Examples
    --------
    >>> seqs = ['ATCCAGCT', 'GGGCAACT', 'ATGGATCT', 'AAGCAACC', 'TTGGAACT', 'ATGCCATT', 'ATGGCACT']
    >>> cs, pm = profile_matrix(seqs)
    >>> print(cs)
    ATGCAACT
    >>> for i, b in enumerate('ACGT'):
    ...     print(f'{c}: ' + ' '.join(map(str, pm[i])))
    A: 5 1 0 0 5 5 0 0
    C: 0 0 1 4 2 0 6 1
    G: 1 1 6 3 0 1 0 0
    T: 1 5 0 0 0 1 1 6
    """
    # Take any sequence, but we can use the initial one to set
    # the (column) width of the matrix.
    n = len(seqs[0])

    # The initial ``4 x n`` profile matrix of zeros
    profile_matrix = [list([0] * n) for i in range(4)]

    # A list of blanks to store the characters of the consensus string
    consensus_str = list([''] * n)

    # The bases string
    base_str = 'ACGT'

    # The main loop
    for j in range(n):
        for i, b in enumerate(base_str):
            profile_matrix[i][j] = sum(1 for s in seqs if s[j] == b)
            if profile_matrix[i][j] == max(profile_matrix[x][j] for x in range(4)):
                consensus_str[j] = b

    return ''.join(consensus_str), profile_matrix


def signed_permutations(n: int, /) -> typing.Generator[tuple[int], None, None]:
    """:py:class:`typing.Generator` : Returns a generator of signed permutations of length ``n`` of the integer range ``1..n``.

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
    >>> list(signed_permutations(3))
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
    yield from chain.from_iterable(
        permutations(p)
        for p in product(*([k, -k] for k in range(1, n + 1)))
    )


def lexicographic_kmers(s: typing.Sequence[str], k: int) -> typing.Generator[str, None, None]:
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
    yield from map(lambda s: ''.join(s), product(s, repeat=k))


def kmer_composition(s: str | Bio.Seq.Seq, A: str, k: int) -> typing.Generator[int, None, None]:
    """:py:class:`typing.Generator` : Generates the ``k``-mer composition of a string with respect to an (ordered) alphabet.

    Solution to the k-Mer Composition problem (KMER):

    https://rosalind.info/problems/kmer/

    .. note::

       The function has been implemented as a generator, because ``k``-mers
       can grow very rapidly depending on the size of the alphabet.

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
        ``i``-value is the frequency of the ``i``-th kmer of ``A`` in ``s``, and
        ``i`` is the index of the kmer in the natural ordering of all ``k``-mers
        of ``A``.
    """
    n = len(s)

    def kmer_count(kmer: str, s: str) -> int:
        return sum(1 for i in range(n - k + 1) if s[i: i + k] == kmer)

    yield from (kmer_count(kmer, s) for kmer in lexicographic_kmers(A, k=k))


def variable_length_lexicographic_ordering(s: typing.Sequence[str], k: int) -> typing.Generator[str, None, None]:
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
    # Map characters in the sequence to their (1-indexed array) indices.
    seq_charmap = {char: idx + 1 for idx, char in enumerate(s)}

    # A lexicographic scoring function for substrings (given as strings or list
    # of one-letter strings) that uses the sequence character map on the
    # substring to build a tuple of index positions of its characters. The
    # resulting tuples can be compared and ordered lexicographically.
    def lex_score(t: str | list) -> str:
        return tuple(seq_charmap[char] for char in t)

    # Map substrings of length at most ``k`` to their lex scores.
    subseqs = {
        ''.join(p): lex_score(p)
        for p in chain.from_iterable(product(s, repeat=j) for j in range(1, k + 1))
    }

    # Now sort them by their lex scores, and generate them.
    yield from sorted(subseqs, key=lambda s: subseqs[s])


def linguistic_sequence_complexity(s: str, A: set, /) -> float:
    """:py:class:`float` : Returns the linguistic complexity (LC) of a stringt.

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

    A : set
        The underlying alphabet of the string.

    Returns
    -------
    float
        The linguistic complexity of the string.

    Examples
    --------
    >>> linguistic_sequence_complexity("ATTTGGATT")
    0.875
    """
    n = len(s)
    m = len(A)

    # Calculate the number of observed substrings as the sum of the number of
    # ``k``-length substrings found in the string, for ``k`` in ``1..n``.
    num_observed_substrings = len(set(s)) + sum(
        len(set([s[j: j + k] for j in range(n - k + 1)]))
        for k in range(2, n + 1)
    )

    # Calculate the number of possible substrings as the sum of the maximum
    # number of ``k``-length substrings formed from the alphabet that can fit
    # into the given string, overlapping as necessary. For each ``k`` the
    # maximum number of possible ``k``-length substrings that can fit into
    # the string is equal to ``min(m^k, n - k + 1)``.
    num_possible_substrings = sum(min(m ** k, n - k + 1) for k in range(1, n + 1))

    return num_observed_substrings / num_possible_substrings
