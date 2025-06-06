__all__ = [
    'basecount',
    'count_dna_motif',
    'sequence_distance_matrix',
    'gc_content',
    'lexicographic_kmers',
    'linguistic_sequence_complexity',
    'longest_common_sequence_motif',
    'max_gc_content',
    'MONOISOTOPIC_MASS_TABLE',
    'point_mutations',
    'protein_mass',
    'reverse_complement',
    'RNA_CODON_TABLE',
    'signed_permutations',
    'transcribe_dna_to_rna',
    'translate_rna_to_protein',
    'variable_length_lexicographic_ordering',
]


# -- IMPORTS --

# -- Standard libraries --
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


def basecount(dna_seq: str | Bio.Seq.Seq, /) -> tuple[int]:
    """:py:class:`str` : Returns a tuple of integer counts of the DNA bases (``'A'``, ``'C'``, ``'G'``, ``'T'`` in that order) from a DNA sequence.

    Parameters
    ----------
    dna_seq : str, Bio.Seq.Seq
        The input DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    tuple
        A tuple of integer counts of the DNA bases (``'A'``, ``'C'``,
        ``'G'``, ``'T'`` in that order) from the given sequence.
    """
    counter = Counter(dna_seq)
    
    return counter['A'], counter['C'], counter['G'], counter['T']


def transcribe_dna_to_rna(dna_seq: str | Bio.Seq.Seq, /) -> str:
    """:py:class:`str` : Returns an RNA transcription of a DNA string.

    Parameters
    ----------
    dna_seq : str
        The input DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    A string of the RNA sequence transcribed from the given DNA sequence.
    """
    return str(dna_seq).translate(str.maketrans('T', 'U'))


def reverse_complement(dna_seq: str | Bio.Seq.Seq, /) -> str:
    """:py:class:`str` : Returns the reverse complement of a DNA string, which is the reversed string with the bases complemented.

    Parameters
    ----------
    dna_seq : str
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
    """
    return str(dna_seq)[::-1].translate(str.maketrans("ATCG", "TAGC"))


def gc_content(dna_seq: str | Bio.Seq.Seq, /) -> float:
    """:py:class:`float` : The proportion of bases in the DNA sequence which are either ``'C'`` or ``'G'``.

    .. note::

       The returned value is a proportion (or ratio), not a percentage.

    Parameters
    ----------
    dna_seq : str, Bio.Seq.Seq
        The DNA sequence.

    Returns
    -------
    float
        The proportion of bases in the DNA sequence which are either ``'C'``
        or ``'G'``.

    """
    return sum(1 for c in dna_seq if c in ['C' or 'G']) / len(dna_seq)


def max_gc_content(fasta_records: Bio.SeqIO.FastaIO.FastaIterator | typing.Iterable[Bio.SeqRecord.SeqRecord], /) -> tuple[str, float]:
    """:py:class:`str` : Returns a tuple containing a record ID and GC content percentage for the FASTA record with the highest GC content in a FASTA file.

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
    """
    sorted_gc_contents = sorted(
        {
            record.id: gc_content(record.seq) * 100
            for record in fasta_records
        }.items(),
        key=lambda t: t[1]
    )

    return sorted_gc_contents[-1]


def point_mutations(dna_seq1: str | Bio.Seq.Seq, dna_seq2: str | Bio.Seq.Seq, /) -> int:
    """:py:class:`int` : Returns a (non-negative) integer count of the mutations (or differences) in two DNA sequences of equal length.

    This is an application of the **Hamming distance** metric from computer
    science, as applied to DNA sequences.

    Parameters
    ----------
    dna_seq1 : str, Bio.Seq.Seq
        The first DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    dna_seq2 : str, Bio.Seq.Seq
        The second DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    int
        A (non-negative) integer count of the mutations in the two given DNA
        sequences.
    """
    return sum(1 for i in range(len(dna_seq1)) if dna_seq1[i] != dna_seq2[i])


def translate_rna_to_protein(rna_seq: str | Bio.Seq.Seq, /) -> str:
    """:py:class:`str` : Returns a protein sequence encoded by an RNA seq.

    .. note::

       An assumption is made that the length of the RNA sequence is a multiple
       of 3, as codons in the RNA codon table are all of length 3.

    Parameters
    ----------
    rna_seq : str, Bio.Seq.Seq
        The RNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    str
        The protein sequence encoded by the RNA sequence.

    """
    encoding = ''
    i = 0

    while i < len(rna_seq) - 2:
        codon = rna_seq[i: i + 3]
        encoding += RNA_CODON_TABLE[codon]
        i += 3

    return encoding


def count_dna_motif(dna_seq: str | Bio.Seq.Seq, dna_subseq: str | Bio.Seq.Seq, /) -> tuple[int]:
    """:py:class:`tuple` : A tuple of all starting indices of occurrences of a DNA subseqence in the given DNA sequence.

    .. note::

       The answer is given in terms of indices of 1-indexed arrays - to get
       zero indices just subtract 1.

    Parameters
    ----------
    dna_seq : str, Bio.Seq.Seq
        The DNA sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    dna_subseq : str, Bio.Seq.Seq
        The DNA subsequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    tuple
        Tuple of starting positions of all occurrences of the DNA subseqence
        in the given DNA sequence.

    """
    n = len(dna_seq)
    m = len(dna_subseq)

    return tuple([
         i + 1 for i in range(n - m + 1)
         if dna_seq[i: i + m] == dna_subseq
    ])


def protein_mass(protein_seq: str | Bio.Seq.Seq, /) -> float:
    """:py:class:`float` : Calculates the mass of a protein sequence using the monoisotopic mass table.

    Parameters
    ----------
    protein_seq : str, Bio.Seq.Seq
        The protein sequence, given as a string or :py:class:`Bio.Seq.Seq`.

    Returns
    -------
    float
        The protein string mass.
    """
    return sum(MONOISOTOPIC_MASS_TABLE[c] for c in protein_seq)


def longest_common_sequence_motif(seqs: typing.Iterable[str | Bio.Seq.Seq], /) -> str:
    """:py:class:`str` : Returns a longest common sequence motif from an iterable of sequences.

    Parameters
    ----------
    dna_seqs : typing.Iterable
        An iterable of sequences, which may be given as plain strings or
        Biopython genetic sequences (:py:class:`Bio.Seq.Seq`).

    Returns
    -------
    str
        A longest common sequence motif (substring) from an iterable of
        sequences.
    """
    # Some key initial steps, starting with a sort of the incoming sequence
    # iterable, in order to get a maximum length sequence.
    sorted_seqs = sorted(seqs, key=lambda seq: len(seq))
    max_len_seq = sorted_seqs[-1]
    max_seq_len = len(max_len_seq)
    seq_len = len(max_len_seq)
    longest_common_motif = ''

    # The outermost ``while`` loop on substring length, which descends from the
    # length of a largest substring (not necessarily unique) to minimal
    # substrings (of length 1).
    while seq_len > 0:
        i = 0
        # The innermost ``while`` loop on substrings of length ``substr_len``
        while i < max_seq_len - seq_len + 1:
            cur_seq = max_len_seq[i: i + seq_len]
            # If the current substring isn't common or isn't longer than the
            # current longest common motif then skip to the next substring
            # sliding window
            if any(cur_seq not in seq for seq in sorted_seqs) or len(cur_seq) <= len(longest_common_motif):
                i += 1
                continue

            # The current substring must be the common largest, as it hasn't
            # failed the checks, so return it.
            return cur_seq

        # Otherwise decrement the substring length, and start the next
        # iteration of the outer loop
        seq_len -= 1

    # At this point the LCM must be an empty string, as it wasn't returned
    # in the loop, so just return it.
    return longest_common_motif


def sequence_distance_matrix(seqs: typing.Iterable[str | Bio.Seq.Seq], /) -> list[list[float]]:
    """:py:class:`list` : Returns a (symmetric) matrix of relative Hamming distances for all (ordered) pairs of sequences from an iterable of equal-length genetic sequences.
    
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
            # the Hamming distance of the current sequence pair, and set
            # the transpose value to the same, avoiding another calculation
            else:
                mat[i][j] = point_mutations(seqs[i], seqs[j]) / m
                mat[j][i] = mat[i][j]
            j += 1
        i += 1

    return mat


def signed_permutations(n: int, /) -> typing.Generator[list[int], None, None]:
    """:py:class:`typing.Generator` : Returns a generator of signed permutations of length ``n`` of the integer range ``1..n``.

    Parameters
    ----------
    n : int
        The (common) length of the signed permutations of ``1..n``.

    Yields
    ------
    list
        A signed permutation as a list of ints.
    """
    yield from chain.from_iterable(
        permutations(p)
        for p in product(*([k, -k] for k in range(1, n + 1)))
    )


def lexicographic_kmers(seq: typing.Sequence[str], /) -> typing.Generator[str, None, None]:
    """:py:class:`typing.Generator` : Returns a generator of ``k``-length substrings formed from subsequences of a given sequence of characters.
    
    Parameters
    ----------
    seq : typing.Sequence
        A sequence of characters, which can be any single letter Unicode
        strings. The given ordering is used in building lexicographically
        ordered substrings.

    Yields
    ------
    str
        ``k``-length substrings formed from the sequence generted in
        lexicographic order (with respect to the original sequence).
    """
    yield from map(
        lambda s: ''.join(s),
        product(seq, repeat=k)
    )


def variable_length_lexicographic_ordering(seq: typing.Sequence[str], k: int, /) -> typing.Generator[str, None, None]:
    """:py:class:`typing.Generator` : Returns a generator of lexicographically ordered, variable-length substrings of size at most ``k``, formed from subsequences of a given sequence of characters.
    
    Parameters
    ----------
    seq : typing.Sequence
        A sequence of characters, which can be any single letter Unicode
        strings. The given ordering is used in building lexicographically
        ordered substrings.

    Yields
    ------
    str
        Variable-length substrings of length at most ``k``, formed from
        the sequence, generated in lexicographic order (with respect to
        the original sequence).
    """
    # A map of characters in the sequence to their (1-indexed) indices.
    seq_charmap = dict([(char, idx + 1) for idx, char in enumerate(seq)])

    # A lexicographic scoring function for substrings (given as strings or list
    # of one-letter strings) that uses the sequence character map on the
    # substring to build a tuple of index positions of its characters. The
    # resulting tuples can be compared and ordered lexicographically.
    def lex_score(subseq: str | list) -> str:
        return tuple(seq_charmap[char] for char in subseq)

    # A dict of all possible substrings of length at most ``k`` mapped to their
    # lex scores.
    subseqs = {
        ''.join(p): lex_score(p)
        for p in chain.from_iterable(product(seq, repeat=j) for j in range(1, k + 1))
    }

    # Now sort them by the lexicographic scoring function, and generate them.
    yield from sorted(
        subseqs,
        key=lambda subseq: lex_score(subseq)
    )


def linguistic_sequence_complexity(string: str, alphabet: set, /) -> float:
    """:py:class:`float` : Returns the linguistic complexity (LC) of a string formed over an alphabet.

    The linguistic complexity (LC) measure is defined in the ROSALIND problem page:

        https://rosalind.info/problems/ling/

    but more clearly defined in the following paper (page W630):

        Orlov YL, Potapov VN. Complexity: an internet resource for analysis
        of DNA sequence complexity. Nucleic Acids Res. 2004 Jul 1;32(Web Server
        issue):W628-33. doi: 10.1093/nar/gkh466. PMID: 15215465; PMCID:
        PMC441604.

    This implementation uses the LC formula from this paper, which is called
    CL.

    Parameters
    ----------
    string : str
        The input string.

    alphabet : set
        The underlying alphabet of the string.

    Returns
    -------
    float
        The linguistic complexity of the string.
    """
    n = len(string)
    m = len(alphabet)

    # Calculate the number of observed substrings as the sum of the number of
    # ``k``-length substrings found in the string, for ``k`` in ``1..n``.
    num_observed_substrings = len(set(string)) + sum(
        len(set([string[j: j + k] for j in range(n - k + 1)]))
        for k in range(2, n + 1)
    )

    # Calculate the number of possible substrings as the sum of the maximum
    # number of ``k``-length substrings formed from the alphabet that can fit
    # into the given string, overlapping as necessary. For each ``k`` the
    # maximum number of possible ``k``-length substrings that can fit into
    # the string is equal to ``min(m^k, n - k + 1)``.
    num_possible_substrings = sum(min(m ** k, n - k + 1) for k in range(1, n + 1))

    return num_observed_substrings / num_possible_substrings
