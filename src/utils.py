__all__ = [
    'find_subsequence',
    'find_substring',
    'hamming_difference',
    'hamming_distance',
    'levenshtein_distance',
    'longest_common_substring',
    'signed_permutations',
    'word_grams',
    'word_k_grams',
]


"""
Generic utilities.
"""


# -- IMPORTS --

# -- Standard libraries --
import functools
import typing

from itertools import batched, chain, compress, groupby, permutations, product

# -- 3rd party libraries --

# -- Internal libraries --


@functools.cache
def find_subsequence(s: str, t: str) -> tuple[int] | typing.Literal[()]:
    """:py:class:`tuple` : Returns a tuple of ``s`-indices of elements of a string ``t`` occurring in ``s`` as a subsequence, or an empty tuple if it is not the case.

    .. note::

       The indices are given in terms of 0-indexed arrays. All elements of
       ``t`` must occur in ``s`` and in order, if not an empty tuple is
       returned.

    Parameters
    ----------
    s : str
        The string in which to search for ``t`` as a subsequence.

    t : str
        The string which has to be found in ``s`` as a subsequence.

    Returns
    -------
    tuple
        A tuple of indices of elements of ``t`` occurring in ``s`` as a
        subsequence, or the empty tuple ``()`` if not found.

    Examples
    --------
    >>> find_subsequence("AAABAAC", "ABAC")
    (0, 3, 4, 6)
    >>> find_subsequence("AAABAAC", "AD")
    ()
    >>> find_subsequence("", "AB")
    ()
    >>> find_subsequence("ABC", "")
    ()
    """
    # Exclude the edge cases where one (or both) of the strings is empty
    if not (s and t):
        return ()

    hits = []
    start = 0

    for c in t:
        i = s.find(c, start)
        if i != -1:
            hits.append(i)
            start = i + 1
        else:
            break

    return tuple(hits) if len(hits) == len(t) else ()


@functools.cache
def find_substring(s: str, t: str) -> tuple[int] | typing.Literal[()]:
    """:py:class:`tuple` : Returns a tuple of starting ``s``-indices of substring matches of a string ``t`` in a string ``s``, or an empty tuple if not found.

    .. note::

       The indices are given in terms of 0-indexed arrays.

    Parameters
    ----------
    s : str
        The string in which to search for the substring.

    t : str
        The substring to search for in the main string.

    Returns
    -------
    tuple
        A tuple of ints of starting ``s``-indices of substring matches, or the
        empty tuple ``()`` if none are found.

    Examples
    --------
    >>> find_substring("AAABAAC", "AA")
    (0, 1, 4)
    >>> find_substring("AAABAAC", "AD")
    ()
    >>> find_substring("", "AB")
    ()
    >>> find_substring("ABC", "")
    ()
    """
    # Exclude the edge cases where one (or both) of the strings is empty
    if not (s and t):
        return ()

    hits = []
    start = 0
    i = 0

    while True:
        i = s.find(t, start)
        if i != -1:
            hits.append(i)
            start = i + 1
        else:
            break

    return tuple(hits) if hits else ()


def hamming_difference(s: str, t: str, /) -> typing.Generator[tuple[int, tuple[str, str]], None, None]:
    """:py:class:`tuple` : Generates pairs of indices and corresponding elements that determine the Hamming distance between two equal-length strings.

    Utility function for the solution to the Transitions and Transversions
    problem (TRAN):

    https://rosalind.info/problems/tran/

    The "Hamming difference" is here defined as a sequence of the indices
    and the corresponding element pairs of the two sequences where differences
    occur. So that the Hamming distance of the sequences is simply the
    length of their Hamming difference.

    Parameters
    ----------
    s : str
        The first string.

    t : str
        The second string, equal in length to the first.

    Yields
    -------
    tuple
        Pairs of indices and corresponding element pairs where differences
        occur in the two strings.

    Raises
    ------
    ValueError
        If the two strings are not equal in length.

    Examples
    --------
    >>> list(hamming_difference("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"))
    [(0, ('G', 'C')),
     (2, ('G', 'T')),
     (4, ('C', 'G')),
     (7, ('C', 'A')),
     (9, ('A', 'G')),
     (14, ('G', 'C')),
     (15, ('A', 'C'))]
    
    """
    if len(s) != len(t):
        raise ValueError("Two equal length strings are required for the Hamming difference")

    yield from (
        (i, (a, b)) for i, (a, b) in enumerate(zip(s, t))
        if a != b
    )


@functools.cache
def hamming_distance(s: str, t: str, /) -> int:
    """:py:class:`int` : Returns the Hamming distance between two equal-length strings.

    Utility function for the solution to the Counting Point Mutations problem
    (HAMM):

    https://rosalind.info/problems/hamm/

    Calls on :py:func:`utils.hamming_distance`.

    Parameters
    ----------
    s : str
        The first string.

    t : str
        The second string, equal in length to the first.

    Returns
    -------
    int
        The Hamming distance between the two strings.

    Examples
    --------
    >>> hamming_distance("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"))
    7
    >>> hamming_distance("ACGT", "ACGT")
    0
    """
    return sum(1 for d in hamming_difference(s, t))


@functools.cache
def levenshtein_distance(
    s: str,
    t: str,
    /,
    *,
    insertion_cost: int | float  = 1,
    deletion_cost: int | float = 1,
    substitution_cost: int | float = 1
) -> float:
    """:py:class:`int` : Returns the Levenshtein distance between two strings.

    This is a utility function for the Edit Distance problem (EDIT):

    https://rosalind.info/problems/edit/

    This is a cached recursive implementation, with customisable costs for
    the insertion, deletion, and substitution operations, where these
    operations are to be performed on the first string, ``s``, to make it
    equal to the second string ``t``.

    .. note::

       This implementation accepts ``int`` or ``float``-valued cost parameters.

    Parameters
    ----------
    s : str
        The first string.

    t : str
        The second string, equal in length to the first.

    insertion_cost : int, float, default=1
        An optional insertion cost, defaulting to ``1``.

    deletion_cost : int, float, default=1
        An optional deletion cost, defaulting to ``1``.

    substitution_cost : int, float, default=1
        An optional substition cost, defaulting to ``1``.

    Returns
    -------
    float
        The Levenshtein distance between the two strings.

    Examples
    --------
    # Minimal examples with (a) insertion cost of 1 (for ``s``)
    >>> levenshtein_distance("read", "bread"))
    1
    # ... (b) deletion cost of 1
    >>> levenshtein_distance("bread", "read")
    1
    # ... (c) substitution cost of 1
    >>> levenshtein_distance("bread", "tread")
    1
    """
    # If ``s`` is empty then ``|t|`` insertions (into ``s``) are required to
    # make it equal to ``t``.
    if len(s) == 0:
        return len(t)

    # If ``t`` is empty then ``|s|`` deletions (from ``s``) are required to
    # make it equal to ``s``.
    if len(t) == 0:
        return len(s)

    # If the heads of both strings are equal there are no costs for these,
    # so compute the distance for the tails
    if s[0] == t[0]:
        return levenshtein_distance(s[1:], t[1:])

    return min(
        insertion_cost + levenshtein_distance(s, t[1:]),        # insertion of ``t[0]`` at ``s[0]``
        deletion_cost + levenshtein_distance(s[1:], t),         # deletion of ``s[0]``
        substitution_cost + levenshtein_distance(s[1:], t[1:])  # substitution of ``t[0]`` for ```s[0]``
    )


@functools.cache
def longest_common_substring(strs: tuple[str], /) -> str:
    """:py:class:`str` : Returns a longest common substring among an iterable of strings.

    Utility function for the solution to the Finding a Shared Motif problem
    (LCSM):

    https://rosalind.info/problems/lcsm/

    Parameters
    ----------
    strs : tuple
        A tuple of strings. The tuple requirement is due to the fact that
        the function uses a cache, which requires arguments to be hashable.

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
    # The longest common substring in a set of strings cannot exceed the
    # length of a smallest string in the set. So any smallest length
    # substring can be used to search for substrings common to all the
    # rest.

    # Some key initial steps, starting with getting a smallest length
    # substring.
    min_len_str = min(strs, key=len)
    min_str_len = len(min_len_str)
    str_len = len(min_len_str)

    # The outermost ``while`` loop on substring length, which descends from the
    # length of a smallest length string (not necessarily unique) to minimal
    # length strings (of length 1).
    while str_len > 0:
        i = 0
        # The innermost ``while`` loop on substrings of length ``substr_len``
        while i < min_str_len - str_len + 1:
            cur_substr = min_len_str[i: i + str_len]
            # If the current substring isn't common then skip to the next
            # substring
            if any(cur_substr not in s for s in strs):
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


def signed_permutations(n: int, /) -> typing.Generator[tuple[int], None, None]:
    """:py:class:`typing.Generator` : Returns a generator of signed permutations of length ``n`` of the integer range ``1..n``.

    Utility function for the solution to the Enumerating Oriented Gene
    Orderings problem (SIGN):

    https://rosalind.info/problems/sign/

    .. note::

       The number of signed permutations of :math:`n` is given by:

       .. math::

          (2^n)n! 

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


def word_grams(w: str, k: int) -> typing.Generator[str, None, None]:
    """:py:class:`typing.Generator` : Returns a generator of all 1-grams, 2-grams, ..., k-grams of a word ``w`` in lexicographic order.
    
    Utility function for the Solution to the Ordering Strings of Varying Length
    Lexicographically problem (LEXV):

    https://rosalind.info/problems/lexv/

    Parameters
    ----------
    w : str
        The word/string from which to generate the 1-, 2-,..., k-grams.

    k : int
        The maximum length of the word grams.

    Yields
    ------
    str
        The 1-, 2-, ... , k-grams in (lexicographic) order.

    Examples
    --------
    >>> list(word_grams('DNA', k=3))
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
    # Map characters in ``w`` to their (0-indexed array) indices.
    w_charmap = {c: i for i, c in enumerate(w)}

    # A lexicographic scoring function for the word grams which uses the
    # sequence character map of ``s`` to build a tuple of indices of a given
    # word gram. The resulting tuples can be compared and ordered
    # lexicographically.
    def lex_score(wg: str | list) -> tuple[int]:
        return tuple(w_charmap[c] for c in wg)

    # Map all word grams of length at most ``k`` to their lex scores.
    subs = {
        ''.join(p): lex_score(p)
        for p in chain.from_iterable(product(w, repeat=j) for j in range(1, k + 1))
    }

    # Now sort them by their lex scores, and generate them.
    yield from sorted(subs, key=lambda wg: subs[wg])


def word_k_grams(w: str, k: int) -> typing.Generator[str, None, None]:
    """:py:class:`typing.Generator` : Returns a generator of all ``k``-grams of a given word ``w`` in lexicographic order.
    
    Utility function for the solution to the Enumerating k-mers
    Lexicographically problem (LEXF):

    https://rosalind.info/problems/lexf/

    Parameters
    ----------
    w : str
        The word/string from which to generate the ``k``-grams.

    k : int
        The (fixed) word gram length.

    Yields
    ------
    str
        ``k``-grams of the given word ``w`` generated in (lexicographic)
        order.

    Examples
    --------
    >>> list(word_k_grams('ACGT', k=2))
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
    yield from map(lambda p: ''.join(p), product(w, repeat=k))
