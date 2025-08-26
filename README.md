Solutions to the ROSALIND problem set
=====================================

Selected Python solutions to the [ROSALIND](https://rosalind.info/) bioinformatics problem set, including some data structures and utility functions.

Notes
-----

* The solutions aren't complete yet, but will added over time. All solutions are original, and the focus is not just on correctness, but conciseness, readability and speed/performance, although it is not always possible to achieve all three at once, and there is no absolute guarantee of performance in case of arbitrarily large inputs.

* [solutions.py](https://github.com/sr-murthy/rosalind/blob/main/src/solutions.py) is the main solution set, while [utils.py](https://github.com/sr-murthy/rosalind/blob/main/src/utils.py) contains generic utilities which are used in the solutions, as required.

* A basic set of [tests]((https://github.com/sr-murthy/rosalind/blob/main/src/tests/test_solutions.py)) has been added, and they mostly use the example problems from the ROSALIND problems as test cases.

* Solutions always produce raw values, and don't depend on formatting, e.g. [GC](https://rosalind.info/problems/gc/) (Computing GC Content), where the marking on ROSALIND depends on formatting the answer in a particular way.

* The function docstrings are written using the [Numpy docstring style](https://numpydoc.readthedocs.io/en/latest/format.html).

* For more background on [linguistic complexity (LC)](https://rosalind.info/problems/ling/) [see](https://pmc.ncbi.nlm.nih.gov/articles/PMC441604/pdf/gkh466.pdf), page W630.

* The counting of kmers in the [KMER](https://rosalind.info/problems/kmer/) (k-Mer Composition) problem must take overlapping substrings into account: so a [custom function](https://github.com/sr-murthy/rosalind/blob/main/src/utils.py#L155) has been used for this purpose, as [`str.count`](https://docs.python.org/3/library/stdtypes.html#str.count) only counts non-overlapping occurences.

* The solutions to several problems, including [SSEQ](https://rosalind.info/problems/sseq/) (Finding a Spliced Motif), involve finding and returning arrays of indices of a matching subsequence or substring, in terms of 1-indexed arrays, as required by the problems. They convert the 0-indexed array indices returned by some generic utility functions that they call on.

* The solution to [EDIT](https://rosalind.info/problems/edit/) (Edit Distance) is (cached) recursive, which is slower than equivalent iterative implementations, but is definitely more readable and easier to understand. It also allows for insertion, deletion, and substitution costs to be customised.
