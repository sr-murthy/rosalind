Solutions to the ROSALIND problem set
=====================================

Selected Python solutions to the [ROSALIND](https://rosalind.info/) bioinformatics problem set, including some data structures and generic utilities.

Notes
-----

* The solutions are a work in progress, and will be added as time permits. All solutions are original, and the focus is not just on correctness, but also conciseness, speed/performance (or scalability), and readability. Although it is not always possible to achieve all three at once, and there is no absolute guarantee of performance for arbitrarily large inputs.

* [solutions.py](https://github.com/sr-murthy/rosalind/blob/main/src/solutions.py) is the main solution set, while [utils.py](https://github.com/sr-murthy/rosalind/blob/main/src/utils.py) contains generic utilities which are used in the solutions, as required.

* A basic set of [tests]((https://github.com/sr-murthy/rosalind/blob/main/src/tests/test_solutions.py)) has been added, and they are mostly based on solutions to the example problems in ROSALIND.

* Solutions always produce raw values, and don't depend on formatting, e.g. [GC](https://rosalind.info/problems/gc/) (Computing GC Content), where the marking in ROSALIND depends on formatting the answer in a particular way.

* In problems with numerical solutions [`decimal.Decimal`](https://docs.python.org/3/library/decimal.html#decimal.Decimal) objects are returned instead of `float` values, where possible, to ensure that results are as exact as possible.

* Several functions (in `utils.py` and `solutions.py`) use caching with [`functools.cache`](https://docs.python.org/3/library/functools.html#functools.cache), which requires that the arguments and parameters passed to the functions are hashable: in particular, all array-valued arguments and parameters to these functions must be tuples, because Python tuples are immutable and hence hashable.

* The function docstrings are written using the [Numpy docstring style](https://numpydoc.readthedocs.io/en/latest/format.html).

* For more background on [linguistic complexity (LC)](https://rosalind.info/problems/ling/) refer [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC441604/pdf/gkh466.pdf) (page W630).

* In Python the counting of k-mers in the [KMER](https://rosalind.info/problems/kmer/) (k-Mer Composition) problem must take overlapping substrings into account, because the Python standard library [`str.count`](https://docs.python.org/3/library/stdtypes.html#str.count) function only counts non-overlapping occurences: so a [custom function](https://github.com/sr-murthy/rosalind/blob/main/src/utils.py#L155) has been used for this purpose.

* The solutions to several problems, including [SSEQ](https://rosalind.info/problems/sseq/) (Finding a Spliced Motif), involve finding and returning arrays of indices of a matching subsequence or substring, in terms of 1-indexed arrays, as required by the problems. They convert the 0-indexed array indices returned by some generic utility functions that they call on.

* The solution to [EDIT](https://rosalind.info/problems/edit/) (Edit Distance) is (cached) recursive, which is slower than equivalent iterative implementations, but more readable and easier to understand. Also, and separately, the solution allows for insertion, deletion, and substitution costs to be customised, with default values of 1, 1, 1 respectively.
