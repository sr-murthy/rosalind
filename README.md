Solutions to the ROSALIND problem set
=====================================

Selected Python solutions to the [ROSALIND](https://rosalind.info/) bioinformatics problem set, including some data structures and utility functions.

Notes
-----

* The marking of solutions to [GC](https://rosalind.info/problems/gc/) (Computing GC Content) seems to depend, unfortunately, on formatting the answer in a particular way. Even with the correct numerical answer, the solution can be marked wrong, e.g. if a percentage symbol is missing.

* Linguistic complexity (LC) is not as clearly defined on the [problem page](https://rosalind.info/problems/ling/). A clearer definition is [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC441604/pdf/gkh466.pdf), page W630.

* The solution to [PDST](https://rosalind.info/problems/pdst/) (Creating a Distance Matrix) uses the solution to the [HAMM](https://rosalind.info/problems/gc/) (Creating a Distance Matrix problem).

* The solution to [KMER](https://rosalind.info/problems/kmer/) (k-Mer Composition) uses the solution to [LEXF](https://rosalind.info/problems/lexf/) (Enumerating k-mers Lexicographically).

* The counting of kmers in the [KMER](https://rosalind.info/problems/kmer/) (k-Mer Composition) problem needs to take overlapping substrings into account, as [`str.count`](https://docs.python.org/3/library/stdtypes.html#str.count) only counts non-overlapping occurences.
