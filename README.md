# Finding all variations of input DNA sequences in other sequences

### Inspiration:
In molecular biology you frequently need to cut out a specific piece of DNA using molecular "scissors". These "scissors" recognise specific sequences in long stretches of DNA. Unfortunalyt the recognition is ambigious, unlike seraching for an english word in a large document. For example the scissor may recognise the sequence ATTGCC which is very specific or it may recognise ANTGCC which expands to any of the following sequences: ATTGCC, AGTGCC, ACTGCC or AATGCC. Recognition sequences get much more complex than this simple example wiht th recogtion sequence expanding out to 20-30 different sequences all recognised by the same "pair of scissors"<p>
The aim was to develop a tool which was capable of finding all expanded input sequences at all locations in a long DNA sequence.

Techniques used:
- Tree based data strcutre with a depth first recursive search

Languages/Packages:
- Python 3.9 (PyCharm 2021 Community Editon IDE from JetBrains)
