# Finding all variations of input DNA sequences in other sequences

### Inspiration:
In molecular biology, it is frequently necessary to “cut out” a specific piece of DNA using molecular "scissors", called Restriction Enzymes. These "scissors" recognise specific sequences in long stretches of DNA. Unfortunately, the recognition process is ambiguous so unlike searching for an English word in a large document, there may be many spellings which are correct.
Using variations in spellings as an example, searching for the English and American spelling of colour means that both colour and color would be considered a match. Whilst this is simple in these examples, the biological world is many, many times more complex. Not only can we have variations in length we can have the same length but with different spellings.
In biology we use the [IUPAC ambiguity codes](https://www.bioinformatics.org/sms/iupac.html) to represent the various “spellings”. For example, the unique, non-ambiguous, sequence ATTGCC is very specific and has no variations or ambiguity. However, the letter N in a sequence means it can be any of the 4 bases A, T, G or C so ANTGCC would expand to any of the following sequences: ATTGCC, AGTGCC, ACTGCC or AATGCC.
Again not to complex, however, the following are a few examples of real world Restriction Enzyme recognition sequences and the number of combinations they expand into.
- NNCASTGNN -> 512 (S can be either G or C and N can be any of A, T, G or C)
- CGANNNNNNTGC -> 4096 combinations
- CCANNNNNNNNNTGG -> 262,144 combinations
To manage this particular data I decided to use a tree data structure (as many sequences would “share” most of their sequence with only a few differences between them and then to use a recursive depth first search to find the matches. I chose this approach because linked lists are memory efficient, can be searched quickly and technically it was a new challenge.
This tool is capable of searching DNA sequences many thousands of characters long and listing all the matches it finds from a file of over 100 restriction enzyme sequences, which vary in length and combinatorial complexity.

### Techniques used/Learnings:
- Tree based data structure with a depth first recursive search

### Languages/Packages:
- Python 3.9 (PyCharm 2021 Community Editon)
