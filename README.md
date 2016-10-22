# deBruijnGraphAssembler

**Usage:**
```
>> python contigGenerator.py <directory of input file> <kmer size> <coverage filter> <weighted edge filter>
```
For example:
```
>> python contigGenerator.py tests\example.data.fasta 31 0 0
```
Note that on Unix machines you would use '/' instead of '\' for the directory. Also note how I used a coverage filter and a weighted edge filter of 0; because the example.data.fasta file has no reads there is no need to filter any kmers/edges.

**Remember that the kmer size MUST be less than or equal to the read size otherwise the output will be empty!**
