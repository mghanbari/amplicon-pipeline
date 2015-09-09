1. join paired ends fastqjoin (hopefully with API)
2. Collapse identical sequences (100% otu clustering)
  a. VSEARCH or by using python dictionaries (DAWG)  
    {'atgtcagcat': {'sample1':43},
                   {'sample2': 3}}
3. remove sequences found in less than x% of sample
4. "Assign taxonomy"

    for seq in seqs:
      if seq in dictionary:
        add it to the sample
      else:
        dictionary[seq] = new entry

```python

```
