# A Rapid Introduction to Molecular Biology
[detailed introduction of this question](https://rosalind.info/problems/dna/)

A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains.

An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."

Given: A DNA string s of length at most 1000 nt.

Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.

Sample Dataset: 

AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC

Sample Output: 

20 12 17 21


```python
def countn(DNA_str):
    DNA_str=DNA_str.upper() # Sometime, our DNA string is not completely consists of capital letters
    i=countA=countT=countG=countC=0
    while i<len(DNA_str):
        if DNA_str[i]=='A':
            countA=countA+1
        elif DNA_str[i]=='T':
            countT=countT+1
        elif DNA_str[i]=='G':
            countG=countG+1
        else:
            countC=countC+1
        i=i+1
    print("Acount is %d \n"%(countA))
    print("Tcount is %d \n"%(countT))
    print("Gcount is %d \n"%(countG))
    print("Ccount is %d \n"%(countC))

def countn2(DNA_str):
    DNA_str=DNA_str.upper()
    DNA_list=list(DNA_str)
    print("Acount is %d \n"%(DNA_list.count('A')))
    print("Tcount is %d \n"%(DNA_list.count('T')))
    print("Gcount is %d \n"%(DNA_list.count('G')))
    print("Ccount is %d \n"%(DNA_list.count('C')))

```


```python
countn('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC')
```

    Acount is 20 
    
    Tcount is 21 
    
    Gcount is 17 
    
    Ccount is 12 
    
    


```python

```
