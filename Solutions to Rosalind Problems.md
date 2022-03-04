# Counting DNA Nucleotides
[detailed introduction to this question](https://rosalind.info/problems/dna/)

A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains.

An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."

Given: A DNA string s of length at most 1000 nt.

Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.

## Sample
Sample Dataset: 

AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC

Sample Output: 

20 12 17 21


```python
def Countn(DNA_str):
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

def Countn2(DNA_str):
    DNA_str=DNA_str.upper()
    DNA_list=list(DNA_str)
    print("Acount is %d \n"%(DNA_list.count('A')))
    print("Tcount is %d \n"%(DNA_list.count('T')))
    print("Gcount is %d \n"%(DNA_list.count('G')))
    print("Ccount is %d \n"%(DNA_list.count('C')))

```


```python
Countn('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC')
```

    Acount is 20 
    
    Tcount is 21 
    
    Gcount is 17 
    
    Ccount is 12 
    
    

# Complementing a Strand of DNA 
[detailed introduction to this question](https://rosalind.info/problems/revc/)

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.

The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s, then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

Given: A DNA string s of length at most 1000 bp.

Return: The reverse complement sc of s.

## Sample
Sample Dataset

AAAACCCGGT

Sample Output

ACCGGGTTTT


```python
def Complement(DNA_str):
    DNA_str=DNA_str.upper()
    DNA_str=DNA_str.replace('A','t').replace('T','a').replace('G','c').replace('C','g')
    Com_DNA_str=DNA_str[::-1]
    Com_DNA_str=Com_DNA_str.upper()
    return Com_DNA_str
```


```python
Complement('AAAACCCGGT')
```




    'ACCGGGTTTT'


