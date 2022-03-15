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
Sample Dataset:

AAAACCCGGT

Sample Output:

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



# Rabbits and Recurrence Relations 
[detailed introduction to this question](https://rosalind.info/problems/fib/)

A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences can be finite or infinite. Two examples are the finite sequence $(π,-\sqrt{2},0,π)$ and the infinite sequence of odd numbers $(1,3,5,7,9,…)$. We use the notation $a_n$ to represent the $n-th$ term of a sequence.

A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal to the number of rabbits that were alive two months prior. As a result, if Fn represents the number of rabbit pairs alive after the n-th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence relation $Fn=F_{n−1}+F{n−2}$ (with $F_1=F_2=1$ to initiate the sequence). Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.

When finding the $n-th$ term of a sequence defined by a recurrence relation, we can simply use the recurrence relation to generate terms for progressively larger values of n. This problem introduces us to the computational technique of dynamic programming, which successively builds up solutions by using the answers to smaller cases.

Given: Positive integers $n≤40$ and $k≤5$.

Return: The total number of rabbit pairs that will be present after $n$ months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of $k$ rabbit pairs (instead of only $1$ pair).

## Sample
Sample Dataset:

5 3

Sample Output:

19   


```python
def Rabbit_rec(n,k):
    if n==1 or n==2:
        return 1
    else:
        return Rabbit_rec(n-1,k)+k*Rabbit_rec(n-2,k)
def Rabbit_dp(n,k):
    r=[0]*(n+1)
    r[1]=r[2]=1
    for i in range(3,n+1):
        r[i]=r[i-1]+k*r[i-2]
    return r[n]
```


```python
Rabbit_dp(5,3)
```




    19



# Computing GC Content 
[detailed introduction to this question](https://rosalind.info/problems/gc/)

The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.

DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.

In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.

## Sample
Sample Dataset:

\>Rosalind_6404

CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG

\>Rosalind_5959

CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC

\>Rosalind_0808

CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT

Sample Output:

Rosalind_0808

60.919540


```python
def GC_count(file_path=""):
    f=open(file_path)
    DNA_dict={}
    DNA_name=" "
    while DNA_name!="":
        DNA_name=f.readline().strip()[1::]
        DNA_dict[DNA_name]=f.readline().strip()
    DNA_dict.pop("")
    GC_dict={}
    for DNA_name in DNA_dict.keys():
        GC_dict[DNA_name]=round((DNA_dict[DNA_name].count('G')+DNA_dict[DNA_name].count('C'))/len(DNA_dict[DNA_name])*100,6)
    max_DNA=""
    max_GC=0
    for DNA_name in GC_dict:
        if GC_dict[DNA_name]>max_GC:
            max_DNA=DNA_name
            max_GC=GC_dict[DNA_name]
    print(max_DNA)
    print(max_GC)
```


```python
GC_count("gccontent.fa")
```

    Rosalind_0808
    60.91954
    

# Evolution as a Sequence of Mistakes 
[detailed introduction to this question](https://rosalind.info/problems/hamm/)

Given two strings s and t of equal length, the Hamming distance between s and t, denoted $d_H(s,t)$, is the number of corresponding symbols that differ in s and t. See the figure.

![point mistake](https://rosalind.info/media/problems/hamm/Hamming_distance.png "The Hamming distance between these two strings is 7. Mismatched symbols are colored red.")

<center><font color=gray size=0.5>Figure: The Hamming distance between these two strings is 7. Mismatched symbols are colored red.</font></center>

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

Return: The Hamming distance $d_H(s,t)$.

## Sample
Sample Dataset:

GAGCCTACTAACGGGAT
CATCGTAATGACGGCCT

Sample Output:

7


```python
def Countpm(s,t):
    i=0
    pmsum=0
    while i<len(s):
       if s[i]!=t[i]:
           pmsum=pmsum+1
       i=i+1
    print(pmsum)
```


```python
Countpm("GAGCCTACTAACGGGAT","CATCGTAATGACGGCCT")
```

    7
    
