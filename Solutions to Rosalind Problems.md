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

A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences can be finite or infinite. Two examples are the finite sequence $(π,-\sqrt{2},0,π)$ and the infinite sequence of odd numbers (1,3,5,7,9,…). We use the notation a<sub>n</sub> to represent the n-th term of a sequence.

A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal to the number of rabbits that were alive two months prior. As a result, if F<sub>n</sub> represents the number of rabbit pairs alive after the n-th month, then we obtain the Fibonacci sequence having terms F<sub>n</sub> that are defined by the recurrence relation F<sub>n</sub>=F<sub>n−1</sub>+F<sub>n−2</sub> (with F<sub>1</sub>=F<sub>2</sub>=1 to initiate the sequence). Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.

When finding the n-th term of a sequence defined by a recurrence relation, we can simply use the recurrence relation to generate terms for progressively larger values of n. This problem introduces us to the computational technique of dynamic programming, which successively builds up solutions by using the answers to smaller cases.

Given: Positive integers n≤40 and k≤5.

Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).

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
    

# Counting Point Mutations
[detailed introduction to this question](https://rosalind.info/problems/hamm/)

Given two strings s and t of equal length, the Hamming distance between s and t, denoted d<sub>H</sub>(s,t), is the number of corresponding symbols that differ in s and t. See the figure.

![point mistake](https://rosalind.info/media/problems/hamm/Hamming_distance.png "The Hamming distance between these two strings is 7. Mismatched symbols are colored red.")

<center><font color=gray size=0.5>Figure: The Hamming distance between these two strings is 7. Mismatched symbols are colored red.</font></center>

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

Return: The Hamming distance d<sub>H</sub>(s,t).

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
    

# Mendel's First Law
[detailed introduction to this question](https://rosalind.info/problems/iprb/)

Probability is the mathematical study of randomly occurring phenomena. We will model such a phenomenon with a random variable, which is simply a variable that can take a number of different distinct outcomes depending on the result of an underlying random process.

For example, say that we have a bag containing 3 red balls and 2 blue balls. If we let X represent the random variable corresponding to the color of a drawn ball, then the probability of each of the two outcomes is given by Pr(X=red)=3/5 and Pr(X=blue)=2/5.

Random variables can be combined to yield new random variables. Returning to the ball example, let Y model the color of a second ball drawn from the bag (without replacing the first ball). The probability of Y being red depends on whether the first ball was red or blue. To represent all outcomes of X and Y, we therefore use a probability tree diagram. This branching diagram represents all possible individual probabilities for X and Y, with outcomes at the endpoints ("leaves") of the tree. The probability of any outcome is given by the product of probabilities along the path from the beginning of the tree; see Figure for an illustrative example.

![Mendelian Inheritance](https://rosalind.info/media/problems/iprb/balls_tree.thumb.png "The probability of any outcome (leaf) in a probability tree diagram is given by the product of probabilities from the start of the tree to the outcome. For example, the probability that X is blue and Y is blue is equal to (2/5)(1/4), or 1/10.")

<center><font color=gray size=0.5>Figure: The probability of any outcome (leaf) in a probability tree diagram is given by the product of probabilities from the start of the tree to the outcome. For example, the probability that X is blue and Y is blue is equal to (2/5)(1/4), or 1/10.</font></center>

An event is simply a collection of outcomes. Because outcomes are distinct, the probability of an event can be written as the sum of the probabilities of its constituent outcomes. For our colored ball example, let A be the event "Y is blue." Pr(A) is equal to the sum of the probabilities of two different outcomes: Pr(X=blue and Y=blue)+Pr(X=red and Y=blue), or 3/10+1/10=2/5 (see Figure above).

Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.

## Sample
Sample Dataset:

2 2 2

Sample Output:

0.78333


```python
def Mendelian(k,m,n):
    indsun=k+m+n
    return 1-m/indsun*(m-1)/(indsun-1)*1/4-n/indsun*(n-1)/(indsun-1)-m/indsun*n/(indsun-1)*2*1/2
```


```python
Mendelian(2,2,2)
```




    0.7833333333333333



# Translating RNA into Protein
[detailed introduction to this question](https://rosalind.info/problems/prot/)

The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.

The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.

Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).

Return: The protein string encoded by s.


## Sample
Sample Dataset:

AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

Sample Output:

MAMAPRTEINSTRING


```python
def Translate(RNA_seq):
    codonTable = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'', 'UAG':'',
    'UGC':'C', 'UGU':'C', 'UGA':'', 'UGG':'W',
    }
    proteinsequence = ''
    for n in range(0,len(RNA_seq),3):
        if RNA_seq[n:n+3] in codonTable.keys():
            proteinsequence += codonTable[RNA_seq[n:n+3]]
    return proteinsequence
```


```python
Translate("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")
```




    'MAMAPRTEINSTRING'



# Finding a Motif in DNA
[detailed introduction to this question](https://rosalind.info/problems/subs/)

Given two strings s and t, t is a substring of s if t is contained as a contiguous collection of symbols in s (as a result, t must be no longer than s).

The position of a symbol in a string is the total number of symbols found to its left, including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position i of s is denoted by s[i].

A substring of s can be represented as s[j:k], where j and k represent the starting and ending positions of the substring in s; for example, if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".

The location of a substring s[j:k] is its beginning position j; note that t will have multiple locations in s if it occurs more than once as a substring of s (see the Sample below).

Given: Two DNA strings s and t (each of length at most 1 kbp).

Return: All locations of t as a substring of s.

## Sample
Sample Dataset:

GATATATGCATATACTT

ATAT

Sample Output:

2 4 10

## Note
Different programming languages use different notations for positions of symbols in strings. Above, we use 1-based numbering, as opposed to 0-based numbering, which is used in Python. For s = "AUGCUUCAGAAAGGUCUUACG", 1-based numbering would state that s[1] = 'A' is the first symbol of the string, whereas this symbol is represented by s[0] in 0-based numbering. The idea of 0-based numbering propagates to substring indexing, so that s[2:5] becomes "GCUU" instead of "UGCU".

Note that in some programming languages, such as Python, s[j:k] returns only fragment from index j up to but not including index k, so that s[2:5] actually becomes "UGC", not "UGCU".


```python
def Motif_search(sequence,motif):
    index_list=[]
    for i in range(len(sequence)-len(motif)+1):
        if sequence[i:i+len(motif)]==motif:
            index_list.append(i+1)
    return index_list
```


```python
Motif_search("GATATATGCATATACTT","ATAT")
```




    [2, 4, 10]



# Consensus and Profile
[detailed introduction to this question](https://rosalind.info/problems/cons/)

A matrix is a rectangular table of values divided into rows and columns. An m×n matrix has m rows and n columns. Given a matrix A, we write Ai,j to indicate the value found at the intersection of row i and column j.

Say that we have a collection of DNA strings, all having the same length n. Their profile matrix is a 4×n matrix P in which P1,j represents the number of times that 'A' occurs in the jth position of one of the strings, P2,j represents the number of times that C occurs in the jth position, and so on (see below).

A consensus string c is a string of length n formed from our collection by taking the most common symbol at each position; the jth symbol of c therefore corresponds to the symbol having the maximum value in the j-th column of the profile matrix. Of course, there may be more than one most common symbol, leading to multiple possible consensus strings.

| | |
| :---: | :---: |
| | A T C C A G C T |
| |G G G C A A C T|
| |A T G G A T C T|
|DNA Strings|A A G C A A C C|
| |T T G G A A C T|
| |A T G C C A T T|
| |A T G G C A C T|

| | | |
| :---: | :---: | :---: |
| |A|5 1 0 0 5 5 0 0|
|Profile|C|0 0 1 4 2 0 6 1|
| |G|1 1 6 3 0 1 0 0|
| |T|1 5 0 0 0 1 1 6|
|Consensus| |A T G C A A C T|

Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)

## Sample
Sample Dataset:

\>Rosalind_1

ATCCAGCT

\>Rosalind_2

GGGCAACT

\>Rosalind_3

ATGGATCT

\>Rosalind_4

AAGCAACC

\>Rosalind_5

TTGGAACT

\>Rosalind_6

ATGCCATT

\>Rosalind_7

ATGGCACT

Sample Output:

ATGCAACT

A: 5 1 0 0 5 5 0 0

C: 0 0 1 4 2 0 6 1

G: 1 1 6 3 0 1 0 0

T: 1 5 0 0 0 1 1 6


```python
def read_fasta(file_path=""): # We will use this function frequently 
    """
    Loading FASTA file and return a iterative object
    """

    line = ""

    try:
        fasta_handle = open(file_path,"r")
    except:
        raise IOError("Your input FASTA file is not right!")

    # make sure the file is not empty
    while True:
        line = fasta_handle.readline()
        if line == "":
            return
        if line[0] == ">":
            break

    # when the file is not empty, we try to load FASTA file
    while True:
        if line[0] != ">":
            raise ValueError("Records in Fasta files should start with '>' character")
        title = line[1:].rstrip()
        lines = []
        line = fasta_handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = fasta_handle.readline()

        yield title,"".join(lines).replace(" ","").replace("\r","")

        if not line:
            return

    fasta_handle.close()
    assert False, "Your input FASTA file have format problem."
def Findanc(fa):
    seq_list=[]
    for key,value in read_fasta(fa):
        seq_list.append(value)
    names=locals()
    for bp in 'ATGC':
        bp_count=[]
        for i in range(len(seq_list[0])):
            col=[x[i] for x in seq_list]
            num=col.count(bp)
            bp_count.append(num)
        names[str(bp)+'_list']=bp_count
        print("{0}:{1}".format(bp,bp_count))
    anc_seq=""
    for i in range(len(seq_list[0])):
        q=0
        for bp in 'ATGC':
            if names[str(bp)+'_list'][i]>q:
                q=names[str(bp)+'_list'][i]
                base=bp
        anc_seq=anc_seq+base
    print(anc_seq)
```


```python
Findanc('profile.fa')
```

    A:[5, 1, 0, 0, 5, 5, 0, 0]
    T:[1, 5, 0, 0, 0, 1, 1, 6]
    G:[1, 1, 6, 3, 0, 1, 0, 0]
    C:[0, 0, 1, 4, 2, 0, 6, 1]
    ATGCAACT
    

# Mortal Fibonacci Rabbits
[detailed introduction to this question](https://rosalind.info/problems/fibd/)

Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation F<sub>n</sub>=F<sub>n−1</sub>+F<sub>n−2</sub> and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed number of months. See Figure 4 for a depiction of a rabbit tree in which rabbits live for three months (meaning that they reproduce only twice before dying).

Given: Positive integers n≤100 and m≤20.

Return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.

## Sample
Sample Dataset:

6 3

Sample Output:

4


```python
def Rabbit_Mortal(n,m):
    r=[0]*(n+1)
    r[1]=r[2]=1
    for i in range(3,n+1):
        r[i]=r[i-1]+r[i-2]
    rsum=r[n]
    for i in range(1,n-m+1):
        rsum=rsum-r[i]
    return rsum
```


```python
Rabbit_Mortal(6,3)
```




    4



# Overlap Graphs
[detailed introduction to this question](https://rosalind.info/problems/grph/)

A graph whose nodes have all been labeled can be represented by an adjacency list, in which each row of the list contains the two node labels corresponding to a unique edge.

A directed graph (or digraph) is a graph containing directed edges, each of which has an orientation. That is, a directed edge is represented by an arrow instead of a line segment; the starting and ending nodes of an edge form its tail and head, respectively. The directed edge with tail v and head w is represented by (v,w) (but not by (w,v)). A directed loop is a directed edge of the form (v,v).

For a collection of strings and a positive integer k, the overlap graph for the strings is a directed graph O<sub>k</sub> in which each string is represented by a node, and string s is connected to string t with a directed edge when there is a length k suffix of s that matches a length k prefix of t, as long as s≠t; we demand s≠t to prevent directed loops in the overlap graph (although directed cycles may be present).

Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.

Return: The adjacency list corresponding to O<sub>3</sub>. You may return edges in any order.

## Sample
Sample Dataset:

\>Rosalind_0498

AAATAAA

\>Rosalind_2391

AAATTTT

\>Rosalind_2323

TTTTCCC

\>Rosalind_0442

AAATCCC

\>Rosalind_5013

GGGTGGG

Sample Output:

Rosalind_0498 Rosalind_2391

Rosalind_0498 Rosalind_0442

Rosalind_2391 Rosalind_2323


```python
def read_fa(fa):
    dic_fa={}    
    with open(fa) as f:
        for line in f:
            line=line.strip()
            if line[0]=='>':
                header=line[1:]
                dic_fa[header]=""
            else:
                seq=line
                dic_fa[header]=seq
    return dic_fa
def Find_overlap(f,k):
    dic_fasta=read_fa(f)
    for key in dic_fasta.keys():
        for key2 in dic_fasta.keys():
                if key != key2 and dic_fasta[key][-k:] == dic_fasta[key2][:k]:
                        print(key+'\t'+key2)
```


```python
Find_overlap("overlap.fa",3)
```

    Rosalind_0498	Rosalind_2391
    Rosalind_0498	Rosalind_0442
    Rosalind_2391	Rosalind_2323
    

# Calculating Expected Offspring
[detailed introduction to this question](https://rosalind.info/problems/iev/)

For a random variable X taking integer values between 1 and n, the expected value of X is $E(X)=\sum_{k=1}^n k×Pr(X=k)$. The expected value offers us a way of taking the long-term average of a random variable over a large number of trials.

As a motivating example, let X be the number on a six-sided die. Over a large number of rolls, we should expect to obtain an average of 3.5 on the die (even though it's not possible to roll a 3.5). The formula for expected value confirms that $E(X)=\sum_{k=1}^6 k×Pr(X=k)=3.5$.

More generally, a random variable for which every one of a number of equally spaced outcomes has the same probability is called a uniform random variable (in the die example, this "equal spacing" is equal to 1). We can generalize our die example to find that if X is a uniform random variable with minimum possible value a and maximum possible value b, then $E(X)=\frac{a+b}{2}$. You may also wish to verify that for the dice example, if Y is the random variable associated with the outcome of a second die roll, then $E(X+Y)=7$.

Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor. In order, the six given integers represent the number of couples having the following genotypes:

1.AA-AA

2.AA-Aa

3.AA-aa

4.Aa-Aa

5.Aa-aa

6.aa-aa

Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.

## Sample
Sample Dataset:

1 0 0 1 0 1

Sample Output:

3.5


```python
def Count_do(a,b,c,d,e,f):
    return a*2+b*2+c*2+d*1.5+e*0.5
```


```python
Count_do(1,0,0,1,0,1)
```




    3.5



# Finding a Shared Motif
[detailed introduction to this question](https://rosalind.info/problems/lcsm/)

A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGTATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".

Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.

Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)

## Sample
Sample Dataset:

\>Rosalind_1

GATTACA

\>Rosalind_2

TAGACCA

\>Rosalind_3

ATACA

Sample Output:

AC


```python
def read_fa(fa):
    dic_fa={}    
    with open(fa) as f:
        for line in f:
            line=line.strip()
            if line[0]=='>':
                header=line
                dic_fa[header]=""
            else:
                seq=line
                dic_fa[header]=seq
    return dic_fa
def Find_sharedmotif(fa):
    dict_fa={}
    shortest=""
    shortestlen=100
    dict_fa=read_fa(fa)
    for key in dict_fa.keys():
        if len(dict_fa[key])<shortestlen:
            shortest=dict_fa[key]
            shortestlen=len(dict_fa[key])
    motif=""
    paired=0
    for i in range(shortestlen):
        j=0
        while j<=i:
            motif=shortest[j:j+shortestlen-i]
            for key in dict_fa.keys():
                if motif in dict_fa[key]:
                    paired=paired+1
            if paired==len(dict_fa.keys()):
                return motif
            else:
                paired=0
                j=j+1
```


```python
Find_sharedmotif("motif1.fa")
```




    'TA'



# Independent Alleles
[detailed introduction to this question](https://rosalind.info/problems/lia/)

Two events A and B are independent if Pr(A and B) is equal to Pr(A)×Pr(B). In other words, the events do not influence each other, so that we may simply calculate each of the individual probabilities separately and then multiply.

More generally, random variables X and Y are independent if whenever A and B are respective events for X and Y, A and B are independent (i.e., Pr(A and B)=Pr(A)×Pr(B)).

As an example of how helpful independence can be for calculating probabilities, let X and Y represent the numbers showing on two six-sided dice. Intuitively, the number of pips showing on one die should not affect the number showing on the other die. If we want to find the probability that X+Y is odd, then we don't need to draw a tree diagram and consider all possibilities. We simply first note that for X+Y to be odd, either X is even and Y is odd or X is odd and Y is even. In terms of probability, Pr(X+Y is odd)=Pr(X is even and Y is odd)+Pr(X is odd and Y is even). Using independence, this becomes \[Pr(X is even)×Pr(Y is odd)\]+\[Pr(X is odd)×Pr(Y is even)\], or $(\frac12)^2+(\frac12)^2=\frac12$.

Given: Two positive integers k (k≤7) and N (N≤2<sup>k</sup>). In this problem, we begin with Tom, who in the 0th generation has genotype Aa Bb. Tom has two children in the 1st generation, each of whom has two children, and so on. Each organism always mates with an organism having genotype Aa Bb.

Return: The probability that at least N Aa Bb organisms will belong to the k-th generation of Tom's family tree (don't count the Aa Bb mates at each level). Assume that Mendel's second law holds for the factors.

## Sample
Sample Dataset:

2 1

Sample Output:

0.684


```python
import math
def Com(a,b):
    return math.factorial(a)/(math.factorial(b)*math.factorial(a-b))
def Prob_AaBb(k,N):
    children_sum=2**k
    AaBb_pro=0
    for i in range(0,N):
        AaBb_pro=AaBb_pro+((1/4)**i)*((3/4)**(children_sum-i))*Com(children_sum,i)
    return 1-AaBb_pro
```


```python
Prob_AaBb(2,1)
```




    0.68359375



# Finding a Protein Motif
[detailed introduction to this question](https://rosalind.info/problems/mprt/)

To allow for the presence of its varying forms, a protein motif is represented by a shorthand as follows: [XY] means "either X or Y" and {X} means "any amino acid except X." For example, the N-glycosylation motif is written as N{P}[ST]{P}.

You can see the complete description and features of a particular protein by its access ID "uniprot_id" in the UniProt database, by inserting the ID number into

http://www.uniprot.org/uniprot/uniprot_id

Alternatively, you can obtain a protein sequence in FASTA format by following

http://www.uniprot.org/uniprot/uniprot_id.fasta

For example, the data for protein B5ZC00 can be found at http://www.uniprot.org/uniprot/B5ZC00.

Given: At most 15 UniProt Protein Database access IDs.

Return: For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of locations in the protein string where the motif can be found.

## Sample
Sample Dataset:

A2Z669

B5ZC00

P07204_TRBM_HUMAN

P20840_SAG1_YEAST

Sample Output:

B5ZC00

85 118 142 306 395

P07204_TRBM_HUMAN

47 115 116 382 409

P20840_SAG1_YEAST

79 109 135 248 306 348 364 402 485 501 614


```python
import urllib
from urllib.request import urlopen
import re
def Find_pro_motif(pro_list):
    for name in list:
        url = 'http://www.uniprot.org/uniprot/'+name+'.fasta'
        req = urllib.request.Request(url)
        response = urlopen(req)
        the_page = response.read()
        start = the_page.find(b'\n')
        seq = the_page[start+1:].replace(b'\n',b'')
        seq = str(seq)
        regex = re.compile(r'N(?=[^P][ST][^P])')
        index = 0
        out=[]
        while(index<len(seq)):
            if re.search(regex,seq[index:]) == None:
                break
            if re.match(regex,seq[index:]) != None:
                out.append(index-1)
            index += 1
        if out != []:
            print (name)
            print (' '.join([ str(i) for i in out]))
```


```python
pro_list = ['A2Z669','B5ZC00','P07204_TRBM_HUMAN','P20840_SAG1_YEAST']
Find_pro_motif(pro_list)
```

    B5ZC00
    85 118 142 306 395
    P07204_TRBM_HUMAN
    47 115 116 382 409
    P20840_SAG1_YEAST
    79 109 135 248 306 348 364 402 485 501 614
    

# Inferring mRNA from Protein
[detailed introduction to this question](https://rosalind.info/problems/mrna/)

For positive integers a and n, a modulo n (written a\mathchoice mod n in shorthand) is the remainder when a is divided by n. For example, 29\mathchoice mod 11=7 because 29=11×2+7.

Modular arithmetic is the study of addition, subtraction, multiplication, and division with respect to the modulo operation. We say that a and b are congruent modulo n if a\mathchoice mod n=b\mathchoice mod n; in this case, we use the notation a≡b\mathchoice mod n.

Two useful facts in modular arithmetic are that if a≡b\mathchoice mod n and c≡d\mathchoice mod n, then a+c≡b+d\mathchoice mod n and a×c≡b×d\mathchoice mod n. To check your understanding of these rules, you may wish to verify these relationships for a=29, b=73, c=10, d=32, and n=11.

As you will see in this exercise, some Rosalind problems will ask for a (very large) integer solution modulo a smaller number to avoid the computational pitfalls that arise with storing such large numbers.

Given: A protein string of length at most 1000 aa.

Return: The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)

## Sample
Sample Dataset:

MA

Sample Output:

12


```python
translate_table={
    "UUU":"F", "CUU":"L", "AUU":"I", "GUU":"V",
    "UUC":"F", "CUC":"L", "AUC":"I", "GUC":"V",
    "UUA":"L", "CUA":"L", "AUA":"I", "GUA":"V",
    "UUG":"L", "CUG":"L", "AUG":"M", "GUG":"V",
    "UCU":"S", "CCU":"P", "ACU":"T", "GCU":"A",
    "UCC":"S", "CCC":"P", "ACC":"T", "GCC":"A",
    "UCA":"S", "CCA":"P", "ACA":"T", "GCA":"A",
    "UCG":"S", "CCG":"P", "ACG":"T", "GCG":"A",
    "UAU":"Y", "CAU":"H", "AAU":"N", "GAU":"D",
    "UAC":"Y", "CAC":"H", "AAC":"N", "GAC":"D",
    "UAA":"*",  "CAA":"Q", "AAA":"K", "GAA":"E",
    "UAG":"*",  "CAG":"Q", "AAG":"K", "GAG":"E",
    "UGU":"C", "CGU":"R", "AGU":"S", "GGU":"G",
    "UGC":"C", "CGC":"R", "AGC":"S", "GGC":"G",
    "UGA":"*",  "CGA":"R", "AGA":"R", "GGA":"G",
    "UGG":"W", "CGG":"R", "AGG":"R", "GGG":"G"
}
translate_count_table={}
for codon,aa in translate_table.items():
    if aa =="*":
        pass
    elif aa not in translate_count_table.keys():
        translate_count_table[aa] = 1
    else:
        translate_count_table[aa] += 1
def Count_mRNA(protein):
    rna_count=1
    for aa in protein:
        rna_count *= translate_count_table[aa]
    return 3*rna_count
```


```python
Count_mRNA("MA")
```




    12



# Open Reading Frames
[detailed introduction to this question](https://rosalind.info/problems/orf/)

Either strand of a DNA double helix can serve as the coding strand for RNA transcription. Hence, a given DNA string implies six total reading frames, or ways in which the same region of DNA can be translated into amino acids: three reading frames result from reading the string itself, whereas three more result from reading its reverse complement.

An open reading frame (ORF) is one which starts from the start codon and ends by stop codon, without any other stop codons in between. Thus, a candidate protein string is derived by translating an open reading frame into amino acids until a stop codon is reached.

Given: A DNA string s of length at most 1 kbp.

Return: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in any order.

## Sample
Sample Dataset:

AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG

Sample Output:

MLLGSFRLIPKETLIQVAGSSPCNLS

M

MGMTPRLGLESLLE

MTPRLGLESLLE


```python
codonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
    'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
}
def Reverse_DNA(dna_seq):
    return dna_seq.upper().replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()[::-1]
def DNA_to_pro(dna_seq):
    Pro_set=set([])
    i=0    
    while i<len(dna_seq):
        pro_seq=''
        if dna_seq[i:i+3]=='ATG':
            pro_seq+='M'
            j=i+3
            while j<=len(dna_seq)-3:
                if dna_seq[j:j+3] in ['TGA','TAA','TAG']:
                    Pro_set.add(pro_seq)
                    break
                else:
                    pro_seq+=codonTable[dna_seq[j:j+3]]
                    j+=3
            i+=1
        else:
            i+=1     
    return Pro_set
def DNA_ORF_pro(dna_seq):
    return DNA_to_pro(dna_seq)|DNA_to_pro(Reverse_DNA(dna_seq))
```


```python
DNA_ORF_pro('AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG')
```




    {'M', 'MGMTPRLGLESLLE', 'MLLGSFRLIPKETLIQVAGSSPCNLS', 'MTPRLGLESLLE'}



# Enumerating Gene Orders
[detailed introduction to this question](https://rosalind.info/problems/perm/)

A permutation of length n is an ordering of the positive integers {1,2,…,n}. For example, π=(5,3,2,1,4) is a permutation of length 5.

Given: A positive integer n≤7.

Return: The total number of permutations of length n, followed by a list of all such permutations (in any order).

## Sample
Sample Dataset:

3

Sample Output:

6

1 2 3

1 3 2

2 1 3

2 3 1

3 1 2

3 2 1


```python
import itertools
import math
def Perm(n):
    print(math.factorial(n))
    num_list = []
    for i in range(1,n+1):
        num_list.append(i)
    perm_list=list(itertools.permutations(num_list,n))
    for i in perm_list:
        print(" ".join(map(str,i)))
```


```python
Perm(3)
```

    6
    1 2 3
    1 3 2
    2 1 3
    2 3 1
    3 1 2
    3 2 1
    
