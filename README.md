# Edison
Given an assembly of an organism with a known reference genome, how accuracy is any given scaffolding of that assembly? While the notion of accuracy is elusive in genome assembly, we attempt to mathematically define it in this context. We developed Edison (<ins>E</ins>dit <ins>Di</ins>stance <ins>S</ins>caff<ins>o</ins>ldi<ins>n</ins>g) to compute this scaffolding accuracy. 

We implement the following functions - 

1) **Edit Distance** - The number of edits needed to alter a scaffolding such that it would most closely resemble a reference genome. 
2) **Accuracy** - Length weighted version of edit distance. It represents the number of bases pairs that would not need to move during the edit process. 
3) **Grouping** - Length weighted accuracy of clustering contigs into chromosomes. 
4) **Ordering** - Length weighted accuracy of contig adjacencies without regard to orientation. 
5) **Orientation** - Length weighted accuracy of contig adjacencies with regard to orientation. 

# Overview
We start with a reference FASTA file and scaffolded assembly FASTA file. The scaffolded assembly is split at any location with more than a specified number of Ns in a row (Default: 10) and an AGP file is created to annotate the positions of these inferred contigs. These contigs are then mapped using Mummer to assign them coordinates to the best matching reference chromosome. The mapping is used to create an AGP file which details where the contigs _ought_ to go. Finally the AGP files are compared to procude the scaffolding accuracy metrics. 
 
## Double Cut and Join Edit Distance
During the late 90s and early 2000s, scientists were attempting to mathematically formulate the most parsimonious explanations for genomic evolution, equivalently known as the edit distance between two genomes. In 2005, Yancopoulos et al. described one of the more straight forward renditions of the edit function - the Double Cut and Join (DCJ) model. In 2006, Bergeron et al. provided an elegant simplification to the DCJ model by utilizing an adjacency graph, and it is this formulation we follow in this work. 

## Edit Distance
Let us examine the central equation in Bergeron et al., - 

d = N - (C + I/2)

Where d is the edit distance, N is the number of contigs, C is the number of cycles formed in the adjacency graph, and I is the number of odd cycles. 

## Accuracy
When two scaffoldings are equivalent, the edit distance is 0 and the number of cycles plus half the number of odd cycles is equal to the number of contigs. Thus contributions reducing edit distance only come from cycles and odd paths. A length weighted version of these contributions is our notion of accuracy. If the length of an adjacency is defined as the sum of the lengths comprising that adjacency, then the length of a cycle or path is simply the maximum length of its component adjacencies. 

## Grouping 
For each chromosome in the reference, we find the scaffold which maximizes the length weighted Jaccard index. Then we take a length weighted average of these values.

## Ordering
We enumerate the adjacencies without regard to orientation in the scaffolding and the reference and then calculated the length weighted percent of scaffolding adjacencies found in the reference. 

## Orientation
Calculated in a similar fashion as the ordering metric, but also keeping track of orientation of contigs. 

# Running
Install [Mummer](https://github.com/mummer4/mummer) and have it available in `$PATH`. 

```
python edit_distance.py \
    -a assembly.fasta \
    -r reference.fasta
```
## Output

```
Generating contigs from scaffolds.
Running Mummer.

     Contig Statistics     
-----------------------------
            Starting:    121
  Missing Alignments:     10
      Low Alignments:      0
   Missing Reference:      0
   Minimum Alignment: 41.48%
   Average Alignment: 81.01%

     Assembly Statistics     
-----------------------------
            Grouping: 71.64%
            Ordering: 78.60%
         Orientation: 74.55%
            Accuracy: 77.46%
       Edit Distance:     29
```
