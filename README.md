This repository contains simulator used in the work "NestedBD: Bayesian Inference of Phylogenetic Trees From Single-Cell Copy Number Profile Data Under a Birth-Death Model". The simulator is based on https://github.com/compbiofan/SingleCellCNABenchmark, with several modifications. 

## Usage 

For the purpose of the study, we set X=8, W=0, m=20M and e=100M. The value of c under different simulation setting is described in the method section. 

### Step 1

python main.par.py -r $dir -n $n -X $X -t $ref -W $W -C $C -m $m -e $e -amp $amp -SP $SP

- $dir: the folder where the simulated data will be put. It could be a relative path. For example: large_dataset/. Default: test.

- $n: number of cells in the tree.

- $X: how much more CNAs on the edge to the root than other edges. For example, 8.

- $ref (required): reference fasta file in an absolute path.

- $W: if there are whole chromosomal amplifications, 1 (yes) or 0 (no).

- $C: the probability that a chromosome may be amplified if $W is 1.

- $m: minimum copy number size.

- $e: parameter p in exponential distribution for the copy number size that will be added to $m.


##### New parameters

- $amp: control the non-uniformness when selecting position of CNA on genome during simulation. When setting to 0 the CNAs are sampled randomly. Default is 0.
 
- $SP (required): Tree option. When setting SP = 0, a prompt will require user enter the path to a file contain a tree in newick format to be used for simulation; SP = 1 simulate a tree with birth death process before adding CNAs along the branches; SP = 2 simulate a beta splitting tree and CNA along branches simulataneously. 

- $c: event multipler, the number of events added to each branch is sampled from a poisson distribution with mean eqaul to the product of %c and branch length.

### Step 2
Same as the original simulator.
