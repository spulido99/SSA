# SSA
Lighweigth framework for Small Subnetwork Analysis (SSA)

## How to run SSA.ME

* Download and extract the software from [this link](http://bioinformatics.intec.ugent.be/ssame/SSA.zip)
* Running SSA.ME consist of 2 steps:
 

### Generate the binary matrix input file

Mutual exclusivity tool run using a genomic alteration matrix that is Zero is a gene is not altered in a sample and One if the gene is altered in that sample. Alterations can represent anything: single nucleotide mutations, INDELS, deletions, amplifications, methilation, etc.

```
java -jar SSA.jar ME_input -m <maf file> -e cnv_peaks=<file>,exp=<file>,cnv_thresholds=<file>
Usage: SSA.ME. [options]

  -o <value> | --outputPrefix <value>
        The name to be used in the output files (XXX.m2 and XXX.glst).
  -s <value> | --seedGenesMutations <value>
        The number of mutated samples required in a gene to be included (default: 1)
  -m <value> | --maf <value>
        Mutation .maf file
  -e <value> | --expression <value>
        expression file and GISTIC files (... -e cnv_peaks=<file1>,exp=<file2>,cnv_thresholds=<file3> ...).
```
We provide a binary alteration matrix generation sub-tool that uses somatic MAF files and expression profiles with amplication/deletion information (GISTIC output) to generate this matrix. 

To generate the binary alteration matrix, for example, download the publicly available [TCGA Breast Cancer data from 2012](https://tcga-data.nci.nih.gov/docs/publications/brca_2012/) .maf file [Somatic MAF archive](http://tcga-data.nci.nih.gov/docs/publications/brca_2012/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.tar.gz) and run SSA.ME input creation step:
```
java -jar SSA.jar ME_input -m genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf
```

It will create 2 files called SSAME_input.m2 and SSAME_input.glst. The m2 represent a binary matrix file containing which samples have each gene mutated. The gene list is a of all the mutated genes. This 2 files follow the *de facto* standard followed by other mutual exclusivity tools and can be used with them.

### Run SSA.ME

SSA.ME uses the .m2 file as input for the Mutual Exclusivity (ME) step.

```
java -jar SSA.jar ME -m SSAME_input.m2 -s 5 -i 10000
```

### Visualizing the Output

The output contain 4 files:

+ __ME_network.html__ : An interactive html page showing the selected network. [See example](http://bioinformatics.intec.ugent.be/ssame/ME_network.html)
+ __ME_edges__ : The interactions between the genes selected forming mutual exclusivity
+ __ME_nodes__ : The genes in the network (nodes) with additional information as the convergence iteration and the best observed Small Subnetwork detected for that gene.
+ __ME_pattern__ : A matrix showing for each gene and sample if there was a mutation present.

## Solve other problems using SSA

Do you want to program new applications using SSA?

* Checkout the ''development'' branch.
* Install [SBT](http://www.scala-sbt.org/0.13/tutorial/Setup.html)
* To program using Eclipse, run the ''[eclipse](https://github.com/typesafehub/sbteclipse/wiki/Using-sbteclipse)'' sbt command
* Follow the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) to create your new biological applications.

