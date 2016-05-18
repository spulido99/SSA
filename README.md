# SSA
Lighweigth framework for Small Subnetwork Analysis (SSA)

## How to run SSA-ME

* Download and extract the software from [this link](http://bioinformatics.intec.ugent.be/ssame/SSA.zip)
* Running SSA.ME consist of 2 steps:
 

### Generate the binary matrix input file

Mutual exclusivity tools run using a genomic alteration matrix that is Zero if a gene is not altered in a sample and One if the gene is altered in that sample. Alterations can represent anything: single nucleotide mutations, INDELS, deletions, amplifications, methilation, etc.

We provide a binary alteration matrix generation tool that uses somatic MAF files and expression profiles with amplication/deletion information (GISTIC output) to generate this matrix. 

To generate the binary alteration matrix, for example, download and decompress the publicly available [Firehose](http://firebrowse.org/?cohort=BRCA) BRCA Data:
```
wget http://gdac.broadinstitute.org/runs/analyses__2015_08_21/data/BRCA-TP/20150821/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0.tar.gz 
wget http://gdac.broadinstitute.org/runs/analyses__2015_08_21/data/BRCA-TP/20150821/gdac.broadinstitute.org_BRCA-TP.Correlate_CopyNumber_vs_mRNA.Level_4.2015082100.0.0.tar.gz
wget http://gdac.broadinstitute.org/runs/analyses__2015_08_21/data/BRCA-TP/20150821/gdac.broadinstitute.org_BRCA-TP.Mutation_Assessor.Level_4.2015082100.0.0.tar.gz

tar xf gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0.tar.gz
tar xf gdac.broadinstitute.org_BRCA-TP.Correlate_CopyNumber_vs_mRNA.Level_4.2015082100.0.0.tar.gz
tar xf gdac.broadinstitute.org_BRCA-TP.Mutation_Assessor.Level_4.2015082100.0.0.tar.gz

```

and run SSA.ME input creation step:
```
java -jar SSA.jar ME_input -m gdac.broadinstitute.org_BRCA-TP.Mutation_Assessor.Level_4.2015082100.0.0/BRCA-TP.maf.annotated -e corr=gdac.broadinstitute.org_BRCA-TP.Correlate_CopyNumber_vs_mRNA.Level_4.2015082100.0.0/BRCA-TP.CORS.tsv,gistic=gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0/ -o BRCA
```

It will create several files called BRCA.m2, KIRC.tbs, BRCA.glst, BRCA.byGene.stats and BRCA.bySample.stats. The .m2 and .tbs file represent a binary matrix containing which samples have which genes mutated (the .tbs file can be used directly in [Gitools](http://www.gitools.org/)). The gene list is a of all the mutated genes present in the matrix. The .stats file show the number of mutations by gene or by sample in the dataset.

### Run SSA-ME

SSA-ME uses the .m2 file as the input for the Mutual Exclusivity (ME) step.

```
java -Xmx30g -jar SSA.jar ME -m BRCA -o SSAME_BRCA -i 5000 -r 0.0002 -f 0.9998 -p 200 -s 3 --processors 60 -n HT,hiII14,reactome
```
(change the number of processors to those you have available)

## Run statistical analysis (Bootstraap)

SSA-ME provides a statistical analysis based on bootstraap to select only those genes supported by random sampling with replacement from the data.

```
java -Xmx30g -jar SSA.jar ME_btstrp -m BRCA -o SSAME_BRCA -i 500 -r 0.0002 -f 0.9998 -p 200 -s 3 --processors 60 -n HT,hiII14,reactome --bootstraapExperiments 1000 --useNCG false
```

### Visualizing the Output

The output contain 4 files (twice, one set from the original run, a second set from the Bootstraaping marked as .selected):

+ __network.html__ : An interactive html page showing the selected network. [See example](http://bioinformatics.intec.ugent.be/ssame/ME_network.html)
+ __edges__ : The interactions between the genes selected forming mutual exclusivity
+ __nodes__ : The genes in the network (nodes) with additional information as the convergence iteration and the best observed Small Subnetwork detected for that gene.
+ __pattern__ : A matrix showing for each gene and sample if there was a mutation present.

## Solve other problems using SSA

Do you want to program new applications using SSA?

* Checkout using any git tool (e.g. [SourceTree](https://www.sourcetreeapp.com/), [GitHub Desktop](https://desktop.github.com/)) the ''develop'' branch.
* Install [SBT](http://www.scala-sbt.org/0.13/tutorial/Setup.html)
* To program using Eclipse, run the ''[eclipse](https://github.com/typesafehub/sbteclipse/wiki/Using-sbteclipse)'' sbt command
* Follow the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) to create your new biological applications.

