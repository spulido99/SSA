#!/bin/bash
for cancer in BLCA BRCA COADREAD GBM HNSC KIRC LAML LUAD LUSC OV UCEC STAD
do
	mkdir $cancer
	cd $cancer
	wget http://gdac.broadinstitute.org/runs/analyses__2015_08_21/data/${cancer}-TP/20150821/gdac.broadinstitute.org_${cancer}-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0.tar.gz 
	wget http://gdac.broadinstitute.org/runs/analyses__2015_08_21/data/${cancer}-TP/20150821/gdac.broadinstitute.org_${cancer}-TP.Correlate_CopyNumber_vs_mRNA.Level_4.2015082100.0.0.tar.gz
	wget http://gdac.broadinstitute.org/runs/analyses__2015_08_21/data/${cancer}-TP/20150821/gdac.broadinstitute.org_${cancer}-TP.Mutation_Assessor.Level_4.2015082100.0.0.tar.gz

	tar xf gdac.broadinstitute.org_${cancer}-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0.tar.gz
	tar xf gdac.broadinstitute.org_${cancer}-TP.Correlate_CopyNumber_vs_mRNA.Level_4.2015082100.0.0.tar.gz
	tar xf gdac.broadinstitute.org_${cancer}-TP.Mutation_Assessor.Level_4.2015082100.0.0.tar.gz

	rm gdac.broadinstitute.org_${cancer}-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0.tar.gz
	rm gdac.broadinstitute.org_${cancer}-TP.Correlate_CopyNumber_vs_mRNA.Level_4.2015082100.0.0.tar.gz
	rm gdac.broadinstitute.org_${cancer}-TP.Mutation_Assessor.Level_4.2015082100.0.0.tar.gz
	cd ..
done
