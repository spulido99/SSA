wget http://gdac.broadinstitute.org/runs/analyses__2015_08_21/data/BLCA-TP/20150821/gdac.broadinstitute.org_BLCA-TP.Correlate_CopyNumber_vs_mRNAseq.Level_4.2015082100.0.0.tar.gz
tar xf gdac.broadinstitute.org_BLCA-TP.Correlate_CopyNumber_vs_mRNAseq.Level_4.2015082100.0.0.tar.gz
rm gdac.broadinstitute.org_BLCA-TP.Correlate_CopyNumber_vs_mRNAseq.Level_4.2015082100.0.0.tar.gz

mkdir gdac.broadinstitute.org_BLCA-TP.Correlate_CopyNumber_vs_mRNA.Level_4.2015082100.0.0.tar.gz

awk '{print $2,"\t",$4,"\t",$5,"\t",$6}' gdac.broadinstitute.org_BLCA-TP.Correlate_CopyNumber_vs_mRNAseq.Level_4.2015082100.0.0/BLCA-TPcors.txt > gdac.broadinstitute.org_BLCA-TP.Correlate_CopyNumber_vs_mRNA.Level_4.2015082100.0.0/BLCA-TP.CORS.tsv
