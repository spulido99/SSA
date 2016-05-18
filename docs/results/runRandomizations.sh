java -Xmx10g -jar SSA.jar ME_pval -m BRCA -o SSAME -i 500 -r 0.002 -f 0.998 -p 100 -s 1 --randomDistance 50000

java -Xmx10g -jar SSA.jar ME_pval -m BRCA -o SSAME_MEMo -i 500 -r 0.002 -f 0.998 -p 100 -s 1 --randomDistance 50000 -n MEMo/hrn2 

java -Xmx30g -jar SSA.jar ME_pval -m BRCA.s10 -n MEMo/hrn2 -o SSAME_MEMo_sameinput -i 50000 -r 0.002 -f 0.998 -p 100 -s 10 -g 10 --randomDistance 50000

java -Xmx10g -jar SSA.jar ME_pval -m GBM -o SSAME_GBM -i 500 -r 0.002 -f 0.998 -p 100 -s 1 --randomDistance 50000

java -Xmx10g -jar SSA.jar ME_pval -m OV -o SSAME_OV -i 500 -r 0.002 -f 0.998 -p 100 -s 1 --randomDistance 50000

