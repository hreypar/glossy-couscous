
# this would be the command line
for i in {1..5}; do java -Xmx32G -jar /home/hreyes/bin/ARACNe-AP/dist/aracne.jar -e data/counts_matrix/healthy/tcga-brca-lumA-healthy.txt --tfs data/counts_matrix/healthy/tcga-brca-lumA-healthy_features.txt -o results/aracne_networks/healthy/ --pvalue 1E-8 --seed $i; done

