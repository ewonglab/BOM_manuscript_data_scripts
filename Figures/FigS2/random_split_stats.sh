my_seeds=(2549889 2221565  856018 3812659 1343338)
for my_seed in ${my_seeds[@]}; do
for celltype in ${celltypes[@]}; do
Rscript save_stats.R ./random_samples/${celltype}*_seed${my_seed}_pred.txt \
./random_samples/BOM_binary_seed${my_seed}_stats
done
done
