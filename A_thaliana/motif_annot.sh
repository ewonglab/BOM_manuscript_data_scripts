motifsdb='Arabidopsis_thaliana.meme'

for fa_file in *_distal_trim600bp.fa
do
base_name=$(basename -s _distal_trim600bp.fa ${fa_file})
/g/data/zk16/cc3704/meme/meme-5.4.0/src/fimo --thresh 0.0001 \
--o trimmed600bp/atha_CisBP2.0/${base_name} ${motifsdb} \
${fa_file}
done
