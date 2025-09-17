declare -a model_li=("B" "CMP.LMPP" "CD14.Mono.1" "Early.Baso" "CD14.Mono.2" "Early.Eryth" "CD4.M" "GMP" "CD4.N1"
                     "GMP.Neut" "CD4.N2" "HSC" "CD8.CM" "Late.Eryth" "CD8.EM" "NK" "CD8.N" "pDC" "cDC"
                     "Plasma" "CLP.1" "Pre.B" "CLP.2" "Unk_26")
declare -a target_li=("B" "CMP.LMPP" "CD14.Mono.1" "Early.Baso" "CD14.Mono.2" "Early.Eryth" "CD4.M" "GMP" "CD4.N1"
                      "GMP.Neut" "CD4.N2" "HSC" "CD8.CM" "Late.Eryth" "CD8.EM" "NK" "CD8.N" "pDC" "cDC"
                      "Plasma" "CLP.1" "Pre.B" "CLP.2" "Unk_26")

for model in ${!model_li[@]}  
do
for target in ${!target_li[@]}
do
if [[ "${model_li[$model]}" != "${target_li[$target]}" ]]; then
Rscript blood_predict_other.R  "${model_li[$model]}" "${target_li[$target]}"
fi
done
done
