
# Rscript run_addanno.r

feature_list="MNDNApos6566 CRC0106 PanC0106"
for s in $feature_list
do
    cat result/GenomicAnnotationIn_${s}_addtoRegion.log|cut -f1 > ${s}.bed
    cp ${s}.bed query.list
    python3 annotation.py
    python3 group_annotation.py
    mv output.txt result/${s}.output.txt
    mv grouped_output.txt result/${s}.grouped_output.txt
    rm -rf query.list ${s}.bed
done
