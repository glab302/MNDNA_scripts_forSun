
### peak calling by same sample
for((s=1;s<=10;s++))
do
    macs2 callpeak -t ../01.rbcdna/rbcDNA_GLRHD000${s}.60m.nodup.q30.bam -c ../02.gdna/gDNA_GLGHD000${s}.60m.nodup.q30.bam -f BAMPE --nomodel --broad -g hs -n broad_t_rbcDNA_GLRHD000${s}.60m.nodup.q30.bam.c_gDNA_GLGHD000${s}.60m.nodup.q30.bam
done

### merge results
cat broad_*_peaks.broadPeak > merge_t_rbcDNA_c_gDNA.10samples.broadPeak.bed
sort -k1,1 -k2,2n merge_t_rbcDNA_c_gDNA.10samples.broadPeak.bed > merge_t_rbcDNA_c_gDNA.10samples.sorted.broadPeak.bed

echo -en 'chromosome\tstart\tend\tcount\tfoldenrichment\n' > result_d_0_merge_t_rbcDNA_c_gDNA.10samples.sorted_count_mean.broadPeak.bed
bedtools merge -i merge_t_rbcDNA_c_gDNA.10samples.sorted.broadPeak.bed -d 0 -c 1,7 -o count,mean >> result_d_0_merge_t_rbcDNA_c_gDNA.10samples.sorted_count_mean.broadPeak.bed
cat result_d_0_merge_t_rbcDNA_c_gDNA.10samples.sorted_count_mean.broadPeak.bed|awk -F '\t' '$1!="X"&&$1!="Y"&&$1!="M"{print $0}' > result_d_0_merge_t_rbcDNA_c_gDNA.10samples.broadPeak.bed

# echo -en 'chromosome\tstart\tend\tcount\tfoldenrichment\n' > result_d_1000_merge_t_rbcDNA_c_gDNA.10samples.sorted_count_mean.broadPeak.bed
# bedtools merge -i merge_t_rbcDNA_c_gDNA.10samples.sorted.broadPeak.bed -d 1000 -c 1,7 -o count,mean >> result_d_1000_merge_t_rbcDNA_c_gDNA.10samples.sorted_count_mean.broadPeak.bed
# cat result_d_1000_merge_t_rbcDNA_c_gDNA.10samples.sorted_count_mean.broadPeak.bed|awk -F '\t' '$1!="X"&&$1!="Y"&&$1!="M"{print $0}' > result_d_1000_merge_t_rbcDNA_c_gDNA.10samples.broadPeak.bed

