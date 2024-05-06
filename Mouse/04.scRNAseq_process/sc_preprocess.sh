#!/bin/bash
#SBATCH -J BM4
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -q normal
#SBATCH -c 40
#SBATCH --mem=400G

module load cellranger/6.1.1
module load R/4.2.1
module load gcc/11.2.0
module load fftw/3.3.10

### 1.preprocess
### for only single-cell RNA-seq:
rawdata_path='Data/CleanData/'
sample="A002 A004 A005 A006 W001 W002 W004 W005 W006"
for s in $sample
do
cellranger count --id=$s --fastqs=$rawdata_path/$s --sample=$s --transcriptome=refdata-gex-mm10-2020-A/ --localcores=40
done
### for single-cell RNA-seq and TCR-seq:
# cellranger multi --id=D01_multi_result --csv=D01_multi-config-template.csv --localcores=40

### 2.quality control (single sample)
R_script="Rscript/"
sample="A002 A004 A005 A006 W001 W002 W004 W005 W006"
for s in $sample
do
subpath="${s}/outs/filtered_feature_bc_matrix/"
mkdir $s
Rscript ${R_script}/1.Seurat3_QualityControl.r -i 1.rawdata/${subpath} -w 2.cellanno/$s/ -o ${s} -s Mus > 2.cellanno/${s}/1.QualityControl.log
done


### 3.combined
combined="Final_20w_WTAPC"
mkdir $combined
cd $combined
ln -s ../../2.cellanno/W002/*_raw.rds W_20w_1_raw.rds
ln -s ../../2.cellanno/A002/*_raw.rds A_20w_1_raw.rds
ln -s ../../2.cellanno/W004/*_raw.rds W_20w_2_raw.rds
ln -s ../../2.cellanno/A004/*_raw.rds A_20w_2_raw.rds
ln -s ../../2.cellanno/A005/*_raw.rds W_20w_3_raw.rds
ln -s ../../2.cellanno/A006/*_raw.rds A_20w_3_raw.rds
ln -s ../../2.cellanno/W005/*_raw.rds W_20w_4_raw.rds
ln -s ../../2.cellanno/W006/*_raw.rds A_20w_4_raw.rds
cd ..

sample_names="W_20w_1,A_20w_1,W_20w_2,A_20w_2,W_20w_3,A_20w_3,W_20w_4,A_20w_4"
rds_files="W_20w_1_raw.rds,A_20w_1_raw.rds,W_20w_2_raw.rds,A_20w_2_raw.rds,W_20w_3_raw.rds,A_20w_3_raw.rds,W_20w_4_raw.rds,A_20w_4_raw.rds"
filter_fea="7500,7500,7500,7500,7500,7500,7500,7500"
filter_mt="10,10,10,10,10,10,10,10"
Rscript ${R_script}/2.IntegratedSamples.Seurat3_PCAselection.r -w 3.combination/$combined -l $sample_names -f $rds_files -u $filter_fea -m $filter_mt -o ${combined} > 3.combination/2.PCAselection.log
Rscript ${R_script}/CCA_3.Seurat3_tSNEorUMAP.r -w 3.combination/$combined -p 30 -f ${combined}_Origin_Integrated.rds -o ${combined} > 3.combination/3.CCA_p30.log
