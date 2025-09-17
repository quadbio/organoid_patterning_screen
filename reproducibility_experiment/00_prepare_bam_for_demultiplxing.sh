## merge all libraries
samtools view -h s2/process/barcode_headAligned_anno.bam | sed -E 's/(CB:Z:[0-9]{2}_[0-9]{2}_[0-9]{2})/\1__s1/g' | samtools view -hbS > exp_fatima_e1_comb/for_vireo/sublibrary1.bam
samtools view -h s20/process/barcode_headAligned_anno.bam | sed -E 's/(CB:Z:[0-9]{2}_[0-9]{2}_[0-9]{2})/\1__s2/g' | samtools view -hbS > exp_fatima_e1_comb/for_vireo/sublibrary2_1.bam
samtools view -h s17/process/barcode_headAligned_anno.bam | sed -E 's/(CB:Z:[0-9]{2}_[0-9]{2}_[0-9]{2})/\1__s3/g' | samtools view -hbS > exp_fatima_e1_comb/for_vireo/sublibrary3.bam
samtools view -h s19/process/barcode_headAligned_anno.bam | sed -E 's/(CB:Z:[0-9]{2}_[0-9]{2}_[0-9]{2})/\1__s4/g' | samtools view -hbS > exp_fatima_e1_comb/for_vireo/sublibrary4.bam


samtools merge -@ 20 exp_fatima_e1_comb/for_vireo/combined_bam_files.bam exp_fatima_e1_comb/for_vireo/*.bam


#subset for pools wibj/wtc and H9/H1
conda activate spipe_new
samtools view -D CB:exp_fatima_e1_comb/for_vireo/barcode_WIBJ_WTC.lis exp_fatima_e1_comb/for_vireo/combined_bam_files.bam -h -b -o exp_fatima_e1_comb/for_vireo/WTC_WIBJ.bam
samtools view -D CB:exp_fatima_e1_comb/for_vireo/barcode_H9_H1.lis exp_fatima_e1_comb/for_vireo/combined_bam_files.bam -h -b -o exp_fatima_e1_comb/for_vireo/H9_H1.bam

# sort and index all of them
cd exp_fatima_e1_comb/for_vireo/

samtools sort -@ 10 -o sorted_WTC_WIBJ.bam WTC_WIBJ.bam
samtools sort -@ 10 -o sorted_H9_H1.bam H9_H1.bam


samtools index sorted_WTC_WIBJ.bam
samtools index sorted_H9_H1.bam


# call SNPs
cd /links/groups/treutlein/DATA/sequencing/20240906_P2911_NADEZHDA_PARSE/processed/exp_fatima_e1_comb
conda activate demuxlet

cellsnp-lite -s for_vireo/sorted_WTC_WIBJ.bam  -b for_vireo/barcode_WIBJ_WTC.lis -O cellSNP_WIBJ_WTC -R /links/groups/treutlein/USERS/nazbukina/General/variants/WIBJ2_WTC/WIBJ2_WTC_hg38_parse.vcf -p 10  --genotype --UMItag pN --gzip --minMAF 0.1 

cellsnp-lite -s for_vireo/sorted_H9_H1.bam  -b for_vireo/barcode_H9_H1.lis -O cellSNP_H9_H1 -R /links/groups/treutlein/USERS/nazbukina/General/variants/H9_H1/H9_H1_hg38_diversed_parse_true.vcf -p 10  --genotype --UMItag pN --gzip --minMAF 0.1 

#demultiplex
vireo -c cellSNP_WIBJ_WTC -d /links/groups/treutlein/USERS/nazbukina/General/variants/WIBJ2_WTC/WIBJ2_WTC_hg38_parse.vcf -o vireo_WIBJ_nodoub -t GT --noDoublet 
vireo -c cellSNP_H9_H1 -d /links/groups/treutlein/USERS/nazbukina/General/variants/H9_H1/H9_H1_hg38_diversed_parse_true.vcf -o vireo_H9_nodoub -t GT --noDoublet 








