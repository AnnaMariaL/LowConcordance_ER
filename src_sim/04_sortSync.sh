#!/bin/bash
########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: sort sync files for same ordering as empirical data

########################
########################

#sync file: DsimFl_20200109-SNPSet01_hemizygousX-unzipped-ds80x.sync

#F10
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $11 "\t" $18 "\t" $25 "\t" $32 "\t" $39 "\t" $46 "\t" $53 "\t" $60 "\t" $67 "\t" $5 "\t" $12 "\t" $19 "\t" $26 "\t" $33 "\t" $40 "\t" $47 "\t" $54 "\t" $61 "\t" $68 }' DsimFl_20200109-SNPSet01_hemizygousX-unzipped-ds80x.sync > DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F10-ordered.sync

#F20
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $11 "\t" $18 "\t" $25 "\t" $32 "\t" $39 "\t" $46 "\t" $53 "\t" $60 "\t" $67 "\t" $6 "\t" $13 "\t" $20 "\t" $27 "\t" $34 "\t" $41 "\t" $48 "\t" $55 "\t" $62 "\t" $69 }' DsimFl_20200109-SNPSet01_hemizygousX-unzipped-ds80x.sync > DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F20-ordered.sync

#F30
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $11 "\t" $18 "\t" $25 "\t" $32 "\t" $39 "\t" $46 "\t" $53 "\t" $60 "\t" $67 "\t" $7 "\t" $14 "\t" $21 "\t" $28 "\t" $35 "\t" $42 "\t" $49 "\t" $56 "\t" $63 "\t" $70 }' DsimFl_20200109-SNPSet01_hemizygousX-unzipped-ds80x.sync > DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F30-ordered.sync

#F40
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $11 "\t" $18 "\t" $25 "\t" $32 "\t" $39 "\t" $46 "\t" $53 "\t" $60 "\t" $67 "\t" $8 "\t" $15 "\t" $22 "\t" $29 "\t" $36 "\t" $43 "\t" $50 "\t" $57 "\t" $64 "\t" $71 }' DsimFl_20200109-SNPSet01_hemizygousX-unzipped-ds80x.sync > DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F40-ordered.sync

#F50
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $11 "\t" $18 "\t" $25 "\t" $32 "\t" $39 "\t" $46 "\t" $53 "\t" $60 "\t" $67 "\t" $9 "\t" $16 "\t" $23 "\t" $30 "\t" $37 "\t" $44 "\t" $51 "\t" $58 "\t" $65 "\t" $72 }' DsimFl_20200109-SNPSet01_hemizygousX-unzipped-ds80x.sync > DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F50-ordered.sync

#F60
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $11 "\t" $18 "\t" $25 "\t" $32 "\t" $39 "\t" $46 "\t" $53 "\t" $60 "\t" $67 "\t" $10 "\t" $17 "\t" $24 "\t" $31 "\t" $38 "\t" $45 "\t" $52 "\t" $59 "\t" $66 "\t" $73 }' DsimFl_20200109-SNPSet01_hemizygousX-unzipped-ds80x.sync > DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F60-ordered.sync

#calculate the AF + Coverage
python compute_AF_cov_Ne.py --input ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F10-ordered.sync --output ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F10-ordered_AFrising_cov.txt --gen 10

python compute_AF_cov_Ne.py --input ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F20-ordered.sync --output ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F20-ordered_AFrising_cov.txt --gen 20

python compute_AF_cov_Ne.py --input ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F30-ordered.sync --output ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F30-ordered_AFrising_cov.txt --gen 30

python compute_AF_cov_Ne.py --input ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F40-ordered.sync --output ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F40-ordered_AFrising_cov.txt --gen 40

python compute_AF_cov_Ne.py --input ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F50-ordered.sync --output ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F50-ordered_AFrising_cov.txt --gen 50

python compute_AF_cov_Ne.py --input ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F60-ordered.sync --output ../data_sim/DsimFl_20200109-SNPSet01_hemizygousX-ds80x_F0F60-ordered_AFrising_cov.txt --gen 60
