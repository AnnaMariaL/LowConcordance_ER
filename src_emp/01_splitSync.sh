########################
########################

#Anna Maria Langmüller 
#2019, Vetmed University Vienna
#Analysis to: "Early selection signatures do not reflect long-term adaptation responses in experimental Drosophila populations"
#original code: Neda Barghi "Genetic redundancy fuels polygenic adaptation in Drosophila", adapted by Anna Maria Langmüller 

########################
########################

#separate F0 and F60 samples from F0-F60SNP_CMH_FET_blockID.sync file (Dryad Digital Repository: https://doi.org/10.5061/dryad.rr137kn) 
cd ../data_emp

#F10
cut -f 1-23 Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > F0F10SNP.sync
python ../src_emp/compute_AF_cov_Ne.py --input F0F10SNP.sync --output F0F10SNP_AFrising_cov.txt --gen 10
#F20
cut -f 1-13,24-33 Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > F0F20SNP.sync
python ../src_emp/compute_AF_cov_Ne.py --input F0F20SNP.sync --output F0F20SNP_AFrising_cov.txt --gen 20
#F30
cut -f 1-13,34-43 Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > F0F30SNP.sync
python ../src_emp/compute_AF_cov_Ne.py --input F0F30SNP.sync --output F0F30SNP_AFrising_cov.txt --gen 30
#F40
cut -f 1-13,44-53 Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > F0F40SNP.sync
python ../src_emp/compute_AF_cov_Ne.py --input F0F40SNP.sync --output F0F40SNP_AFrising_cov.txt --gen 40
#F50
cut -f 1-13,54-63 Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > F0F50SNP.sync
python ../src_emp/compute_AF_cov_Ne.py --input F0F50SNP.sync --output F0F50SNP_AFrising_cov.txt --gen 50
#F60
cut -f 1-13,64-73 Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > F0F60SNP.sync
python ../src_emp/compute_AF_cov_Ne.py --input F0F60SNP.sync --output F0F60SNP_AFrising_cov.txt --gen 60

#calculate joint file 
python 0../src_emp/compute_AF_cov_Ne_all.py --input cut -f 1-73 Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.AFrising.txt