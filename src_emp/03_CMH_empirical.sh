#!/bin/bash

########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: CMH test, empirical data 

########################
########################
cd ../data_emp
perl ../src_emp/popoolation2/cmh-test.pl --input F0F10SNP.sync --output F0F10SNP.sync.cmh --min-count 10 --min-coverage 1 --max-coverage 1000 --population 1-11,2-12,3-13,4-14,5-15,6-16,7-17,8-18,9-19,10-20 --remove-temp

perl ../src_emp/popoolation2/cmh-test.pl --input F0F20SNP.sync --output F0F20SNP.sync.cmh --min-count 10 --min-coverage 1 --max-coverage 1000 --population 1-11,2-12,3-13,4-14,5-15,6-16,7-17,8-18,9-19,10-20 --remove-temp

perl ../src_emp/popoolation2/cmh-test.pl --input F0F30SNP.sync --output F0F30SNP.sync.cmh --min-count 10 --min-coverage 1 --max-coverage 1000 --population 1-11,2-12,3-13,4-14,5-15,6-16,7-17,8-18,9-19,10-20 --remove-temp

perl ../src_emp/popoolation2/cmh-test.pl --input F0F40SNP.sync --output F0F40SNP.sync.cmh --min-count 10 --min-coverage 1 --max-coverage 1000 --population 1-11,2-12,3-13,4-14,5-15,6-16,7-17,8-18,9-19,10-20 --remove-temp

perl ../src_emp/popoolation2/cmh-test.pl --input F0F50SNP.sync --output F0F50SNP.sync.cmh --min-count 10 --min-coverage 1 --max-coverage 1000 --population 1-11,2-12,3-13,4-14,5-15,6-16,7-17,8-18,9-19,10-20 --remove-temp

perl ../src_emp/popoolation2/cmh-test.pl --input F0F60SNP.sync --output F0F60SNP.sync.cmh --min-count 10 --min-coverage 1 --max-coverage 1000 --population 1-11,2-12,3-13,4-14,5-15,6-16,7-17,8-18,9-19,10-20 --remove-temp
