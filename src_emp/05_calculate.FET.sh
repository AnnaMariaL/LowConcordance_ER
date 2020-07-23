#!/bin/bash

########################
########################
#Author: Anna Maria Langmüller 
#Year: 2020
#Department: Vetmed University Vienna, Department of Population Genetics
#Project: Low concordance of short-term and long-term selection responses in experimental Drosophila populations
#Purpose: Calculate FET

########################
########################

perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep1.sync  --output F0F10SNP_Rep1.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep2.sync  --output F0F10SNP_Rep2.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep3.sync  --output F0F10SNP_Rep3.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep4.sync  --output F0F10SNP_Rep4.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep5.sync  --output F0F10SNP_Rep5.sync.fet --min-count 5 --min-coverage 2 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep6.sync  --output F0F10SNP_Rep6.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep7.sync  --output F0F10SNP_Rep7.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep8.sync  --output F0F10SNP_Rep8.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep9.sync  --output F0F10SNP_Rep9.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F10SNP_Rep10.sync  --output F0F10SNP_Rep10.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep1.sync  --output F0F20SNP_Rep1.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep2.sync  --output F0F20SNP_Rep2.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep3.sync  --output F0F20SNP_Rep3.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep4.sync  --output F0F20SNP_Rep4.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep5.sync  --output F0F20SNP_Rep5.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep6.sync  --output F0F20SNP_Rep6.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep7.sync  --output F0F20SNP_Rep7.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep8.sync  --output F0F20SNP_Rep8.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep9.sync  --output F0F20SNP_Rep9.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F20SNP_Rep10.sync  --output F0F20SNP_Rep10.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep1.sync  --output F0F30SNP_Rep1.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep2.sync  --output F0F30SNP_Rep2.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep3.sync  --output F0F30SNP_Rep3.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep4.sync  --output F0F30SNP_Rep4.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep5.sync  --output F0F30SNP_Rep5.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep6.sync  --output F0F30SNP_Rep6.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep7.sync  --output F0F30SNP_Rep7.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep8.sync  --output F0F30SNP_Rep8.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep9.sync  --output F0F30SNP_Rep9.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F30SNP_Rep10.sync  --output F0F30SNP_Rep10.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep1.sync  --output F0F40SNP_Rep1.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep2.sync  --output F0F40SNP_Rep2.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep3.sync  --output F0F40SNP_Rep3.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep4.sync  --output F0F40SNP_Rep4.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep5.sync  --output F0F40SNP_Rep5.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep6.sync  --output F0F40SNP_Rep6.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep7.sync  --output F0F40SNP_Rep7.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep8.sync  --output F0F40SNP_Rep8.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep9.sync  --output F0F40SNP_Rep9.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F40SNP_Rep10.sync  --output F0F40SNP_Rep10.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep1.sync  --output F0F50SNP_Rep1.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep2.sync  --output F0F50SNP_Rep2.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep3.sync  --output F0F50SNP_Rep3.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep4.sync  --output F0F50SNP_Rep4.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep5.sync  --output F0F50SNP_Rep5.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep6.sync  --output F0F50SNP_Rep6.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep7.sync  --output F0F50SNP_Rep7.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep8.sync  --output F0F50SNP_Rep8.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep9.sync  --output F0F50SNP_Rep9.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F50SNP_Rep10.sync  --output F0F50SNP_Rep10.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep1.sync  --output F0F60SNP_Rep1.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep2.sync  --output F0F60SNP_Rep2.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep3.sync  --output F0F60SNP_Rep3.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep4.sync  --output F0F60SNP_Rep4.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep5.sync  --output F0F60SNP_Rep5.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep6.sync  --output F0F60SNP_Rep6.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep7.sync  --output F0F60SNP_Rep7.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep8.sync  --output F0F60SNP_Rep8.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep9.sync  --output F0F60SNP_Rep9.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
perl ../src_emp/popoolation2/fisher-test.pl --input F0F60SNP_Rep10.sync  --output F0F60SNP_Rep10.sync.fet --min-count 5 --min-coverage 1 --max-coverage 1000
