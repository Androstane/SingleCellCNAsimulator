#!/usr/bin/env bash

if [ $# -lt 2 ]; then
    echo "location of simulated results from first step required\n";
    exit
fi

path=$1
name=$2



python main.par.py -S wgsim-master/ -p 5 -t ref/hg19.fa -k 1 -r ${path}
cd ~/simulator_test/SingleCellCNABenchmark_modified/${path}
mkdir fastq
mkdir bam 
mkdir bed 
mv *fq fastq/
mv *fa fastq/
mv *.fai fastq/
cp ~/simulator_test/SingleCellCNABenchmark_modified/merge_allele.py ~/simulator_test/SingleCellCNABenchmark_modified/${path}fastq/merge_allele.py
cd fastq/
python merge_allele.py 
bash merge_allele.sh
rm *allele0_1.fq
rm *allele0_2.fq
rm *allele1_1.fq
rm *allele1_2.fq
cp ~/simulator_test/SingleCellCNABenchmark_modified/generate_bam.py ~/simulator_test/SingleCellCNABenchmark_modified/${path}fastq/generate_bam.py
python generate_bam.py

sed -i "s&path&${path}/fastq&g" ~/simulator_test/SingleCellCNABenchmark_modified/${path}fastq/gen_bam.sh
mv ~/simulator_test/SingleCellCNABenchmark_modified/${path}fastq/gen_bam.sh ~/simulator_test/SingleCellCNABenchmark_modified/gen_bam_${name}.sh
cd ~/simulator_test/SingleCellCNABenchmark_modified/
bash gen_bam_${name}.sh
cd ~/simulator_test/SingleCellCNABenchmark_modified/${path}fastq

for file in ~/simulator_test/SingleCellCNABenchmark_modified/${path}fastq/*.bam 
do 
	bedtools bamtobed -i $file > ${file}.bed 
done

mv ~/simulator_test/SingleCellCNABenchmark_modified/${path}fastq/*.bed ~/simulator_test/SingleCellCNABenchmark_modified/${path}bed/
cd ~/simulator_test/SingleCellCNABenchmark_modified/
ginkgo-master/cli/ginkgo.sh --input ~/simulator_test/SingleCellCNABenchmark_modified/${path}bed/ --genome hg19 --binning variable_175000_48_bwa




