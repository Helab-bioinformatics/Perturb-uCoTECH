##Perturb uCoTECH single-cell deconvolution code developed by Yaxi Liu, 2023.

#DNA
mkdir raw
mv *fq.gz raw/
cd raw/
for i in ./*_1.fq.gz; do base=$(basename $i "_1.fq.gz");  umi_tools whitelist --stdin $i --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCTCCACGC){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>ATAAGAGACAG[e<=1]).*'  --extract-method=regex --set-cell-number 1000 --ed-above-threshold=correct  --error-correct-threshold=2 --log2stderr --knee-method=distance > ${base}_whitelist.txt --plot-prefix ${base} ; done
cd ..
mkdir 01_extract
cd raw/
for i in ./*_whitelist.txt; do base=$(basename $i "_whitelist.txt");  umi_tools extract --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCTCCACGC){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>ATAAGAGACAG[e<=1]).*' --stdin ${base}_1.fq.gz  --stdout ../01_extract/${base}_1.extract.fq.gz  --read2-in  ${base}_2.fq.gz  --read2-out ../01_extract/${base}_2.extract.fq.gz  --error-correct-cell  --extract-method=regex  --whitelist ${base}_whitelist.txt; done

#RNA
for i in ./*_2.fq.gz; do base=$(basename $i "_2.fq.gz"); cutadapt -a AAAAAAAA  -A AAAAAAAA -q 20 -O 8  --trim-n  -m 6  --max-n 0.1 -o ${base}_1.trim.fq.gz -p ${base}_2.trim.fq.gz ${base}_1.fq.gz ${base}_2.fq.gz ; done
for i in ./*_1.fq.gz; do base=$(basename $i "_1.fq.gz");  umi_tools whitelist --stdin $i --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCT){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>TTTTTT){e<=1}.*'  --extract-method=regex --ed-above-threshold=correct --error-correct-threshold=2 --set-cell-number 1000 --knee-method=distance --log2stderr > ${base}_whitelist.txt --plot-prefix ${base}_expect_whitelist; done
cd ..
mkdir 01_extract
cd raw/
for i in ./*_whitelist.txt; do base=$(basename $i "_whitelist.txt");  umi_tools extract --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCT){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>TTTTTT){e<=1}.*' --stdin ${base}_1.fq.gz  --stdout ../01_extract/${base}_1.extract.fq.gz  --read2-in  ${base}_2.fq.gz  --read2-out ../01_extract/${base}_2.extract.fq.gz   --error-correct-cell --extract-method=regex  --whitelist ${base}_whitelist.txt; done
cd ../01_extract/
for i in ./*_2.extract.fq.gz; do base=$(basename $i "_2.extract.fq.gz"); cutadapt -a AAAAAAAA  -A AAAAAAAA -q 20 -O 8  --trim-n  -m 6  --max-n 0.1 -o ${base}_1.trim.fq.gz -p ${base}_2.trim.fq.gz ${base}_1.extract.fq.gz ${base}_2.extract.fq.gz ; done
for i in ./*_2.trim.fq.gz; do base=$(basename $i "_2.trim.fq.gz"); cutadapt -a  CTGTCTCTTATACAC  -q 20 -O 8  --trim-n  -m 20  --max-n 0.1 -o ${base}_2.extract.trim.fq.gz ${base}_2.trim.fq.gz; done

#sgRNA
for i in ./*_1.fq.gz; do base=$(basename $i "_1.fq.gz");  umi_tools whitelist --stdin $i --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCT){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>TTTTTT){e<=1}.*'  --extract-method=regex --ed-above-threshold=correct --error-correct-threshold=2 --set-cell-number 1000 --knee-method=distance --log2stderr > ${base}_whitelist.txt --plot-prefix ${base}_expect_whitelist; done
cd ..
mkdir 01_extract
cd raw/
for i in ./*_whitelist.txt; do base=$(basename $i "_whitelist.txt");  umi_tools extract --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCT){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>TTTTTT){e<=1}.*' --stdin ${base}_1.fq.gz  --stdout ../01_extract/${base}_1.extract.fq.gz  --read2-in  ${base}_2.fq.gz  --read2-out ../01_extract/${base}_2.extract.fq.gz   --error-correct-cell --extract-method=regex  --whitelist ${base}_whitelist.txt; done


