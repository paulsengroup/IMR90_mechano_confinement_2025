#!/bin/bash
name=$(echo "$1" | cut -f 1 -d '.')
awk -v filename=$1 -v bin_size_arg=$2 '{ bin_size = bin_size_arg; for (counter_bin1 = $2; counter_bin1 < $3; counter_bin1 += bin_size) {
if (counter_bin1 + bin_size > $3){second_boundry = $3;} else {second_boundry = counter_bin1 + bin_size;}
for (counter_bin2 = $5; counter_bin2 < $6; counter_bin2 += bin_size){
if (counter_bin2 + bin_size > $6){fourth_boundry = $6;} else {fourth_boundry = counter_bin2 + bin_size;}
print $1  "\t" counter_bin1 "\t" second_boundry "\t" $4 "\t" counter_bin2 "\t" fourth_boundry;
}
}
}' $1 > tmp.bedpe

cut -f 1,2,3,4,5,6 tmp.bedpe | awk '{if(!(($1 == $4) && ($2 == $5) && ($3 == $6))) print $0}' > tmp1.bedpe
cut -f 1,2,3 tmp1.bedpe > tmp1_left.bedpe
cut -f 4,5,6 tmp1.bedpe  > tmp1_right.bedpe
cat tmp1_left.bedpe tmp1_right.bedpe | sort -u | bedtools sort -g hg38.chrom.sizes > tmp1_beads.bedpe
bedtools intersect -a $3 -b tmp1_beads.bedpe -v | cat - tmp1_beads.bedpe | bedtools sort -g hg38.chrom.sizes > tmp1_beads_complemented1.bedpe
bedtools makewindows -b  tmp1_beads_complemented1.bedpe -w $2 > tmp1_beads_complemented.bedpe
curl -s "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz" | gunzip -c | grep acen |  grep '^chr[XY0-9]\\+[[:space:]]' | bedtools pairtobed -a tmp1.bedpe -b stdin -type neither > tmp1_noncen.bedpe
curl -s "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz" | gunzip -c | grep '^chr[XY0-9]\\+[[:space:]]' |  bedtools pairtobed -a tmp1_noncen.bedpe -b stdin -type neither > tmp1_noncen_nongap.bedpe
python3 makeGtrack.py tmp1_noncen_nongap.bedpe tmp1_beads_complemented.bedpe > $name"_"$2_t.gtrack

awk 'BEGIN{print "###seqid\tstart\tend\tid\tradius\tperiphery\tedges"}1' $name"_"$2_t.gtrack > $name"_"$2.gtrack

rm tmp.bedpe tmp1.bedpe tmp1_left.bedpe tmp1_right.bedpe tmp1_beads.bedpe tmp1_beads_complemented1.bedpe tmp1_beads_complemented.bedpe tmp1_noncen.bedpe tmp1_noncen_nongap.bedpe $name"_"$2_t.gtrack
