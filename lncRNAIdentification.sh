#!/bin/bash
export PATH=/home/dklly/miniconda3/envs/bioinformatics/bin:$PATH
workDir="/home/dklly/mammaryGland"
index="~/genomeIndex/hisat2Index/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel"
genome="/home/dklly/genomes/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa"
refgtf="/home/dklly/genomes/Bos_taurus.ARS-UCD1.2.96.gtf"
knonwSplitSite="~/genomes/Bos_taurus.ARS-UCD1.2.96.ss"
gtfTools="/home/dklly/inHoseGTFtools/inHouseGTFtools"

## mapping to genome
mkdir sam
cd fastq

ls *.fastq|while read id
do hisat2 -p 4 --dta -x $index \
--known-splicesite-infile $knonwSplitSite $id \
-S ../sam/${id%.*}.sam &>> ../sam/alignmentSummary.txt
done

## covert sam to bam
cd ../sam
ls *.sam|while read id
do 
samtools view -bS $id > ${id%.*}.bam
done

## sort bam
ls *.bam|while read id
do 
samtools sort -@ 4 $id ${id%.*}.sorted
rm ${id%.*}.bam
done

## reconstruct transcripts using Stringtie
mkdir ../gtfs
ls *.sorted.bam|while read id
do stringtie $id -p 4 -G $refgtf  -o ../gtfs/${id}.gtf -l ${id%.*}
done

## merge transcript
cd ../gtfs
ls *.gtf > merge.list
stringtie --merge -G $refgtf merge.list -o merged.gtf

## quantity
cd ../sam

ls *sorted.bam|while read id
do stringtie -e -B -p 4 -G ../gtfs/merged.gtf -o ../ballgown/${id%%.*}/${id%%.*}.bg.gtf $id
done

## obtain FPKM matrix
cd ..
mkdir quantity

prepDE.py ballgown
mv *.csv quantity


## obtain FPKM
cd quantity

ls $workDir/ballgown/|while read id
do 
less $workDir/ballgown/$id/$id.bg.gtf| \
awk '$3=="transcript"'|grep -Eo 'transcript_id \"\w+|transcript_id \"\w+\.\w+\.\w+'|cut -d\" -f 2 > $id.t
done
ls $workDir/ballgown/|while read id 
do
less $workDir/ballgown/$id/$id.bg.gtf|awk '$3=="transcript"'|grep -Eo 'FPKM \"\w+\.\w+'|awk -F\" '{print$2}' > $id.value
done
ls $workDir/ballgown/|while read id; do paste $id.t $id.value > $id.FPKM;done
ls *.FPKM|while read id;do less $id|sort -k 1 > ${id%.*}.sorted.FPKM;done
echo "transcript" > FPKM
less `ls *.sorted.FPKM|shuf -n1`|cut -f 1 >> FPKM
ls *.sorted.FPKM|while read id;do echo ${id%%.*} > ${id%%.*}.tmp;less $id|cut -f 2 >> ${id%%.*}.tmp;done
paste FPKM *.tmp > all.FPKM
ls *.value *.t *.FPKM *.tmp FPKM|grep -v all.FPKM |xargs rm

## idenfication
### compare to reference transcripts, parameter M refer to remove single exon transcripts
mkdir ../Identification
cd ../Identification
gffcompare -R -M -r $refgtf -s $genome ../gtfs/merged.gtf


## class code filteration

$gtfTools STscreen -classCode u -i gffcmp.annotated.gtf > u.gtf
## read sequence
gffread -w u.fa -g $genome u.gtf
samtools faidx u.fa
less u.fa.fai |awk '$2>200{print$1}' > uL200.MSTRG

## length filteration
less u.fa.fai |awk '$2>200{print$1}' > uL200.MSTRG

## FPKM filteraion
python ../scripts/extractFPKMtable.py ../quautity/all.FPKM uL200.MSTRG |python ../scripts/filterFPKM.py 1 > uL200F0.5.MSTRG

## coding ability prediction
....





















