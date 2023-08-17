. ./common.sh

set -e

# bsmooth won't let you specify abs path...
## INDEX
#mkdir bsmooth
#perl $BSMOOTH/bin/bswc_bowtie2_index.pl --name=bsmooth/$(basename $REF .fa) $REF

mkdir -p results/
prog=bsmooth

BAM=results/bsmooth/$prog-$name.bam

rm -f logs/$prog-$name.err logs/$prog-$name.out
mkdir -p results/$prog/$name/ 
echo "
perl $BSMOOTH/bin/bswc_bowtie2_align.pl \
    --metrics results/$prog/$name-metrics.txt \
    --out results/$prog/$name/ \
    --bam $BAM \
    -- bsmooth/$(basename $REF .fa) -- $REF \
    -- --very-sensitive -p 8 -- $FQ1 -- $FQ2

samtools view -h $BAM.crick.bam | python src/bsmooth-adjust-crick.py | samtools view -bS - > $BAM.crick.fix.bam
java -jar $PICARD/MergeSamFiles.jar I=$BAM.watson.bam I=$BAM.crick.fix.bam O=$BAM AS=false " | bsub -J $prog-$name -e logs/$prog-$name.err -o logs/$prog-$name.out -n 12 -R "span[hosts=1]"

mkdir -p results/trim/bsmooth/
BAM=results/trim/bsmooth/$prog-$name.bam
rm -f logs/trim-$prog-$name.err logs/trim-$prog-$name.out
mkdir -p results/trim/$prog/$name/ 
echo "
perl $BSMOOTH/bin/bswc_bowtie2_align.pl \
    --metrics results/trim/bsmooth/metrics.txt \
    --out results/trim/$prog/$name/ \
    --bam $BAM \
    -- bsmooth/$(basename $REF .fa) -- $REF \
    -- --very-sensitive -p 8 -- $TRIM_FQ1 -- $TRIM_FQ2

samtools view -h $BAM.crick.bam | python src/bsmooth-adjust-crick.py | samtools view -bS - > $BAM.crick.fix.bam
java -jar $PICARD/MergeSamFiles.jar I=$BAM.watson.bam I=$BAM.crick.fix.bam O=$BAM AS=false " | bsub -J trim-$prog-$name -e logs/trim-$prog-$name.err -o logs/trim-$prog-$name.out -n 12 -R "span[hosts=1]"
