. ./common.sh


REF=/data/Schwartz/brentp/mm10/ # bowtie

method=bis1

mkdir -p $OUT/$method/trim/

rm -f $OUT/$method/$name*.bam $OUT/trim/$method/$name*.bam

cmd="bismark -n 2 -l 24 --bam --temp_dir $TEMP --output_dir $OUT/$method/ --prefix $name $REF -1 $FQ1 -2 $FQ2;
samtools sort $OUT/$method/$name*.bam $OUT/$method-$name;
samtools index $OUT/$method-$name.bam"

rm -f logs/$method-$name.err logs/$method-$name.out
echo $cmd  | bsub -J $method-$name -e logs/$method-$name.err -o logs/$method-$name.out -n 2 -R "span[hosts=1]"

mkdir -p $TEMP/trim/

cmd="bismark -n 2 -l 24 --bam --temp_dir $TEMP/trim/ --output_dir $OUT/$method/trim/ --prefix $name $REF -1 $TRIM_FQ1 -2 $TRIM_FQ2;
samtools sort $OUT/$method/trim/$name*.bam $OUT/trim/$method-$name;
samtools index $OUT/trim/$method-$name.bam"

rm -f logs/trim-$method-$name.err logs/trim-$method-$name.out
echo $cmd | bsub -J trim-$method-$name -e logs/trim-$method-$name.err -o logs/trim-$method-$name.out -n 2 -R "span[hosts=1]"
