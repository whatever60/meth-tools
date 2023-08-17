fa_dir=~/data/genome
input_fa=$fa_dir/GRCm39.primary_assembly.genome.fa
converted_fa_c2t=$fa_dir/GRCm39.primary_assembly.genome.c2t.fa
converted_fa_c2t_f=$fa_dir/GRCm39.primary_assembly.genome.c2t.f.fa
converted_fa_c2t_r=$fa_dir/GRCm39.primary_assembly.genome.c2t.r.fa
converted_fa_t2c=$fa_dir/GRCm39.primary_assembly.genome.t2c.fa
converted_fa_t2c_f=$fa_dir/GRCm39.primary_assembly.genome.t2c.f.fa
converted_fa_t2c_r=$fa_dir/GRCm39.primary_assembly.genome.t2c.r.fa

input_gff=~/data/genome/gencode.vM32.annotation.gff3
converted_gff_f=~/data/genome/gencode.vM32.annotation.f.gff3
converted_gff_r=~/data/genome/gencode.vM32.annotation.r.gff3

star_index_dir_prefix=~/s3_wang/data/index_files/STAR_2.7.10b/gencode_m32
unconverted_suffix=""
# for DNA methylation (BS-seq)
c2t_forward_suffix="_c2t_f"
c2t_reverse_suffix="_c2t_r"
# for RNA methylation (GLORI-seq)
t2c_forward_suffix="_t2c_f"
t2c_reverse_suffix="_t2c_r"

# convert genome with bwa-meth (make sure bwa-meth dependency is installed)
python3 ~/src/bwa-meth/bwameth.py c2t_fa $input_fa --config bs-seq --out_dir $fa_dir
python3 ~/src/bwa-meth/bwameth.py c2t_fa $input_fa --config glori --out_dir $fa_dir
# split the fasta into two, since bwa-meth will duplicate the genome by converting both strands
awk '/^>/{ if ($0 ~ /^>f/) filename="'$converted_fa_c2t_f'"; else filename="'$converted_fa_c2t_r'"; } { print > filename; }' $converted_fa_c2t
awk '/^>/{ if ($0 ~ /^>f/) filename="'$converted_fa_t2c_f'"; else filename="'$converted_fa_t2c_r'"; } { print > filename; }' $converted_fa_t2c

# convert genome annotation, i.e., changing chromosome names according to the new genome
awk 'BEGIN {FS="\t"; OFS="\t"}
     !/^#/ { $1 = "f" $1; print; next }
     /^##sequence-region/ { split($0, a, " "); a[2] = "f" a[2]; printf "%s %s %s %s\n", a[1], a[2], a[3], a[4]; next }
     1' $input_gff >$converted_gff_f
awk 'BEGIN {FS="\t"; OFS="\t"}
     !/^#/ { $1 = "r" $1; print; next }
     /^##sequence-region/ { split($0, a, " "); a[2] = "r" a[2]; printf "%s %s %s %s\n", a[1], a[2], a[3], a[4]; next }
     1' $input_gff >$converted_gff_r

# build STAR index for three genomes separately (unconverted, forward-strand converted, reverse-strand converted)
~/software/STAR/Linux_x86_64/STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir $star_index_dir_prefix$unconverted_suffix  \
    --genomeFastaFiles $input_fa \
    --sjdbGTFfile $input_gff \
    --sjdbOverhang 75
~/software/STAR/Linux_x86_64/STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir $star_index_dir_prefix$c2t_forward_suffix  \
    --genomeFastaFiles $converted_fa_c2t_f \
    --sjdbGTFfile $converted_gff_f \
    --sjdbOverhang 75
~/software/STAR/Linux_x86_64/STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir $star_index_dir_prefix$c2t_reverse_suffix  \
    --genomeFastaFiles $converted_fa_c2t_r \
    --sjdbGTFfile $converted_gff_r \
    --sjdbOverhang 75
~/software/STAR/Linux_x86_64/STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir $star_index_dir_prefix$t2c_forward_suffix  \
    --genomeFastaFiles $converted_fa_t2c_f \
    --sjdbGTFfile $converted_gff_f \
    --sjdbOverhang 75
~/software/STAR/Linux_x86_64/STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir $star_index_dir_prefix$t2c_reverse_suffix  \
    --genomeFastaFiles $converted_fa_t2c_r \
    --sjdbGTFfile $converted_gff_r \
    --sjdbOverhang 75

# perform alignment
fq_dir=~/s3_wang/20230814_glori/fastq_trim_galore
align_dir_prefix=~/s3_wang/20230814_glori/star_out
# S37, S38, S39 are samples without GLORI protocol, i.e., no conversion.

declare -A samples
samples=( ["S37"]="YQRNAseqModPilotS1" ["S38"]="YQRNAseqModPilotS2" ["S39"]="YQRNAseqModPilotS3" )
for sample in "${!samples[@]}"; do
    prefix=${samples[$sample]}
    
    ~/software/STAR/Linux_x86_64/STAR \
        --genomeDir $star_index_dir_prefix$unconverted_suffix \
        --readFilesIn $fq_dir/${prefix}_${sample}_R1_001_trimmed.fq $fq_dir/${prefix}_${sample}_R2_001_trimmed.fq \
        --runThreadN 16 \
        --outFileNamePrefix $align_dir_prefix$unconverted_suffix
done

# convert fastq using bwa-meth
fq_dir_conv=~/s3_wang/20230814_glori/fastq_trim_galore_t2c
declare -A samples
samples=( ["S40"]="YQRNAseqModPilotS4" ["S41"]="YQRNAseqModPilotS5" ["S42"]="YQRNAseqModPilotS6" )

for sample in "${!samples[@]}"; do
    prefix=${samples[$sample]}
    
    python3 bwameth.py c2t_fq \
        $fq_dir/${prefix}_${sample}_R1_001_trimmed.fq \
        $fq_dir/${prefix}_${sample}_R2_001_trimmed.fq \
        --config glori \
        --out_dir $fq_dir_conv
done

suffix="_001_trimmed.t2c.fq"
# align to converted forward genome
for sample in "${!samples[@]}"; do
    prefix=${samples[$sample]}
    
    ~/software/STAR/Linux_x86_64/STAR \
        --genomeDir $star_index_dir_prefix$t2c_forward_suffix \
        --readFilesIn $fq_dir_conv/${prefix}_${sample}_R1$suffix $fq_dir_conv/${prefix}_${sample}_R2$suffix \
        --runThreadN 16 \
        --outFileNamePrefix $align_dir_prefix$t2c_forward_suffix
done

# align to converted reverse genome
for sample in "${!samples[@]}"; do
    prefix=${samples[$sample]}
    
    ~/software/STAR/Linux_x86_64/STAR \
        --genomeDir $star_index_dir_prefix$t2c_reverse_suffix \
        --readFilesIn $fq_dir_conv/${prefix}_${sample}_R1$suffix $fq_dir_conv/${prefix}_${sample}_R2$suffix \
        --runThreadN 16 \
        --outFileNamePrefix $align_dir_prefix$t2c_reverse_suffix
done

