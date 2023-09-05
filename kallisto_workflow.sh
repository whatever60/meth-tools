# do this for all of S1_S37, S2_S38, S3_S39 (without GLORI protocol)
for i in {1..3}; do
    kallisto quant -i /home/ubuntu/s3_wang/data/index_files/kallisto_0.50.0/gencode.vM32.transcripts.idx \
        -o ~/s3_wang/20230814_glori/kallisto_out/s${i} \
        /home/ubuntu/s3_wang/20230814_glori/fastq_trim_galore/YQRNAseqModPilotS${i}_S$((i+36))_R1_001.fq \
        /home/ubuntu/s3_wang/20230814_glori/fastq_trim_galore/YQRNAseqModPilotS${i}_S$((i+36))_R2_001.fq \
        -t 16
done

# do this for all samples (with or without GLORI protocol)
for i in {1..6}; do
    kallisto quant -i /home/ubuntu/s3_wang/data/index_files/kallisto_0.50.0/gencode.vM32.transcripts.t2c.combine.idx \
        -o ~/s3_wang/20230814_glori/kallisto_out/s${i}_t2c \
        /home/ubuntu/s3_wang/20230814_glori/fastq_trim_galore_t2c/YQRNAseqModPilotS${i}_S$((i+36))_R1_001.t2c.fq \
        /home/ubuntu/s3_wang/20230814_glori/fastq_trim_galore_t2c/YQRNAseqModPilotS${i}_S$((i+36))_R2_001.t2c.fq \
        -t 16
done