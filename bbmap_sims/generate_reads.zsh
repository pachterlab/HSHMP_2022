# Mature reads for nuclear simulation
/home/kristjan/bbmap/randomreads.sh seed=42 ref=/home/kristjan/kallisto_bf_analysis/partial_transcriptomes/f1_starsolo_nopolyA.fa paired=f out=/home/kristjan/kallisto_bf_analysis/simulated_reads/nuclear/f1_hg38_mature_only_1000000_reads_20221126.fq q=36 adderrors=f minlen=150 maxlen=150 reads=1000000

# Nascent reads for nuclear simulation
/home/kristjan/bbmap/randomreads.sh seed=42 ref=/home/kristjan/kallisto_bf_analysis/partial_transcriptomes/nascent_starsolo_v2.fa paired=f out=/home/kristjan/kallisto_bf_analysis/simulated_reads/nuclear/f1_hg38_nascent_only_4000000_reads_20221126.fq q=36 adderrors=f minlen=150 maxlen=150 reads=4000000

# Fake 10xV3 format for reads
scripts/bbmap_to_10xv3.py --nascent-in /home/kristjan/kallisto_bf_analysis/simulated_reads/nuclear/f1_hg38_nascent_only_4000000_reads_20221126.fq --nascent-out /home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/nuclear/Nascent --mature-in /home/kristjan/kallisto_bf_analysis/simulated_reads/nuclear/f1_hg38_mature_only_1000000_reads_20221126.fq --mature-out /home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/nuclear/Mature --r1-template /home/kristjan/cellranger_comparison/samples/5k_pbmc_sims_v3_S1_L001_R1_001.fastq.gz
