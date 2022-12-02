<pre>cd /home/dsullivan/benchmarking/kristjan/</pre>

## Download mouse references

<pre>wget --continue https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget --continue https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.108.gtf.gz</pre>


<pre>kallisto="/home/dsullivan/benchmarking/kristjan/kallisto-bf/build/src/kallisto"
prev_kallisto="/home/dsullivan/kallisto/build/src/kallisto"
bustools="/home/dsullivan/bustools/build/src/bustools"

salmon="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/exe/salmon_1.9.0"
af="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/exe/alevin-fry_0.8.0"

star="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript//exe/STAR_2.7.9a"

cellranger="/home/kristjan/cellranger/cellranger-7.0.1/cellranger"
</pre>

Paths to indices:

<pre>main_path="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript"
genome_name="human_CR_3.0.0"
genome_file="$main_path/genomes/$genome_name/genome.fa"
mouse_genome_file="Mus_musculus.GRCm39.dna.primary_assembly.fa"
mouse_gtf_file="Mus_musculus.GRCm39.108.gtf"
mouse_genome_name="mouse_CR"
mouse_nascent_genome_file="/home/kristjan/kallisto_bf_analysis/partial_transcriptomes/mus_musculus_nascent_v2.fa"
transcripts_file="$main_path/genomes/$genome_name/transcripts.fa"
gtf_file="$main_path/genomes/$genome_name/annotations.gtf"
n_threads="20"

salmon_index_standard="$main_path/genomes/index/salmon_1.9.0/$genome_name/standard/index"
t2gFile_salmon_standard="$main_path/genomes/$genome_name/transcript_to_gene.2col.txt"
salmon_index_splici="$main_path/genomes/index/salmon_1.9.0/$genome_name/splici/i150"
t2gFile_salmon_splici="$main_path/genomes/index/salmon_1.9.0/$genome_name/splici/salmon_splici_150/splici_fl145_t2g_3col.tsv"
kallisto_index="$main_path/genomes/index/kallisto_0.49.0/$genome_name/standard_1/index.idx"
kallisto_index_offlist="$main_path/genomes/index/kallisto_0.49.0/$genome_name/standard_offlist_1/index.idx"
star_index="$main_path/genomes/index/STAR_2.7.9a/human_CR_3.0.0/fullSA/"</pre>

## Link data

<pre>mkdir -p data
ln -s /home/ggorin/datasets/ ./data/</pre>


## Generate Mouse Indexes

<pre>$cellranger mkref --genome=$mouse_genome_name --fasta=$mouse_genome_file --genes=$mouse_gtf_file --nthreads=$n_threads</pre>

<pre>out_dir="kallisto_index_mouse"
mkdir -p $out_dir
kb ref -i $out_dir/index.idx --kallisto $kallisto --workflow standard --overwrite -f1 $out_dir/f1 -g $out_dir/g $mouse_genome_file $mouse_gtf_file > $out_dir/log.txt 2>&1
$kallisto index -i $out_dir/index_standard.idx $out_dir/f1 # TODO: DELETE THIS ONCE WE FIGURE OUT WHY TF KB ISN'T WORKING!
$kallisto index -t $n_threads -b $mouse_genome_file -i $out_dir/index_offlist.idx $out_dir/f1
$kallisto index -t $n_threads -b $out_dir/f1 -i $out_dir/index_nucleus.idx $mouse_nascent_genome_file</pre>

<pre>out_dir="star_index_mouse"
mkdir -p $out_dir/fullSA
mkdir -p $out_dir/sparseSA3
$star --runMode genomeGenerate --runThreadN $n_threads --genomeDir $out_dir/fullSA --genomeFastaFiles $mouse_genome_file --sjdbGTFfile $mouse_gtf_file > $out_dir/fullSA/log.txt 2>&1
$star --runMode genomeGenerate --runThreadN $n_threads --genomeDir $out_dir/sparseSA3 --genomeSAsparseD 3 --genomeFastaFiles $mouse_genome_file --sjdbGTFfile $mouse_gtf_file > $out_dir/sparseSA3/log.txt 2>&1
</pre>

<pre>out_dir="salmon_index_mouse"
mkdir -p $out_dir/standard
mkdir -p $out_dir/standard_sparse
mkdir -p $out_dir/splici
mkdir -p $out_dir/splici_sparse
pyroe make-splici "$mouse_genome_file" "$mouse_gtf_file" 90 $out_dir/splici/salmon_splici_90 --flank-trim-length 5 --filename-prefix splici
$salmon index --keepDuplicates -t kallisto_index_mouse/f1 -i $out_dir/standard/index -p $n_threads
$salmon index --keepDuplicates -t kallisto_index_mouse/f1 -i $out_dir/standard_sparse/index -p $n_threads --sparse
$salmon index -t $out_dir/splici/salmon_splici_90/splici_fl85.fa -i $out_dir/splici/index -p $n_threads
$salmon index -t $out_dir/splici/salmon_splici_90/splici_fl85.fa -i $out_dir/splici_sparse/index -p $n_threads
</pre>

## CellRanger Run

<pre>/usr/bin/time -v $cellranger count --localcores $n_threads --fastqs data/datasets/brain_10x_5k_fastqs/ --id sc_mouse_brain_cellranger7 --transcriptome $mouse_genome_name  1> sc_mouse_brain_cellranger7_stdout.txt 2> sc_mouse_brain_cellranger7_stderr.txt
/usr/bin/time -v $cellranger count --localcores $n_threads --fastqs data/datasets/brain_nuc_10x_5k_fastqs/ --id sn_mouse_brain_cellranger7 --transcriptome $mouse_genome_name  1> sn_mouse_brain_cellranger7_stdout.txt 2> sn_mouse_brain_cellranger7_stderr.txt</pre>

## Kallisto Run

<pre>paste <(cat kallisto_index_mouse/g|cut -f2|cut -d. -f1) <(cat kallisto_index_mouse/g|cut -f2|cut -d. -f1) > kallisto_index_mouse/g_</pre>

<pre>data_files="data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R1_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R2_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R1_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R2_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R1_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R2_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R1_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R2_001.fastq.gz"

data_files_nuc="data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L001_R1_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L001_R2_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L002_R1_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L002_R2_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L003_R1_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L003_R2_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L004_R1_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L004_R2_001.fastq.gz"</pre>

<pre>/usr/bin/time -v kb count --overwrite --kallisto $kallisto --bustools $bustools -i kallisto_index_mouse/index_standard.idx -g kallisto_index_mouse/g -t $n_threads -x 10XV3 -o sc_mouse_brain_kallisto_standard/ $data_files  1> sc_mouse_brain_kallisto_standard_stdout.txt 2> sc_mouse_brain_kallisto_standard_stderr.txt
/usr/bin/time -v kb count --overwrite --kallisto $kallisto --bustools $bustools -i kallisto_index_mouse/index_standard.idx -g kallisto_index_mouse/g -t $n_threads -x 10XV3 -o sn_mouse_brain_kallisto_standard/ $data_files_nuc  1> sn_mouse_brain_kallisto_standard_stdout.txt 2> sn_mouse_brain_kallisto_standard_stderr.txt</pre>

<pre>/usr/bin/time -v kb count --overwrite --kallisto $kallisto --bustools $bustools -i kallisto_index_mouse/index_offlist.idx -g kallisto_index_mouse/g -t $n_threads -x 10XV3 -o sc_mouse_brain_kallisto_offlist/ $data_files 1> sc_mouse_brain_kallisto_offlist_stdout.txt 2> sc_mouse_brain_kallisto_offlist_stderr.txt
/usr/bin/time -v kb count --overwrite --kallisto $kallisto --bustools $bustools -i kallisto_index_mouse/index_nucleus.idx -g kallisto_index_mouse/g_ -t $n_threads -x 10XV3 -o sn_mouse_brain_kallisto_offlist/ $data_files_nuc 1> sn_mouse_brain_kallisto_offlist_stdout.txt 2> sn_mouse_brain_kallisto_offlist_stderr.txt</pre>

## Salmon run

<pre>data_files="-1 data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R1_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R1_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R1_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R1_001.fastq.gz -2 data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R2_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R2_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R2_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R2_001.fastq.gz"</pre>

<pre>
data_files_nuc="-1 data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L001_R1_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L002_R1_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L003_R1_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L004_R1_001.fastq.gz -2 data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L001_R2_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L002_R2_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L003_R2_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L004_R2_001.fastq.gz"
</pre>



<pre>t2g="salmon_index_mouse/splici/salmon_splici_90/splici_fl85_t2g_3col.tsv"

out="sc_mouse_brain_salmon_splici_cr_like"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt
out="sn_mouse_brain_salmon_splici_cr_like"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files_nuc" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt

out="sc_mouse_brain_salmon_splici_sketch_cr_like"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --sketch --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt
out="sn_mouse_brain_salmon_splici_sketch_cr_like"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --sketch --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files_nuc" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt

out="sc_mouse_brain_salmon_splici_cr_like_em"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like-em -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt
out="sn_mouse_brain_salmon_splici_cr_like_em"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like-em -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files_nuc" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt

out="sc_mouse_brain_salmon_sparse_splici_cr_like"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici_sparse/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt
out="sn_mouse_brain_salmon_sparse_splici_cr_like"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici_sparse/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files_nuc" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt

out="sc_mouse_brain_salmon_sparse_splici_sketch_cr_like"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --sketch --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici_sparse/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt
out="sn_mouse_brain_salmon_sparse_splici_sketch_cr_like"
mkdir -p $out
/usr/bin/time -v sh -c '$0 alevin --sketch --chromiumV3 -p $1 -l IU -i salmon_index_mouse/splici_sparse/index --tgMap $2 --rad -o $3/ $4 && $5 generate-permit-list -d fw -u sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt -i $3/ -o $3/ && $5 collate -t $1 -i $3/ -r $3/ && $5 quant --resolution cr-like -t $1 -i $3/ -o $3/ --use-mtx --tg-map $2' $salmon $n_threads $t2g $out "$data_files_nuc" $af  1> "$out"_stdout.txt 2> "$out"_stderr.txt

</pre>

## STAR run

<pre>data_files="data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R2_001.fastq.gz,data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R2_001.fastq.gz,data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R2_001.fastq.gz,data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R2_001.fastq.gz data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R1_001.fastq.gz,data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R1_001.fastq.gz,data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R1_001.fastq.gz,data/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R1_001.fastq.gz"</pre>

<pre>
data_files_nuc="data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L001_R2_001.fastq.gz,data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L002_R2_001.fastq.gz,data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L003_R2_001.fastq.gz,data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L004_R2_001.fastq.gz data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L001_R1_001.fastq.gz,data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L002_R1_001.fastq.gz,data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L003_R1_001.fastq.gz,data/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L004_R1_001.fastq.gz"
</pre>


<pre>/usr/bin/time -v $star --genomeDir star_index_mouse/fullSA/ --runThreadN $n_threads --readFilesCommand zcat --readFilesIn $data_files --soloCBwhitelist sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt --soloUMIlen 12 --limitIObufferSize 50000000 50000000 --soloType CB_UMI_Simple --outSAMtype None --outFileNamePrefix sc_mouse_brain_star_full/ 1> sc_mouse_brain_star_full_stdout.txt 2> sc_mouse_brain_star_full_stderr.txt

/usr/bin/time -v $star --genomeDir star_index_mouse/sparseSA3/ --runThreadN $n_threads --readFilesCommand zcat --readFilesIn $data_files --soloCBwhitelist sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt --soloUMIlen 12 --limitIObufferSize 50000000 50000000 --soloType CB_UMI_Simple --outSAMtype None --outFileNamePrefix sc_mouse_brain_star_sparse/ 1> sc_mouse_brain_star_sparse_stdout.txt 2> sc_mouse_brain_star_sparse_stderr.txt</pre>

<pre>/usr/bin/time -v $star --genomeDir star_index_mouse/fullSA/ --soloFeatures GeneFull --runThreadN $n_threads --readFilesCommand zcat --readFilesIn $data_files_nuc --soloCBwhitelist sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt --soloUMIlen 12 --limitIObufferSize 50000000 50000000 --soloType CB_UMI_Simple --outSAMtype None --outFileNamePrefix sn_mouse_brain_star_full/ 1> sn_mouse_brain_star_full_stdout.txt 2> sn_mouse_brain_star_full_stderr.txt

/usr/bin/time -v $star --genomeDir star_index_mouse/sparseSA3/ --soloFeatures GeneFull --runThreadN $n_threads --readFilesCommand zcat --readFilesIn $data_files_nuc --soloCBwhitelist sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt --soloUMIlen 12 --limitIObufferSize 50000000 50000000 --soloType CB_UMI_Simple --outSAMtype None --outFileNamePrefix sn_mouse_brain_star_sparse/ 1> sn_mouse_brain_star_sparse_stdout.txt 2> sn_mouse_brain_star_sparse_stderr.txt</pre>

<pre>/usr/bin/time -v $star --genomeDir star_index_mouse/fullSA/ --soloFeatures Gene Velocyto --runThreadN $n_threads --readFilesCommand zcat --readFilesIn $data_files_nuc --soloCBwhitelist sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt --soloUMIlen 12 --limitIObufferSize 50000000 50000000 --soloType CB_UMI_Simple --outSAMtype None --outFileNamePrefix sn_mouse_brain_star_velocyto_full/ 1> sn_mouse_brain_star_velocyto_full_stdout.txt 2> sn_mouse_brain_star_velocyto_full_stderr.txt

/usr/bin/time -v $star --genomeDir star_index_mouse/sparseSA3/ --soloFeatures Gene Velocyto --runThreadN $n_threads --readFilesCommand zcat --readFilesIn $data_files_nuc --soloCBwhitelist sc_mouse_brain_kallisto_offlist/10x_version3_whitelist.txt --soloUMIlen 12 --limitIObufferSize 50000000 50000000 --soloType CB_UMI_Simple --outSAMtype None --outFileNamePrefix sn_mouse_brain_star_velocyto_sparse/ 1> sn_mouse_brain_star_velocyto_sparse_stdout.txt 2> sn_mouse_brain_star_velocyto_sparse_stderr.txt</pre>



## Extract runtimes and memory

### kallisto

<pre>
echo "program_run"$'\t'"map"$'\t'"other"$'\t'"memory" > kallisto_performance.txt
for f in *kallisto*stderr.txt; do
fname=$(echo $f|sed s/_stderr.txt//i)
tmap=$(cat $f|grep "Using index\|Sorting BUS"|head -2|cut -c 13-20|awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }'|awk 'NR > 1 { print $0 - first } { first = $0 }')
ttotal=$(cat $f|grep Elapsed|cut -d' ' -f8|cut -d. -f1|awk -F: '{ print ($1 * 60) + $2 }')
tother=$(($ttotal-$tmap))
memkb=$(cat $f|grep "Maximum resident set size"|cut -d':' -f2|tr -d ' ')
echo "$fname"$'\t'"$tmap"$'\t'"$tother"$'\t'"$memkb" >> kallisto_performance.txt
done
</pre>


<pre>
echo "program_run"$'\t'"map"$'\t'"other"$'\t'"memory" > salmon_performance.txt
for f in *salmon*stderr.txt; do
fname=$(echo $f|sed s/_stderr.txt//i)
tmap=$(cat $f|grep "rad flag\|finished sc_align"|head -2|cut -d'[' -f2|cut -c 12-19|awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }'|awk 'NR > 1 { print $0 - first } { first = $0 }')
ttotal=$(cat $f|grep Elapsed|cut -d' ' -f8|cut -d. -f1|awk -F: '{ print ($1 * 60) + $2 }')
tother=$(($ttotal-$tmap))
memkb=$(cat $f|grep "Maximum resident set size"|cut -d':' -f2|tr -d ' ')
echo "$fname"$'\t'"$tmap"$'\t'"$tother"$'\t'"$memkb" >> salmon_performance.txt
done
</pre>
