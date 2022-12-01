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
$kallisto index -t $n_threads -b $mouse_genome_file -i $out_dir/index_offlist.idx $out_dir/f1</pre>

## CellRanger Run

<pre>/usr/bin/time -v $cellranger count --localcores $n_threads --fastqs data/datasets/brain_10x_5k_fastqs/ --id sc_mouse_brain_cellranger7 --transcriptome $mouse_genome_name  1> sc_mouse_brain_cellranger7_stdout.txt 2> sc_mouse_brain_cellranger7_stderr.txt</pre>

## Kallisto Run

<pre>/usr/bin/time -v kb count --kallisto $kallisto --bustools $bustools -i $out_dir/index_standard.idx -t $n_threads -x 10XV3 --fastqs data/datasets/brain_10x_5k_fastqs/ --id sc_mouse_brain_cellranger7 --transcriptome $mouse_genome_name  1> sc_mouse_brain_kallisto_standard_stdout.txt 2> sc_mouse_brain_kallisto_standard_stderr.txt</pre>

<pre>/usr/bin/time -v kb count --kallisto $kallisto --bustools $bustools -i $out_dir/index_offlist.idx -t $n_threads -x 10XV3 --fastqs data/datasets/brain_10x_5k_fastqs/ --id sc_mouse_brain_cellranger7 --transcriptome $mouse_genome_name  1> sc_mouse_brain_kallisto_standard_stdout.txt 2> sc_mouse_brain_kallisto_standard_stderr.txt</pre>


