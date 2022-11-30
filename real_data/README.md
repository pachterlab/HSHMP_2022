<pre>cd /home/dsullivan/benchmarking/kristjan/</pre>

<pre>kallisto="/home/dsullivan/benchmarking/kristjan/kallisto-bf/build/src/kallisto"
prev_kallisto="/home/dsullivan/kallisto/build/src/kallisto"
bustools="/home/dsullivan/bustools/build/src/bustools"

salmon="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/exe/salmon_1.9.0"
af="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/exe/alevin-fry_0.8.0"

star="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript//exe/STAR_2.7.9a"
</pre>

Paths to indices:

<pre>main_path="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript"
genome_name="human_CR_3.0.0"
genome_file="$main_path/genomes/$genome_name/genome.fa"
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
ln -s /home/kristjan/data/liver_andrews/ ./data/</pre>

## CellRanger7 Index

<pre>cellranger="/home/kristjan/cellranger/cellranger-7.0.1/cellranger"
$cellranger mkref --genome=$genome_name --fasta=$genome_file --genes=$gtf_file --nthreads=$n_threads</pre>

## CellRanger Run

<pre>$cellranger count --fastqs cellranger_cytoplasmic/ --sample Mature --id cellranger_cytoplasmic_mature --transcriptome human_CR_3.0.0  --chemistry SC3Pv3
</pre>
