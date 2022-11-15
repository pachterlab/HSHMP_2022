# Navigate to base current working directory

<pre>cd /home/dsullivan/benchmarking/kristjan/</pre>

# Clone kallisto-bf report

<pre>git clone -b https://github.com/Yenaled/kallisto-bf.git
cd kallisto-bf/ext/htslib && autoheader && autoconf && cd ../../
mkdir -p build && cd build
cmake .. && make && cd ../../</pre>

# Python scripts to extract introns

<pre>cp ../extract_introns.py ./
cp ../utils.py ./
./extract_introns.py --gtf $gtf_file --fa $genome_file --out ./introns_.fa --union --diff
cat ./introns_.fa|awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' > ./introns.fa</pre>

# Define paths

Set up paths to kallisto (current version and previous version), genome FASTA file, and GTF file.

<pre>kallisto="/home/dsullivan/benchmarking/kristjan/kallisto-bf/build/src/kallisto"
prev_kallisto="/home/dsullivan/kallisto/build/src/kallisto"
genome_file="/home/kristjan/kallisto_bf_analysis/hg38_dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf_file="/home/kristjan/kallisto_bf_analysis/hg38_dna/Homo_sapiens.GRCh38.104.gtf.gz"</pre>

# Create kallisto indices

## cDNA-only index

<pre>out_dir="kallisto_index_cdna/"
mkdir -p $out_dir
kb ref -i $out_dir/index.idx --kallisto $kallisto --workflow standard --overwrite -f1 $out_dir/f1 -g $out_dir/g $genome_file $gtf_file > $out_dir/log.txt 2>&1
$kallisto index -i $out_dir/index.idx $out_dir/f1 # TODO: DELETE THIS ONCE WE FIGURE OUT WHY TF KB ISN'T WORKING!</pre>

## cDNA + Off-list introns

<pre>prev_dir="$out_dir"
out_dir="kallisto_index_offlist/"
mkdir -p $out_dir
cp $prev_dir/* $out_dir
$kallisto index -t 44 -b ./introns.fa -i $out_dir/index.idx $prev_dir/f1</pre>

## Intron-only index

<pre>out_dir="kallisto_index_introns/"
mkdir -p $out_dir
$kallisto index -t 44 -i $out_dir/index.idx ./introns.fa</pre>

## Lamanno index (previous kallisto version)

<pre>prev_dir="$out_dir"
out_dir="kallisto_index_lamanno/"
mkdir -p $out_dir
kb ref --kallisto $prev_kallisto --workflow=lamanno -i $out_dir/index.idx -g $out_dir/g -f1 $out_dir/f1 -f2 $out_dir/f2 -c1 $out_dir/c1 -c2 $out_dir/c2 $genome_file $gtf_file</pre>

# Define paths to the four simulations

<pre>f1="/home/kristjan/kallisto_bf_analysis/simulated_reads/hg38_cDNA+introns_5000000_reads_20221110.fq"
f2="/home/kristjan/kallisto_bf_analysis/simulated_reads/hg38_cDNA_only_5000000_reads_20221110.fq"
f3="/home/kristjan/kallisto_bf_analysis/simulated_reads/hg38_introns_only_5000000_reads_20221110.fq"
f4="/home/kristjan/kallisto_bf_analysis/simulated_reads/hg38_nascent_only_5000000_reads_20221110.fq"</pre>

# Run the various kallisto indices on the four simulations

<pre>out_dir="kallisto_index_cdna/"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f1_out/ --single -l 1 -s 1 --single-overhang "$f1"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f2_out/ --single -l 1 -s 1 --single-overhang "$f2"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f3_out/ --single -l 1 -s 1 --single-overhang "$f3"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f4_out/ --single -l 1 -s 1 --single-overhang "$f4"

out_dir="kallisto_index_offlist/"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f1_out/ --single -l 1 -s 1 --single-overhang "$f1"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f2_out/ --single -l 1 -s 1 --single-overhang "$f2"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f3_out/ --single -l 1 -s 1 --single-overhang "$f3"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f4_out/ --single -l 1 -s 1 --single-overhang "$f4"

out_dir="kallisto_index_introns/"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f1_out/ --single -l 1 -s 1 --single-overhang "$f1"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f2_out/ --single -l 1 -s 1 --single-overhang "$f2"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f3_out/ --single -l 1 -s 1 --single-overhang "$f3"
$kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f4_out/ --single -l 1 -s 1 --single-overhang "$f4"

out_dir="kallisto_index_lamanno/"
$prev_kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f1_out/ --single -l 1 -s 1 --single-overhang "$f1"
$prev_kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f2_out/ --single -l 1 -s 1 --single-overhang "$f2"
$prev_kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f3_out/ --single -l 1 -s 1 --single-overhang "$f3"
$prev_kallisto quant -i $out_dir/index.idx -t 16 -o $out_dir/f4_out/ --single -l 1 -s 1 --single-overhang "$f4"</pre>

# Analyze number of reads mapped

Each output folder contains a run_info.json file that tells you how many reads were mapped. For example, to see the number of reads mapped for the f2 FASTQ file (hg38_cDNA_only_5000000_reads_20221110.fq) for the cDNA + Off-list introns kallisto index (kallisto_index_offlist/index.idx), you should go to **kallisto_index_offlist/f2_out/run_info.json** 

