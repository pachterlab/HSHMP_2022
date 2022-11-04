# STARSolo simulation

On tolva, everything is downloaded into: /home/dsullivan/benchmarking/starsolo

Replace the Makefile and Mf* files in the STARsoloManuscript directory with the ones here.

## Download and run simulation (Sep. 27, 2022)

<pre>git clone https://github.com/dobinlab/STARsoloManuscript
cd STARsoloManuscript
git checkout d5cfb6bf80861ccf9d19ccd99026d131c476d095 # (May 25, 2021 commit)

make -C samples

make -C exe # Gives errors: "make: *** [Makefile:26: kbpy_0.25.0] Error 1" but oh well

make -C exe gffread # Make sure we don't skip over installing gffread

cat genomes/Makefile|sed 's/Genomes/genomes/' > genomes/Makefile2
mv genomes/Makefile2 genomes/Makefile
make -C genomes # Needed to edit Makefile (to make Genomes lower-case) above otherwise the files can't be found

# The above command still has error "make: *** [Makefile:23: human_CR_3.0.0/refdata-cellranger-GRCh38-3.0.0] Error 8", so do the following:

mv genomes/human_CR_3.0.0/labshare.cshl.edu/shares/dobin/dobin/STARsolo/Preprint/genomes/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa genomes/human_CR_3.0.0/
mv genomes/human_CR_3.0.0/labshare.cshl.edu/shares/dobin/dobin/STARsolo/Preprint/genomes/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf genomes/human_CR_3.0.0/annotations.gtf

samtools faidx genomes/human_CR_3.0.0/genome.fa

./exe/gffread -w genomes/human_CR_3.0.0/transcripts.fa -g genomes/human_CR_3.0.0/genome.fa genomes/human_CR_3.0.0/annotations.gtf

awk '$3=="transcript" {gene=$0; gsub(/.*gene_id "/,"",gene); gsub(/".*/,"",gene); tr=$0; gsub(/.*transcript_id "/,"",tr); gsub(/".*/,"",tr); print tr "\t" gene }' \
	    genomes/human_CR_3.0.0/annotations.gtf > genomes/human_CR_3.0.0/transcript_to_gene.txt && \
	awk '{print $1 "\t" $2}' genomes/human_CR_3.0.0/transcript_to_gene.txt > genomes/human_CR_3.0.0/transcript_to_gene.2col.txt
awk '{if (!($2 in G)) print $2; G[$2]=0}'  genomes/human_CR_3.0.0/transcript_to_gene.txt > genomes/human_CR_3.0.0/genes.tsv

aa=$(pwd|sed 's/\//\\\//g') && cat Mf.common|sed 's/\/scratch\/dobin\/STAR\/STARsoloManuscript\//'"$aa"'\//' > Mf.common2
mv Mf.common2 Mf.common

cat data/Makefile|sed 's/Data/data/' > data/Makefile2
mv data/Makefile2 data/Makefile

make -C data

rm -f ./genomes/human_CR_3.0.0/genome_transcripts*

make -C sims
</pre>

# Create symlinks to executables

<pre>ln -s /home/dsullivan/kallisto-bf/build/src/kallisto exe/kallisto_0.49.0
ln -s /home/dsullivan/bustools/build/src/bustools exe/bustools_0.41.1
ln -s /home/kristjan/cellranger/cellranger-6.1.2/cellranger exe/CellRanger_6.1.6
ln -s /home/kristjan/cellranger/cellranger-7.0.1/cellranger exe/CellRanger_7.0.1
ln -s /home/dsullivan/salmon-1.9.0_linux_x86_64/bin/salmon exe/salmon_1.9.0
ln -s /home/dsullivan/.cargo/bin/alevin-fry exe/alevin-fry_0.8.0</pre>

# Create indices

<pre>genome_name="human_CR_3.0.0"
genome_file="genomes/$genome_name/genome.fa"
gtf_file="genomes/$genome_name/annotations.gtf"
n_threads="20"</pre>

## kallisto (kb-python)

<pre>out_dir="genomes/index/kallisto_0.49.0/$genome_name/standard_1"
mkdir -p $out_dir
kb ref -i $out_dir/index.idx --kallisto exe/kallisto_0.49.0 --workflow standard --overwrite -f1 $out_dir/f1 -g $out_dir/g $genome_file $gtf_file > $out_dir/log.txt 2>&1
exe/kallisto_0.49.0 index -i $out_dir/index.idx $out_dir/f1 # TODO: DELETE THIS ONCE WE FIGURE OUT WHY TF KB ISN'T WORKING!
</pre>

### off-list (in-progress)

<pre>out_dir="genomes/index/kallisto_0.49.0/$genome_name/standard_offlist_1"
mkdir -p $out_dir
cp ../extract_introns/* ./
./extract_introns.py --gtf <(cat $gtf_file|gzip) --fa <(cat $genome_file|gzip) --out $out_dir/introns.fa --union --diff
exe/kallisto_0.49.0 index -t 4 -b $out_dir/introns.fa -i $out_dir/index.idx genomes/index/kallisto_0.49.0/$genome_name/standard_1/f1
# ^TODO: REPLACE ABOVE WITH A KB REF COMMAND ONCE WE ALLOW OFFLIST IN KB REF
cp genomes/index/kallisto_0.49.0/$genome_name/standard_1/g $out_dir/g
# ^TODO: REPLACE ABOVE WITH A KB REF COMMAND ONCE WE ALLOW OFFLIST IN KB REF
</pre>


## STAR

<pre>out_dir="genomes/index/STAR_2.7.9a/$genome_name"
mkdir -p $out_dir/fullSA
exe/STAR_2.7.9a --runMode genomeGenerate --runThreadN $n_threads --genomeDir $out_dir/fullSA --genomeFastaFiles $genome_file --sjdbGTFfile $gtf_file > $out_dir/fullSA/log.txt 2>&1
mkdir -p $out_dir/sparseSA3
exe/STAR_2.7.9a --runMode genomeGenerate --runThreadN $n_threads --genomeDir $out_dir/sparseSA3 --genomeSAsparseD 3 --genomeFastaFiles $genome_file --sjdbGTFfile $gtf_file > $out_dir/sparseSA3/log.txt 2>&1</pre>

# Run simulations

Clear simulations:

<pre>rm -rf count/</pre>

Run on real data:

<pre>make real</pre>

Run on simulated data:

<pre>make sims_run</pre>

# Generate results from simulations

<pre>cp ../comparisons.py ./</pre>

Run python script with 7 arguments: 1) truth matrix, 2) truth genes list, 3) truth barcode list, 4) program matrix, 5) program genes list, 6) program barcodes list, 7) output file name.

## kallisto comparison

<pre>python3 comparisons.py samples/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/truth.mtx genomes/human_CR_3.0.0/genes.tsv data/whitelists/10Xv3 count/kallisto_0.49.0/human_CR_3.0.0/standard_1/default/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/counts_unfiltered/cells_x_genes.mtx count/kallisto_0.49.0/human_CR_3.0.0/standard_1/default/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/counts_unfiltered/cells_x_genes.genes.txt count/kallisto_0.49.0/human_CR_3.0.0/standard_1/default/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/counts_unfiltered/cells_x_genes.barcodes.txt results_sim_vs_kallisto.txt</pre>

## STAR comparison

<pre>python3 comparisons.py samples/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/truth.mtx genomes/human_CR_3.0.0/genes.tsv data/whitelists/10Xv3 count/STAR_2.7.9a/human_CR_3.0.0/fullSA/10X_noSAM_sims_mult_ENCODE/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/Solo.out/Gene/raw/matrix.mtx count/STAR_2.7.9a/human_CR_3.0.0/fullSA/10X_noSAM_sims_mult_ENCODE/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/Solo.out/Gene/raw/features.tsv count/STAR_2.7.9a/human_CR_3.0.0/fullSA/10X_noSAM_sims_mult_ENCODE/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/Solo.out/Gene/raw/barcodes.tsv results_sim_vs_star.txt --transpose</pre>

The output file will be formatted as follows:

* First line of output file is:
* * "# num_barcodes_in_intersection num_barcodes_in_simulation,num_barcodes_in_program"
* Second line of output file is:
* * "# num_genes_in_intersection num_genes_in_simulation,num_genes_in_program"
* Rest of output file repeats as the following 5 lines:
* * Five lines: 1) barcode, 2) simulation gene counts, 3) program gene counts, 4) false positive gene names w.r.t. simulation, 5) false negative gene names w.r.t. simulation

Here's an example of the first 17 lines of the output file:

<pre># 913780 6794880,913780
# 33538 33538,33538
AAACCCAAGAAACACT
1
1


AAACCCAAGAAACCCA
1,1
1,1


AAACCCAAGAAACCCG
1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1
ENSG00000072274
ENSG00000123975</pre>


