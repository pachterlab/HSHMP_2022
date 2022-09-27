# STARSolo simulation

On tolva, everything is downloaded into: /home/dsullivan/benchmarking/starsolo

## Download and run simulation

<pre>git clone https://github.com/dobinlab/STARsoloManuscript
cd STARsoloManuscript
git checkout d5cfb6bf80861ccf9d19ccd99026d131c476d095 # (May 25, 2021 commit)

make -C samples

make -C exe # Gives errors: "make: *** [Makefile:26: kbpy_0.25.0] Error 1" but oh well

cat genomes/Makefile|sed 's/Genomes/genomes/' > genomes/Makefile2
mv genomes/Makefile2 genomes/Makefile
make -C genomes # Needed to edit Makefile (to make Genomes lower-case) above otherwise the files can't be found

# The above command still has error "make: *** [Makefile:23: human_CR_3.0.0/refdata-cellranger-GRCh38-3.0.0] Error 8", so do the following:

mv genomes/human_CR_3.0.0/labshare.cshl.edu/shares/dobin/dobin/STARsolo/Preprint/genomes/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa genomes/human_CR_3.0.0/
mv genomes/human_CR_3.0.0/labshare.cshl.edu/shares/dobin/dobin/STARsolo/Preprint/genomes/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf genomes/human_CR_3.0.0/annotations.gtf

samtools faidx genomes/human_CR_3.0.0/genome.fa

</pre>
