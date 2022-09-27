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
</pre>
