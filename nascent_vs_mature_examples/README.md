Run the bbmap_sims steps first.

<pre>cd /home/dsullivan/benchmarking/kristjan/</pre>

<pre>kallisto="/home/dsullivan/benchmarking/kristjan/kallisto-bf/build/src/kallisto"
prev_kallisto="/home/dsullivan/kallisto/build/src/kallisto"
bustools="/home/dsullivan/bustools/build/src/bustools"</pre>

## MYC

<pre>out_dir="myc"
mkdir -p $out_dir</pre>

<pre>cat kallisto_index_cdna/f1|grep -A1 "ENSG00000136997"|grep -v ^\- > $out_dir/cdna.fa</pre>

<pre>gtf_file="/home/kristjan/kallisto_bf_analysis/hg38_dna/Homo_sapiens.GRCh38.104.gtf.gz"
zcat < $gtf_file|grep ENSG00000136997|grep $'\t'gene$'\t' > $out_dir/annotation.gtf</pre>

<pre>echo "@exon1
GCTCGCCCAAGTCCTGCGCCTCGCAAGACTCCAGCGCCTTCTCTCCGTCCTCGGATTCTCTGCTCTCCTCGACGGAGTCCTCCCCGCAGGGCAGCCCCGAGCCCCTGGTGCTCCATGAGGAGACACCGCCCACCACCAGCAGCGACTCTG
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exon2
AGGAGGAACAAGAAGATGAGGAAGAAATCGATGTTGTTTCTGTGGAAAAGAGGCAGGCTCCTGGCAAAAGGTCAGAGTCTGGATCACCTTCTGCTGGAGGCCACAGCAAACCTCCTCACAGCCCACTGGTCCTCAAGAGGTGCCACGTCT
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exon1exon2
AGTCCTCCCCGCAGGGCAGCCCCGAGCCCCTGGTGCTCCATGAGGAGACACCGCCCACCACCAGCAGCGACTCTGAGGAGGAACAAGAAGATGAGGAAGAAATCGATGTTGTTTCTGTGGAAAAGAGGCAGGCTCCTGGCAAAAGGTCAG
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exon1intron
AGTCCTCCCCGCAGGGCAGCCCCGAGCCCCTGGTGCTCCATGAGGAGACACCGCCCACCACCAGCAGCGACTCTGGTAAGCGAAGCCCGCCCAGGCCTGTCAAAAGTGGGCGGCTGGATACCTTTCCCATTTTCATTGGCAGCTTATTTA
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
@intron
GGTGCTTGGGAATGTGCTTTGCTTTGGGTGTGTCCAAAGCCTCATTAAGTCTTAGGTAAGAATTGGCATCAATGTCCTATCCTGGGAAGTTGCACTTTTCTTGTCCATGCCATAACCCAGCTGTCTTTCCCTTTATGAGACTCTTACCTT
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
@intronexon2
TACAGCATTAATCTGGTAATTGATTATTTTAATGTAACCTTGCTAAAGGAGTGATTTCTATTTCCTTTCTTAAAGAGGAGGAACAAGAAGATGAGGAAGAAATCGATGTTGTTTCTGTGGAAAAGAGGCAGGCTCCTGGCAAAAGGTCAG
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK" > $out_dir/reads.fq
</pre>

* Mature: exon1exon2
* Nascent: exon1intron, intron, intronexon2
* Ambiguous: exon1, exon2

### Run kallisto

#### Myc cdna-only index

<pre>$kallisto index -i $out_dir/$out_dir.idx $out_dir/cdna.fa</pre>
<pre>$kallisto bus -n -i $out_dir/$out_dir.idx -o $out_dir/quant/ $out_dir/reads.fq</pre>
<pre>$bustools text -pf $out_dir/quant/output.bus</pre>

(Everything except the intron-only read, named intron, should map)

#### Full cdna index

<pre>$kallisto bus -n -i kallisto_index_cdna/index.idx -o $out_dir/quant_cdna/ $out_dir/reads.fq</pre>
<pre>$bustools text -pf $out_dir/quant_cdna/output.bus</pre>

(Everything except the intron-only read, named intron, should map)

#### Full cdna index plus offlist

<pre>$kallisto bus -n -i kallisto_index_offlist/index.idx -o $out_dir/quant_offlist/ $out_dir/reads.fq</pre>
<pre>$bustools text -pf $out_dir/quant_offlist/output.bus</pre>

(Everything except the intron-only read, named intron, should map)

#### Intron-only index

<pre>$kallisto bus -n -i kallisto_index_introns/index.idx -o $out_dir/quant_introns/ $out_dir/reads.fq</pre>
<pre>$bustools text -pf $out_dir/quant_introns/output.bus</pre>

(exon1intron, intron, and intronexon2 should map)

#### Lamanno index

<pre>$prev_kallisto bus -n -i kallisto_index_lamanno/index.idx -o $out_dir/quant_lamanno/ $out_dir/reads.fq</pre>
<pre>$bustools text -pf $out_dir/quant_lamanno/output.bus</pre>

(Everything except exon1intron and intronexon2 should map)

