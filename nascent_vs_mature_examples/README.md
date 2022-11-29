Run the bbmap_sims steps first.

<pre>cd /home/dsullivan/benchmarking/kristjan/</pre>

<pre>kallisto="/home/dsullivan/benchmarking/kristjan/kallisto-bf/build/src/kallisto"
prev_kallisto="/home/dsullivan/kallisto/build/src/kallisto"
bustools="/home/dsullivan/bustools/build/src/bustools"

salmon="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/exe/salmon_1.9.0"
af="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/exe/alevin-fry_0.8.0"</pre>

Paths to indices: (TODO: )

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
kallisto_index_offlist="$main_path/genomes/index/kallisto_0.49.0/$genome_name/standard_offlist_1/index.idx"</pre>

Additional kallisto indices:

<pre>kallisto_index_offlist_mature="/home/kristjan/kallisto_bf_analysis/nascent_starsolo+threshold+offlist+polyA.idx"</pre>

<pre>nascent_fasta="/home/kristjan/kallisto_bf_analysis/partial_transcriptomes/nascent_starsolo_v2.fa"
cdna_fasta="/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/genomes/index/kallisto_0.49.0/human_CR_3.0.0/standard_1/f1"

$kallisto index -i kallisto_nascent_mature.idx -t $n_threads "$cdna_fasta" "$nascent_fasta"</pre>

<pre>out_dir="lamanno/"
mkdir -p $out_dir
kb ref --kallisto $prev_kallisto --workflow=lamanno -i $out_dir/index.idx -g $out_dir/g -f1 $out_dir/f1 -f2 $out_dir/f2 -c1 $out_dir/c1 -c2 $out_dir/c2 $genome_file $gtf_file</pre>

## Simulations (TODO)

### Reformat R1 files to have unique barcodes

<pre>nucleus_mature_r1="/home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/nuclear/Mature_S1_L001_R1_001.fastq.gz"
nucleus_nascent_r1="/home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/nuclear/Nascent_S1_L001_R1_001.fastq.gz"
cytoplasmic_mature_r1="/home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/cytoplasmic/Mature_S1_L001_R1_001.fastq.gz"
cytoplasmic_nascent_r1="/home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/cytoplasmic/Nascent_S1_L001_R1_001.fastq.gz"

paste <(zcat < $nucleus_mature_r1) <(./gen_unique_barcodes.sh $(zcat < $nucleus_mature_r1|wc -l))|awk '{if(NR%2==0) {print $2} else {print $1} }'|gzip > nucleus_mature_r1.fastq.gz
paste <(zcat < $nucleus_nascent_r1) <(./gen_unique_barcodes.sh $(zcat < $nucleus_nascent_r1|wc -l))|awk '{if(NR%2==0) {print $2} else {print $1} }'|gzip > nucleus_nascent_r1.fastq.gz
paste <(zcat < $cytoplasmic_mature_r1) <(./gen_unique_barcodes.sh $(zcat < $cytoplasmic_mature_r1|wc -l))|awk '{if(NR%2==0) {print $2} else {print $1} }'|gzip > cytoplasmic_mature_r1.fastq.gz
paste <(zcat < $cytoplasmic_nascent_r1) <(./gen_unique_barcodes.sh $(zcat < $cytoplasmic_nascent_r1|wc -l))|awk '{if(NR%2==0) {print $2} else {print $1} }'|gzip > cytoplasmic_nascent_r1.fastq.gz</pre>

### Set up paths

<pre>nucleus_mature_r1="nucleus_mature_r1.fastq.gz"
nucleus_mature_r2="/home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/nuclear/Mature_S1_L001_R2_001.fastq.gz"
nucleus_nascent_r1="nucleus_nascent_r1.fastq.gz"
nucleus_nascent_r2="/home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/nuclear/Nascent_S1_L001_R2_001.fastq.gz"

cytoplasmic_mature_r1="cytoplasmic_mature_r1.fastq.gz"
cytoplasmic_mature_r2="/home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/cytoplasmic/Mature_S1_L001_R2_001.fastq.gz"
cytoplasmic_nascent_r1="cytoplasmic_nascent_r1.fastq.gz"
cytoplasmic_nascent_r2="/home/kristjan/kallisto_bf_analysis/simulated_reads/10xV3_format/cytoplasmic/Nascent_S1_L001_R2_001.fastq.gz"</pre>

### Kallisto

<pre>results_file="sim_results_kallisto/results.txt"
out_dir_cytoplasmic="sim_results_kallisto/cytoplasmic"
out_dir_nucleus="sim_results_kallisto/nucleus"
mkdir -p $out_dir_cytoplasmic
mkdir -p $out_dir_nucleus
$kallisto bus -x 10xv3 --unstranded -i $kallisto_index_offlist -o $out_dir_cytoplasmic/mature/ -t $n_threads $cytoplasmic_mature_r1 $cytoplasmic_mature_r2
$kallisto bus -x 10xv3 --unstranded -i $kallisto_index_offlist -o $out_dir_cytoplasmic/nascent/ -t $n_threads $cytoplasmic_nascent_r1 $cytoplasmic_nascent_r2
$kallisto bus -x 10xv3 --unstranded -i $kallisto_index_offlist_mature -o $out_dir_nucleus/mature/ -t $n_threads $nucleus_mature_r1 $nucleus_mature_r2
$kallisto bus -x 10xv3 --unstranded -i $kallisto_index_offlist_mature -o $out_dir_nucleus/nascent/ -t $n_threads $nucleus_nascent_r1 $nucleus_nascent_r2

echo "sim n_reads_mapped n_reads_total" > $results_file
cat $out_dir_cytoplasmic/mature/run_info.json|grep "n_processed\|n_pseudoaligned"|cut -d' ' -f2|tr -d ','|xargs|awk '{print "kallisto_cytoplasmic_mature",$2,$1}' >> $results_file
cat $out_dir_cytoplasmic/nascent/run_info.json|grep "n_processed\|n_pseudoaligned"|cut -d' ' -f2|tr -d ','|xargs|awk '{print "kallisto_cytoplasmic_nascent",$2,$1}' >> $results_file
cat $out_dir_nucleus/mature/run_info.json|grep "n_processed\|n_pseudoaligned"|cut -d' ' -f2|tr -d ','|xargs|awk '{print "kallisto_nucleus_mature",$2,$1}' >> $results_file
cat $out_dir_nucleus/nascent/run_info.json|grep "n_processed\|n_pseudoaligned"|cut -d' ' -f2|tr -d ','|xargs|awk '{print "kallisto_nucleus_nascent",$2,$1}' >> $results_file
</pre>

Extract reads (nucleus mature that are mapped):

<pre>zcat < $nucleus_mature_r1|grep -B1 -f <($bustools text -p $out_dir_nucleus/mature/output.bus|cut -f1,2|tr -d '\t'|sort -u)|grep ^@ > tmp.txt
zcat < $nucleus_mature_r2|grep -A3 -f tmp.txt|grep -v ^\-\- > kallisto_reads_nucleus_mature_mapped.fastq
</pre>

Extract reads (cytoplasmic nascent that are mapped):

<pre>zcat < $cytoplasmic_nascent_r1|grep -B1 -f <($bustools text -p $out_dir_cytoplasmic/nascent/output.bus|cut -f1,2|tr -d '\t'|sort -u)|grep ^@ > tmp.txt
zcat < $cytoplasmic_nascent_r2|grep -A3 -f tmp.txt|grep -v ^\-\- > kallisto_reads_cytoplasmic_nascent_mapped.fastq
</pre>

### Salmon

<pre>results_file="sim_results_salmon/results.txt"
out_dir_cytoplasmic="sim_results_salmon/cytoplasmic"
out_dir_nucleus="sim_results_salmon/nucleus"
mkdir -p $out_dir_cytoplasmic
mkdir -p $out_dir_nucleus
$salmon alevin --chromiumV3 -p $n_threads -i $salmon_index_splici --tgMap $t2gFile_salmon_splici -l IU --rad -o $out_dir_cytoplasmic/mature/ -1 $cytoplasmic_mature_r1 -2 $cytoplasmic_mature_r2
$salmon alevin --chromiumV3 -p $n_threads -i $salmon_index_splici --tgMap $t2gFile_salmon_splici -l IU --rad -o $out_dir_cytoplasmic/nascent/ -1 $cytoplasmic_nascent_r1 -2 $cytoplasmic_nascent_r2
$salmon alevin --chromiumV3 -p $n_threads -i $salmon_index_splici --tgMap $t2gFile_salmon_splici -l IU --rad -o $out_dir_nucleus/mature/ -1 $nucleus_mature_r1 -2 $nucleus_mature_r2
$salmon alevin --chromiumV3 -p $n_threads -i $salmon_index_splici --tgMap $t2gFile_salmon_splici -l IU --rad -o $out_dir_nucleus/nascent/ -1 $nucleus_nascent_r1 -2 $nucleus_nascent_r2

echo "sim n_unspliced n_spliced n_reads_mapped n_reads_total" > $results_file
printf "%s %s %s %s %s\n" "salmon_cytoplasmic_mature" $($af view --rad $out_dir_cytoplasmic/mature/map.rad|grep "\-I"|cut -f4,5|sort -u|wc -l) $($af view --rad $out_dir_cytoplasmic/mature/map.rad|grep -v "\-I"|cut -f4,5|sort -u|wc -l) $(cat $out_dir_cytoplasmic/mature/aux_info/meta_info.json|grep "num_processed\|num_mapped"|cut -d':' -f2|tr -d ','|xargs|awk '{print $2,$1}') >> $results_file
printf "%s %s %s %s %s\n" "salmon_cytoplasmic_nascent" $($af view --rad $out_dir_cytoplasmic/nascent/map.rad|grep "\-I"|cut -f4,5|sort -u|wc -l) $($af view --rad $out_dir_cytoplasmic/nascent/map.rad|grep -v "\-I"|cut -f4,5|sort -u|wc -l) $(cat $out_dir_cytoplasmic/nascent/aux_info/meta_info.json|grep "num_processed\|num_mapped"|cut -d':' -f2|tr -d ','|xargs|awk '{print $2,$1}') >> $results_file
printf "%s %s %s %s %s\n" "salmon_nucleus_mature" $($af view --rad $out_dir_nucleus/mature/map.rad|grep "\-I"|cut -f4,5|sort -u|wc -l) $($af view --rad $out_dir_nucleus/mature/map.rad|grep -v "\-I"|cut -f4,5|sort -u|wc -l) $(cat $out_dir_nucleus/mature/aux_info/meta_info.json|grep "num_processed\|num_mapped"|cut -d':' -f2|tr -d ','|xargs|awk '{print $2,$1}') >> $results_file
printf "%s %s %s %s %s\n" "salmon_nucleus_nascent" $($af view --rad $out_dir_nucleus/nascent/map.rad|grep "\-I"|cut -f4,5|sort -u|wc -l) $($af view --rad $out_dir_nucleus/nascent/map.rad|grep -v "\-I"|cut -f4,5|sort -u|wc -l) $(cat $out_dir_nucleus/nascent/aux_info/meta_info.json|grep "num_processed\|num_mapped"|cut -d':' -f2|tr -d ','|xargs|awk '{print $2,$1}') >> $results_file
</pre>

Extract reads (nucleus mature that are mapped to intronic):

<pre>zcat < $nucleus_mature_r1|grep -B1 -f <($af view --rad $out_dir_nucleus/mature/map.rad|grep "\-I"|cut -f4,5|sed -r 's/:/\t/g'|cut -f2,4|tr -d "\t"|sort -u)|grep ^@ > tmp.txt
zcat < $nucleus_mature_r2|grep -A3 -f tmp.txt|grep -v ^\-\- > salmon_reads_nucleus_mature_unspliced.fastq
</pre>

### Salmon Sketch

<pre>results_file="sim_results_salmon_sketch/results.txt"
out_dir_cytoplasmic="sim_results_salmon_sketch/cytoplasmic"
out_dir_nucleus="sim_results_salmon_sketch/nucleus"
mkdir -p $out_dir_cytoplasmic
mkdir -p $out_dir_nucleus
$salmon alevin --chromiumV3 -p $n_threads -i $salmon_index_splici --tgMap $t2gFile_salmon_splici -l IU --rad -o $out_dir_cytoplasmic/mature/ -1 $cytoplasmic_mature_r1 -2 $cytoplasmic_mature_r2 --sketch
$salmon alevin --chromiumV3 -p $n_threads -i $salmon_index_splici --tgMap $t2gFile_salmon_splici -l IU --rad -o $out_dir_cytoplasmic/nascent/ -1 $cytoplasmic_nascent_r1 -2 $cytoplasmic_nascent_r2 --sketch
$salmon alevin --chromiumV3 -p $n_threads -i $salmon_index_splici --tgMap $t2gFile_salmon_splici -l IU --rad -o $out_dir_nucleus/mature/ -1 $nucleus_mature_r1 -2 $nucleus_mature_r2 --sketch
$salmon alevin --chromiumV3 -p $n_threads -i $salmon_index_splici --tgMap $t2gFile_salmon_splici -l IU --rad -o $out_dir_nucleus/nascent/ -1 $nucleus_nascent_r1 -2 $nucleus_nascent_r2 --sketch

echo "sim n_unspliced n_spliced n_reads_mapped n_reads_total" > $results_file
printf "%s %s %s %s %s\n" "salmon_cytoplasmic_mature" $($af view --rad $out_dir_cytoplasmic/mature/map.rad|grep "\-I"|cut -f4,5|sort -u|wc -l) $($af view --rad $out_dir_cytoplasmic/mature/map.rad|grep -v "\-I"|cut -f4,5|sort -u|wc -l) $(cat $out_dir_cytoplasmic/mature/aux_info/meta_info.json|grep "num_processed\|num_mapped"|cut -d':' -f2|tr -d ','|xargs|awk '{print $2,$1}') >> $results_file
printf "%s %s %s %s %s\n" "salmon_cytoplasmic_nascent" $($af view --rad $out_dir_cytoplasmic/nascent/map.rad|grep "\-I"|cut -f4,5|sort -u|wc -l) $($af view --rad $out_dir_cytoplasmic/nascent/map.rad|grep -v "\-I"|cut -f4,5|sort -u|wc -l) $(cat $out_dir_cytoplasmic/nascent/aux_info/meta_info.json|grep "num_processed\|num_mapped"|cut -d':' -f2|tr -d ','|xargs|awk '{print $2,$1}') >> $results_file
printf "%s %s %s %s %s\n" "salmon_nucleus_mature" $($af view --rad $out_dir_nucleus/mature/map.rad|grep "\-I"|cut -f4,5|sort -u|wc -l) $($af view --rad $out_dir_nucleus/mature/map.rad|grep -v "\-I"|cut -f4,5|sort -u|wc -l) $(cat $out_dir_nucleus/mature/aux_info/meta_info.json|grep "num_processed\|num_mapped"|cut -d':' -f2|tr -d ','|xargs|awk '{print $2,$1}') >> $results_file
printf "%s %s %s %s %s\n" "salmon_nucleus_nascent" $($af view --rad $out_dir_nucleus/nascent/map.rad|grep "\-I"|cut -f4,5|sort -u|wc -l) $($af view --rad $out_dir_nucleus/nascent/map.rad|grep -v "\-I"|cut -f4,5|sort -u|wc -l) $(cat $out_dir_nucleus/nascent/aux_info/meta_info.json|grep "num_processed\|num_mapped"|cut -d':' -f2|tr -d ','|xargs|awk '{print $2,$1}') >> $results_file
</pre>

## MYC

<pre>out_dir="myc"
mkdir -p $out_dir</pre>

<pre>cat $cdna_fasta|grep -A1 "ENSG00000136997"|grep -v ^\- > $out_dir/cdna.fa</pre>

<pre>cat $gtf_file|grep ENSG00000136997|grep $'\t'gene$'\t' > $out_dir/annotation.gtf</pre>

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
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exonoverlapintron
CCGCCACCGCCGGGCCCCGGCCGTCCCTGGCTCCCCTCCTGCCTCGAGAAGGGCAGGGCTTCTCAGAGGCTTGGCGGGAAAAAGAACGGAGGGAGGGATCGCGCTGAGTATAAAAGCCGGTTTTCGGGGCTTTATCTAACTCGCTGTAGT
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK" > $out_dir/reads.fq
</pre>

<pre>echo "@exon1
GGGGGGGGGGGGGGGGAAAAAATTTTTT
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exon2
GGGGGGGGGGGGGGGGAAAAAACCCCCC
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exon1exon2
GGGGGGGGGGGGGGGGAAAAAAGGGGGG
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exon1intron
GGGGGGGGGGGGGGGGTTTTTTAAAAAA
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@intron
GGGGGGGGGGGGGGGGCCCCCCAAAAAA
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@intronexon2
GGGGGGGGGGGGGGGGGGGGGGAAAAAA
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exonoverlapintron
GGGGGGGGGGGGGGGGCCCCCCCCCCCC
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK" > $out_dir/reads_bc_umi.fq
</pre>

<pre>echo "@exon1
GGGGGGGGGGGGGGGGAAAAAATTTTTT
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exon2
GGGGGGGGGGGGGGGGAAAAAACCCCCC
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exon1exon2
GGGGGGGGGGGGGGGGAAAAAAGGGGGG
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exon1intron
GGGGGGGGGGGGGGGGTTTTTTAAAAAA
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@intron
GGGGGGGGGGGGGGGGAAAAAATTTTTT
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@intronexon2
GGGGGGGGGGGGGGGGGGGGGGAAAAAA
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK
@exonoverlapintron
GGGGGGGGGGGGGGGGCCCCCCCCCCCC
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKK" > $out_dir/reads_bc_ambiguous_umi.fq
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

<pre>$kallisto bus -n -i $kallisto_index -o $out_dir/quant_cdna/ $out_dir/reads.fq</pre>
<pre>$bustools text -pf $out_dir/quant_cdna/output.bus</pre>

(Everything except the intron-only read, named intron, should map)

#### Full cdna index plus offlist

<pre>$kallisto bus -n -i $kallisto_index_offlist -o $out_dir/quant_offlist/ $out_dir/reads.fq</pre>
<pre>$bustools text -pf $out_dir/quant_offlist/output.bus</pre>

(Only exon1, exon2, exon1exon2, and exonoverlapintron should map; none of the nascent reads: intron, exon1intron, and intronexon2 should map)

#### Intron-only index (TODO note: uses an index from an old GTF)

<pre>$kallisto bus -n -i kallisto_index_introns/index.idx -o $out_dir/quant_introns/ $out_dir/reads.fq</pre>
<pre>$bustools text -pf $out_dir/quant_introns/output.bus</pre>

(exon1intron, intron, intronexon2, and exonoverlapintron should map)

#### Lamanno index

<pre>$prev_kallisto bus -n -i lamanno/index.idx -o $out_dir/quant_lamanno/ $out_dir/reads.fq</pre>
<pre>$bustools text -pf $out_dir/quant_lamanno/output.bus</pre>

(Everything except exon1intron and intronexon2 should map)

### Run salmon

#### Full cdna index (aka standard)

<pre>$salmon alevin -l ISR --rad -1 $out_dir/reads_bc_umi.fq -2 $out_dir/reads.fq --chromiumV3 -p $n_threads -o $out_dir/salmon_standard/ -i $salmon_index_standard --tgMap $t2gFile_salmon_standard
$af view --rad $out_dir/salmon_standard/map.rad</pre>

(Only exon1, exon2, exon1exon2, and exonoverlapintron should map; none of the nascent reads: intron, exon1intron, and intronexon2 should map)

#### Full cdna index (aka standard) with sketch

<pre>$salmon alevin -l ISR --rad -1 $out_dir/reads_bc_umi.fq -2 $out_dir/reads.fq --chromiumV3 -p $n_threads -o $out_dir/salmon_standard_sketch/ -i $salmon_index_standard --tgMap $t2gFile_salmon_standard --sketch
$af view --rad $out_dir/salmon_standard_sketch/map.rad</pre>

(Everything except the intron-only read, named intron, should map)

#### Splici index

<pre>$salmon alevin -l ISR --rad -1 $out_dir/reads_bc_umi.fq -2 $out_dir/reads.fq --chromiumV3 -p $n_threads -o $out_dir/salmon_splici/ -i $salmon_index_splici --tgMap $t2gFile_salmon_splici
$af view --rad $out_dir/salmon_splici/map.rad
$af generate-permit-list --force-cells 1 -d fw -i $out_dir/salmon_splici/ -o $out_dir/salmon_splici/
$af collate -t $n_threads -i $out_dir/salmon_splici/ -r $out_dir/salmon_splici/
$af quant --resolution cr-like -t $n_threads -i $out_dir/salmon_splici/ -o $out_dir/salmon_splici/quant/ --use-mtx --tg-map $t2gFile_salmon_splici
cat $out_dir/salmon_splici/quant/alevin/quants_mat.mtx</pre>

(exon1, exon2, and exon1exon2 will map to spliced target ENSG00000136997; the nascent reads: intron, exon1intron, and intronexon2 will map to intron target: ENSG00000136997-I1, and be considered unspliced; exonoverlapintron is considered ambiguous)

##### Splici index + an exon/intron ambiguous UMI

<pre>$salmon alevin -l ISR --rad -1 $out_dir/reads_bc_ambiguous_umi.fq -2 $out_dir/reads.fq --chromiumV3 -p $n_threads -o $out_dir/salmon_splici_ambiguous/ -i $salmon_index_splici --tgMap $t2gFile_salmon_splici
$af view --rad $out_dir/salmon_splici_ambiguous/map.rad
$af generate-permit-list --force-cells 1 -d fw -i $out_dir/salmon_splici_ambiguous/ -o $out_dir/salmon_splici_ambiguous/
$af collate -t $n_threads -i $out_dir/salmon_splici_ambiguous/ -r $out_dir/salmon_splici_ambiguous/
$af quant --resolution cr-like -t $n_threads -i $out_dir/salmon_splici_ambiguous/ -o $out_dir/salmon_splici_ambiguous/quant/ --use-mtx --tg-map $t2gFile_salmon_splici
cat $out_dir/salmon_splici_ambiguous/quant/alevin/quants_mat.mtx</pre>

The UMI AAAAAATTTTTT, which is associated with the two reads: exon1 and intron, will be classified as "ambiguous". The read associated with exonoverlapintron will also be classified as "ambiguous". The remaining 4 reads are either spliced (exon2 and exon1exon2) or unspliced (exon1intron and intronexon2).
