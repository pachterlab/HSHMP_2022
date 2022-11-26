#!/usr/bin/env python3

paths = [
    '/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/count/salmon-alevin-fry_1.9.0_0.8.0/human_CR_3.0.0/splici/rad/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/gpl_knee/quant_cr-like-em/alevin/quants_mat_rows.txt',
    '/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/count/salmon-alevin-fry_1.9.0_0.8.0/human_CR_3.0.0/splici/rad/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/gpl_knee/quant_cr-like/alevin/quants_mat_rows.txt',
    '/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/count/STAR_2.7.9a/human_CR_3.0.0/fullSA/10X_noSAM_sims_mult_ENCODE/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/Solo.out/Gene/filtered/barcodes.tsv'
]

bc = set(b.rstrip('\n') for b in open(paths[0]))
for p in paths[1:]:
    bc = bc & set(b.rstrip('\n') for b in open(p))

for b in bc:
    print(b)
