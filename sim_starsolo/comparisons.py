import scipy
import numpy as np
from scipy import io
from scipy.sparse import csr_matrix
from scipy import stats
import sys

# Simulation paths:

sim_truth_mtx_MultiGeneNo = sys.argv[1]
sim_truth_genes = sys.argv[2]
sim_truth_barcodes = sys.argv[3]

# Program run paths:

program_mtx_standard = sys.argv[4]
program_gene_standard = sys.argv[5]
program_barcodes_standard = sys.argv[6]

# Output file path:

outfilename = sys.argv[7]

# Extract genes and barcodes lists

genes = np.array([line.rstrip().split('.', 1)[0] for line in open(sim_truth_genes)])
barcodes = np.array([line.rstrip() for line in open(sim_truth_barcodes)])
genes_ = np.array([line.rstrip().split('.', 1)[0] for line in open(program_gene_standard)])
barcodes_ = np.array([line.rstrip() for line in open(program_barcodes_standard)])

# Intersect genes and barcodes lists

intersection_barcodes = np.intersect1d(barcodes, barcodes_);
sorter = np.argsort(barcodes)
bc1 = sorter[np.searchsorted(barcodes, intersection_barcodes, sorter=sorter)]
sorter = np.argsort(barcodes_)
bc2 = sorter[np.searchsorted(barcodes_, intersection_barcodes, sorter=sorter)]

intersection_genes = np.intersect1d(genes, genes_);
sorter = np.argsort(genes)
genes1 = sorter[np.searchsorted(genes, intersection_genes, sorter=sorter)]
sorter = np.argsort(genes_)
genes2 = sorter[np.searchsorted(genes_, intersection_genes, sorter=sorter)]

# Read matrix files (note: simulation and program have different matrix formats so we read them differently)

program_data = scipy.io.mmread(program_mtx_standard)
program_data_csr = program_data.tocsr()
sim_data = np.array([])
with open(sim_truth_mtx_MultiGeneNo, 'r') as f:
    sim_data = np.array([[int(num) for num in line.split(' ')] for line in f])
sparse_mtx_sim_data = csr_matrix((sim_data[:,2], (sim_data[:,1]-1, sim_data[:,0]-1)), shape = (barcodes.size, genes.size))

# Output file

n = 0
outfile = open(outfilename, "w")
# First line of output file is: "# num_barcodes_in_intersection num_barcodes_in_simulation,num_barcodes_in_program"
# Second line of output file is: "# num_genes_in_intersection num_genes_in_simulation,num_genes_in_program"
outfile.write("# " + str(len(intersection_barcodes)) + " " + str(len(barcodes)) + "," + str(len(barcodes_)) + "\n")
outfile.write("# " + str(len(intersection_genes)) + " " + str(len(genes)) + "," + str(len(genes_)) + "\n")
# Rest of output file repeats as the following 5 lines: 1) barcode, 2) simulation gene counts, 3) program gene counts, 4) false positive gene names w.r.t. simulation, 5) false negative gene names w.r.t. simulation
for i in range(len(intersection_barcodes)):
    x = [program_data_csr[bc2[i],genes2].toarray()[0], sparse_mtx_sim_data[bc1[i],genes1].toarray()[0]]
    line0 = intersection_barcodes[n]
    n += 1
    where_i = np.where(x[0]+x[1] > 0)[0]
    line1 = (",".join(map(str, x[1][where_i].astype(int))))
    line2 = (",".join(map(str, x[0][where_i].astype(int))))
    where_fp = np.where((x[0]+x[1] > 0) & (x[1] == 0))[0]
    where_fn = np.where((x[0]+x[1] > 0) & (x[0] == 0))[0]
    line3 = (",".join(intersection_genes[where_fp]))
    line4 = (",".join(intersection_genes[where_fn]))
    outfile.write(line0 + "\n" + line1 + "\n" + line2 + "\n" + line3 + "\n" + line4 + "\n")
    if n % 10000 == 0:
        print("Processed: " + n + " out of " + len(intersection_barcodes) + " barcodes")
outfile.close()
print("DONE")
