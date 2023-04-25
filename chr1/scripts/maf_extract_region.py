from Bio import AlignIO
from Bio.AlignIO import MafIO
import sys

name = sys.argv[1]
target_seqname = sys.argv[2]
start_coord = int(sys.argv[3])
end_coord = int(sys.argv[4])

# Load mafindex
idx = MafIO.MafIndex(f'{name}.mafindex', f'{name}.maf', target_seqname)
# Parse the alignment
results = idx.search([start_coord], [end_coord])

output_handle = open(f"{name}.region.maf", "w")
AlignIO.write(results, output_handle, "maf")
