###################################################
#       HCA - Euclidean (S0 + CD + MD)            #                  
# Hierarchical Clustering Analysis (HCA)          #
# Protein Hotspot Analysis (FTMap / XDrugPy)      #
#                                                 #
# Author: Genesio Martins de Aguiar Neto (UFBA)   #
#                                                 #
# Purpose:                                        #
# Hierarchical clustering of FTMap-derived        #
# protein hotspots for structure-based            #
# drug design and comparative analysis.           #
#                                                 #
# Usage:                                          #
#   pymol -c <script_name>.py                     #
###################################################

from glob import glob
from os.path import expanduser, expandvars
import time, sys

from pymol import cmd as pm
from matplotlib import pyplot as plt
from xdrugpy.hotspots import load_ftmap, plot_euclidean_hca

BASE_PATH = "/home/jaimilson/Downloads/3B24_B_-_ADRIAN/*.pdb"
S0_CUTOFF = 20
OUTPUT = "3B24_B_S0_20_Euclidean.png"

EUCLIDEAN_SEL = f"*.K15_* AND p.S0>{S0_CUTOFF}"

def progress_bar(i, total, start, label):
    elapsed = time.time() - start
    avg = elapsed / i if i else 0
    eta = avg * (total - i)
    bar = "=" * int(30 * i / total)
    bar = bar.ljust(30, "-")
    sys.stdout.write(
        f"\r{label} [{bar}] {i}/{total} | ETA {eta:5.1f}s"
    )
    sys.stdout.flush()

if hasattr(pm, "undo_disable"):
    pm.undo_disable()

print("\n=== HCA: Euclidean (S0 + CD + MD) ===\n")

files = glob(expandvars(expanduser(BASE_PATH)))
if not files:
    raise RuntimeError("No FTMap PDB files found.")

start = time.time()
for i, f in enumerate(files, 1):
    load_ftmap(f)
    progress_bar(i, len(files), start, "Loading")

print("\n✔ FTMap loaded\n")

print("Running Euclidean clustering...")
start = time.time()

plot_euclidean_hca(
    EUCLIDEAN_SEL,
    properties=["S0", "CD", "MD"],
    normalize=True,
    linkage_method="ward",
    annotate=True
)

plt.title("Euclidean HCA (S0 + CD + MD) | S0 > 20")
plt.savefig(OUTPUT, dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ Euclidean HCA completed in {time.time()-start:.1f}s")
print(f"Saved: {OUTPUT}")
