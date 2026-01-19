###################################################
#       HCA - Residue Jaccard Similarity          #
#                                                 #
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
from xdrugpy.hotspots import load_ftmap, plot_pairwise_hca

BASE_PATH = "/home/jaimilson/Downloads/3B24_B_-_ADRIAN/*.pdb"
S0_CUTOFF = 20
OUTPUT = "3B24_B_S0_20_Jaccard.png"

HOTSPOT_SEL = f"((*.B* *.D*) AND p.S0>{S0_CUTOFF})"

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

print("\n=== HCA: Residue Jaccard ===\n")

files = glob(expandvars(expanduser(BASE_PATH)))
if not files:
    raise RuntimeError("No FTMap PDB files found.")

start = time.time()
for i, f in enumerate(files, 1):
    load_ftmap(f)
    progress_bar(i, len(files), start, "Loading")

print("\n✔ FTMap loaded\n")

objs = pm.get_object_list()
ref = objs[0]
for obj in objs[1:]:
    pm.align(obj, ref)

hotspots = pm.get_object_list(HOTSPOT_SEL)
if len(hotspots) < 2:
    raise RuntimeError("Not enough hotspots for Jaccard analysis")

print("Running Jaccard clustering...")
start = time.time()

plot_pairwise_hca(
    HOTSPOT_SEL,
    function="jaccard",
    align=True,
    radius=5.0,
    annotate=True
)

plt.title("Residue Jaccard Similarity | S0 > 20")
plt.savefig(OUTPUT, dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ Jaccard completed in {time.time()-start:.1f}s")
print(f"Saved: {OUTPUT}")
