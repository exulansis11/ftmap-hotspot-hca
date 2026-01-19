###################################################
#      HCA - Residue Overlap Similarity           #
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

# ================= USER PARAMETERS =================
BASE_PATH = "/home/jaimilson/Downloads/3B24_B_-_ADRIAN/*.pdb"
S0_CUTOFF = 20
OUTPUT = "S0_20_Overlap.png"

HOTSPOT_SEL = f"((*.B* *.D*) AND p.S0>{S0_CUTOFF})"

# ================= Progress bar =================
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

# ================= PyMOL optimization =================
if hasattr(pm, "undo_disable"):
    pm.undo_disable()

print("\n=== HCA: Residue Overlap Similarity ===\n")

# ================= Load FTMap structures =================
files = glob(expandvars(expanduser(BASE_PATH)))
if not files:
    raise RuntimeError("No FTMap PDB files found. Check BASE_PATH.")

start = time.time()
for i, f in enumerate(files, 1):
    load_ftmap(f)
    progress_bar(i, len(files), start, "Loading")

print("\n✔ FTMap structures loaded\n")

# ================= Global alignment =================
objs = pm.get_object_list()
ref = objs[0]

for i, obj in enumerate(objs[1:], 1):
    pm.align(obj, ref)

print("✔ Structures aligned\n")

# ================= Hotspot validation =================
hotspots = pm.get_object_list(HOTSPOT_SEL)
if len(hotspots) < 2:
    raise RuntimeError(
        f"Only {len(hotspots)} hotspots found with S0 > {S0_CUTOFF}. "
        "At least two are required for overlap analysis."
    )

print(f"✔ {len(hotspots)} hotspots selected\n")

# ================= HCA - Overlap =================
print("Running residue overlap clustering...")
start = time.time()

plot_pairwise_hca(
    HOTSPOT_SEL,
    function="overlap",
    align=True,
    radius=5.0,
    annotate=True
)

plt.title("Residue Overlap Similarity | S0 > 20")
plt.savefig(OUTPUT, dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ Overlap HCA completed in {time.time() - start:.1f}s")
print(f"Saved: {OUTPUT}\n")
