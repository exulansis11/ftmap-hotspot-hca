###################################################
#                 HCA_HO                          #
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

# =================================================
# USER PARAMETERS
# =================================================
BASE_PATH = "[path]*.pdb"
S0_CUTOFF = [number]
OUTPUT = "[FILENAME]_HO.png"

HOTSPOT_SEL = f"((*.B* *.D*) AND p.S0>{S0_CUTOFF})"

# =================================================
# Progress bar utility
# =================================================
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

# =================================================
# PyMOL optimization
# =================================================
if hasattr(pm, "undo_disable"):
    pm.undo_disable()

print("\n=== HCA: Hotspot Overlap (HO) ===\n")

# =================================================
# Load FTMap structures
# =================================================
files = glob(expandvars(expanduser(BASE_PATH)))
if not files:
    raise RuntimeError("No FTMap PDB files found. Check BASE_PATH.")

start = time.time()
for i, f in enumerate(files, 1):
    load_ftmap(f)
    progress_bar(i, len(files), start, "Loading FTMap")

print("\n✔ FTMap structures loaded\n")

# =================================================
# Align ONLY protein structures (CRITICAL FIX)
# =================================================
print("Aligning protein structures only...")

proteins = [
    obj for obj in pm.get_object_list()
    if pm.count_atoms(f"{obj} and polymer.protein") > 0
]

if len(proteins) < 2:
    raise RuntimeError("Not enough protein structures for alignment.")

ref = proteins[0]
align_start = time.time()

for i, obj in enumerate(proteins[1:], 1):
    try:
        pm.align(obj, ref)
    except Exception as e:
        print(f"\n⚠️ Alignment failed for {obj}: {e}")
    progress_bar(i, len(proteins) - 1, align_start, "Aligning proteins")

print("\n✔ Protein alignment completed\n")

# =================================================
# Hotspot validation
# =================================================
hotspots = pm.get_object_list(HOTSPOT_SEL)
if len(hotspots) < 2:
    raise RuntimeError(
        f"Only {len(hotspots)} hotspots found with S0 > {S0_CUTOFF}. "
        "At least two are required for HO analysis."
    )

print(f"✔ {len(hotspots)} hotspots selected for analysis\n")

# =================================================
# HCA - Hotspot Overlap (HO)
# =================================================
print("Running Hotspot Overlap (HO) clustering...")
start = time.time()

plot_pairwise_hca(
    HOTSPOT_SEL,
    function="ho",
    align=False,      # already aligned proteins
    radius=1.5,
    annotate=True
)

plt.title("Hotspot Overlap (HO) | S0 > 20")
plt.savefig(OUTPUT, dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ HO clustering completed in {time.time() - start:.1f}s")
print(f"Saved: {OUTPUT}")
print("\n=== HO ANALYSIS FINISHED SUCCESSFULLY ===\n")
