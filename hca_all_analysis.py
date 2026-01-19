###################################################
#     Comprehensive Hotspot Analysis Pipeline     #
#          FTMap / XDrugPy / PyMOL                #
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

from os.path import expanduser, expandvars
from glob import glob
import time
import sys

from pymol import cmd as pm
from matplotlib import pyplot as plt

from xdrugpy.hotspots import (
    load_ftmap,
    plot_pairwise_hca,
    plot_euclidean_hca
)

# =================================================
# USER PARAMETERS
# =================================================
BASE_PATH = "[path]/*.pdb"
S0_CUTOFF = [number]
ALIGN_STRUCTURES = True
OUTPUT_PREFIX = "[FILENAME]"

# =================================================
# Utility: progress bar with ETA
# =================================================
def progress_bar(current, total, start_time, prefix=""):
    elapsed = time.time() - start_time
    avg = elapsed / current if current > 0 else 0
    eta = avg * (total - current)

    bar_len = 30
    filled = int(bar_len * current / total)
    bar = "=" * filled + "-" * (bar_len - filled)

    sys.stdout.write(
        f"\r{prefix} [{bar}] {current}/{total} | "
        f"Elapsed: {elapsed:6.1f}s | ETA: {eta:6.1f}s"
    )
    sys.stdout.flush()

# =================================================
# PyMOL optimization
# =================================================
if hasattr(pm, "undo_disable"):
    pm.undo_disable()

global_start = time.time()
print("\n=== HOTSPOT ANALYSIS PIPELINE STARTED ===\n")

# =================================================
# Load FTMap structures
# =================================================
files = glob(expandvars(expanduser(BASE_PATH)))
n_files = len(files)

if n_files == 0:
    raise RuntimeError("No FTMap PDB files found. Check BASE_PATH.")

print(f"Found {n_files} FTMap structures\n")

load_start = time.time()
for i, file in enumerate(files, 1):
    load_ftmap(file)
    progress_bar(i, n_files, load_start, prefix="Loading FTMap files")

print("\n✔ FTMap structures loaded\n")

# =================================================
# Global structural alignment
# =================================================
if ALIGN_STRUCTURES:
    print("Aligning structures...")
    align_start = time.time()

    objs = pm.get_object_list()
    ref = objs[0]

    for i, obj in enumerate(objs[1:], 1):
        try:
            pm.align(obj, ref)
        except:
            pass
        progress_bar(i, len(objs) - 1, align_start, prefix="Aligning")

    print("\n✔ Alignment completed\n")

# =================================================
# Hotspot selection (PyMOL selector VALID)
# =================================================
HOTSPOT_SEL = f"((*.B* *.D*) AND p.S0>{S0_CUTOFF})"
print(f"Using hotspot selection:\n  {HOTSPOT_SEL}\n")

hotspots = pm.get_object_list(HOTSPOT_SEL)
if len(hotspots) < 2:
    raise RuntimeError(
        f"Only {len(hotspots)} hotspots found with S0 > {S0_CUTOFF}. "
        "At least two are required for HCA."
    )

print(f"✔ {len(hotspots)} hotspots selected for analysis\n")

# =================================================
# 1) Spatial similarity — Hotspot Overlap (HO)
# =================================================
print("Running HCA: Hotspot Overlap (HO)...")
start = time.time()

plot_pairwise_hca(
    HOTSPOT_SEL,
    function="ho",
    align=False,
    radius=1.5,
    color_threshold=0,
    hide_threshold=False,
    annotate=True
)

plt.title("Hotspot Overlap (HO) | S0 > 20")
plt.savefig(f"{OUTPUT_PREFIX}_01_HO_dendrogram.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ HO completed in {time.time() - start:.1f}s\n")

# =================================================
# 2) Functional similarity — Jaccard
# =================================================
print("Running HCA: Residue Jaccard similarity...")
start = time.time()

plot_pairwise_hca(
    HOTSPOT_SEL,
    function="jaccard",
    align=True,
    radius=5.0,
    color_threshold=0,
    hide_threshold=False,
    annotate=True
)

plt.title("Residue Jaccard Similarity | S0 > 20")
plt.savefig(f"{OUTPUT_PREFIX}_02_Jaccard_dendrogram.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ Jaccard completed in {time.time() - start:.1f}s\n")

# =================================================
# 3) Functional core similarity — Overlap
# =================================================
print("Running HCA: Residue Overlap similarity...")
start = time.time()

plot_pairwise_hca(
    HOTSPOT_SEL,
    function="overlap",
    align=True,
    radius=5.0,
    color_threshold=0,
    hide_threshold=False,
    annotate=True
)

plt.title("Residue Overlap Similarity | S0 > 20")
plt.savefig(f"{OUTPUT_PREFIX}_03_Overlap_dendrogram.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ Overlap completed in {time.time() - start:.1f}s\n")

# =================================================
# 4) Multivariate Euclidean HCA (S0 + CD + MD)
# =================================================
print("Running multivariate Euclidean HCA (S0 + CD + MD)...")
start = time.time()

plot_euclidean_hca(
    f"*.K15_* AND p.S0>{S0_CUTOFF}",
    properties=["S0", "CD", "MD"],
    normalize=True,
    linkage_method="ward",
    hide_threshold=False,
    annotate=True
)

plt.title("Euclidean HCA (S0 + CD + MD) | S0 > 20")
plt.savefig(f"{OUTPUT_PREFIX}_04_Euclidean_HCA.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ Euclidean HCA completed in {time.time() - start:.1f}s\n")

# =================================================
# Final report
# =================================================
total_time = time.time() - global_start
print("=== PIPELINE COMPLETED SUCCESSFULLY ===")
print(f"Total execution time: {total_time:.1f} seconds")
print("Generated files:")
print(f"  {OUTPUT_PREFIX}_01_HO_dendrogram.png")
print(f"  {OUTPUT_PREFIX}_02_Jaccard_dendrogram.png")
print(f"  {OUTPUT_PREFIX}_03_Overlap_dendrogram.png")
print(f"  {OUTPUT_PREFIX}_04_Euclidean_HCA.png\n")

