#!/usr/bin/env python3

import random
import time
import os
from datetime import datetime
from PhyNetPy.MPSugar import MP_SUGAR

# Your corrected taxon mapping
taxon_map : dict[str, list[str]] = {
   'Helia': ['Helia'],
    'RS244': ['RS244'],
    'RS247': ['RS247'],
    'RS265': ['RS265'],
    'RS306': ['RS306'],
    'RS416': ['RS416'],
    'RS418': ['RS418'],
    'RS421': ['RS421'],
    'RSC01': ['RSC01', 'RSC01_1', 'RSC01_2', 'RSC01_3', 'RSC01_4', 'RSC01_5', 'RSC01_6', 'RSC01_7'],
    'RSC02': ['RSC02', 'RSC02_1', 'RSC02_2', 'RSC02_3', 'RSC02_4', 'RSC02_5', 'RSC02_6'],
    'RSC03': ['RSC03', 'RSC03_1', 'RSC03_2', 'RSC03_3', 'RSC03_4'],
    'RSC04': ['RSC04', 'RSC04_1', 'RSC04_2', 'RSC04_3', 'RSC04_4', 'RSC04_5'],
    'RSC05': ['RSC05', 'RSC05_1', 'RSC05_2', 'RSC05_3', 'RSC05_4', 'RSC05_5'],
    'RSC06': ['RSC06', 'RSC06_1', 'RSC06_2', 'RSC06_3', 'RSC06_4', 'RSC06_5', 'RSC06_6', 'RSC06_7', 'RSC06_8'],
    'RSC07': ['RSC07', 'RSC07_1', 'RSC07_2', 'RSC07_3'],
    'Senec': ['Senec']
}

num_hill_climbing_chains = 1

for dummy in range(num_hill_climbing_chains):
    
    # Change pathname to trimmed nexus file if doing that version, and ensure taxon map parameter is set to the correct variable 
    # Use an absolute path if J_(un)trimmed.nex is not in current running directory, else just use "J_(un)trimmed.nex"
    output_networks : dict[tuple[DAG, int]] = MP_SUGAR("/groups/itay_mayrose/tomulanovski/gene2net/papers/Ren_2024/gene_trees/mpsugar_test.nex", taxon_map, iter_ct = 5, seed = random.randint(0, 1000))

    #STEP 3: Analyze output
    for net, score in output_networks.items():
        print(net.to_newick())
        print(f"This network scored {score}!")
    
        # Do any other analysis you'd like to do here with the net object:
        #
        # INSERT STUFF HERE. 
        #
