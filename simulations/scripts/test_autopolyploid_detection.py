#!/usr/bin/env python3
"""
Test script to verify autopolyploidization event detection on Ding_2023.tre

Key distinction:
- Autopolyploidization EVENT: The single WGD that duplicated the Rch clade
- Autopolyploid SPECIES: The 7 species that appear twice (Rch, Pla, Jmi, CiP, Jre, Jma, Aro)

Expected result for Ding_2023.tre:
- 7 polyploid species (species appearing >1 time)
- 0 reticulations (H_Strict = 0, no hybridization)
- 1 autopolyploidization event (the Rch clade duplication)
- Total_WGD = 1
"""

import sys
from pathlib import Path

# Add scripts directory to path
scripts_dir = Path(__file__).parent
sys.path.insert(0, str(scripts_dir))

from run_reticulation_stats import analyze_mul_trees

# Test on the networks directory
networks_dir = scripts_dir.parent / "networks"

print("Testing autopolyploidization event detection on Ding_2023.tre...")
print("="*80)

df = analyze_mul_trees(str(networks_dir))

# Filter to Ding_2023.tre
ding = df[df['Filename'] == 'Ding_2023.tre']

if not ding.empty:
    row = ding.iloc[0]
    print(f"\nResults for Ding_2023.tre:")
    print(f"  Num_Species: {row['Num_Species']}")
    print(f"  Num_Polyploids (species with >1 copy): {row['Num_Polyploids']}")
    print(f"  Max_Copies: {row['Max_Copies']}")
    print(f"  H_Strict (reticulations/hybridization): {row['H_Strict']}")
    print(f"  Num_Autopolyploidization_Events: {row['Num_Autopolyploidization_Events']}")
    print(f"  Total_WGD: {row['Total_WGD']}")
    print(f"  Polyploid species: {row['Polyploid_Names']}")

    print("\n" + "="*80)
    print("Expected:")
    print("  Num_Polyploids: 7")
    print("    (Aro, CiP, Jma, Jmi, Jre, Pla, Rch all appear twice)")
    print("  H_Strict: 0")
    print("    (no reticulations/hybridization)")
    print("  Num_Autopolyploidization_Events: 1")
    print("    (ONE event: the Rch clade duplication)")
    print("  Total_WGD: 1")
    print("    (1 autopolyploidization + 0 reticulations)")

    print("\n" + "="*80)
    print("Explanation:")
    print("  The tree has two identical Rch subtrees (siblings).")
    print("  This represents ONE autopolyploidization event that duplicated")
    print("  the entire Rch clade, which contains 7 species.")
    print("  Each of those 7 species now appears twice (= 7 autopolyploid species).")
    print("  But we count it as 1 EVENT, not 7 events!")

    print("\n" + "="*80)
    if (row['Num_Autopolyploidization_Events'] == 1 and
        row['H_Strict'] == 0 and
        row['Total_WGD'] == 1 and
        row['Num_Polyploids'] == 7):
        print("✓ TEST PASSED")
    else:
        print("✗ TEST FAILED - values don't match expected")
        print(f"  Got: auto_events={row['Num_Autopolyploidization_Events']}, H_Strict={row['H_Strict']}, Total_WGD={row['Total_WGD']}, Polyploids={row['Num_Polyploids']}")
        print(f"  Expected: auto_events=1, H_Strict=0, Total_WGD=1, Polyploids=7")
else:
    print("ERROR: Ding_2023.tre not found in results")
