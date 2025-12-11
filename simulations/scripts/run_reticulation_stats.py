import os
import glob
import pandas as pd
# Ensure reticulate_tree.py is in the same directory
from reticulate_tree import ReticulateTree 

# Path to your files
DIR_PATH = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks/"

def analyze_mul_trees(directory, relaxed_threshold=0.2):
    files = glob.glob(os.path.join(directory, "*"))
    results = []

    print(f"Found {len(files)} files. Analyzing...")

    for filepath in files:
        # Filter for valid tree files only
        if os.path.isdir(filepath) or not filepath.endswith(('.tre', '.nwk', '.tree')):
            continue
            
        filename = os.path.basename(filepath)
        
        try:
            with open(filepath, 'r') as f:
                content = f.read().strip()
            
            if not content or not content.endswith(';'):
                continue

            # --- 1. GET BASIC STATS (From Strict Instance) ---
            rt_strict = ReticulateTree(content, is_multree=True)
            stats_strict = rt_strict.measure(printout=False)
            
            # --- 2. CALCULATE LEAF STATS ---
            leaf_counts = stats_strict['leaf_counts'] # Dict like {'SpA': 2, 'SpB': 1}
            
            num_species = len(leaf_counts)
            
            # Count how many species have > 1 copy
            polyploid_species = [name for name, count in leaf_counts.items() if count > 1]
            num_polyploids = len(polyploid_species)
            
            # Find the maximum number of copies any single species has
            max_copies = max(leaf_counts.values()) if leaf_counts else 0
            
            # --- 3. GET RELAXED STATS ---
            rt_relaxed = ReticulateTree(content, is_multree=True, threshold=relaxed_threshold, normalize=True)
            stats_relaxed = rt_relaxed.measure(printout=False)

            # --- 4. COMPILE ROW ---
            h_strict = stats_strict['reticulation_count']
            h_relaxed = stats_relaxed['reticulation_count']

            row = {
                'Filename': filename,
                'Num_Species': num_species,
                'Num_Polyploids': num_polyploids,
                'Max_Copies': max_copies,            # <--- NEW COLUMN
                'H_Strict': h_strict,
                'H_Relaxed': h_relaxed,
                'H_Diff': h_relaxed - h_strict,
                'Polyploid_Names': ", ".join(sorted(polyploid_species))
            }
            results.append(row)

        except Exception as e:
            print(f"Error on {filename}: {e}")

    return pd.DataFrame(results)

if __name__ == "__main__":
    df = analyze_mul_trees(DIR_PATH)

    if not df.empty:
        # Save to CSV
        out_file = os.path.join(DIR_PATH, "mul_tree_final_stats.csv")
        df.to_csv(out_file, index=False)
        
        print("\n--- Summary Statistics ---")
        # Included Max_Copies in the printed summary
        cols = ['Filename', 'Num_Species', 'Num_Polyploids', 'Max_Copies', 'H_Strict', 'H_Relaxed', 'H_Diff']
        print(df[cols].to_string())
        print(f"\nStats saved to: {out_file}")
    else:
        print("No valid data found.")