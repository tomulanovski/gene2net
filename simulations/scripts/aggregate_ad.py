#!/usr/bin/env python3
"""
Aggregate AD (Average Distance) results from all networks into a summary CSV file.
Usage: python aggregate_ad_results.py
"""
import os
import re
import csv
import numpy as np

def extract_ad_from_file(filepath):
    """Extract AD value from rf_distance_results.txt file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            match = re.search(r'Average Distance \(AD\):\s+([\d.]+)\s+\(([\d.]+)%\)', content)
            if match:
                return float(match.group(1)), float(match.group(2))
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return None, None

def main():
    BASE_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
    OUTPUT_CSV = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis/aggregate_AD_low_ILS_low_dup.csv"
    
    networks = [
        "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
        "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019", "Popp_2005", "Wu_2015",
        "Liu_2023", "Ren_2024", "Marcussen_2011", "Marcussen_2012", "Sessa_2012b", "Zhao_2021",
        "Hori_2014", "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014"
    ]
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    
    # Collect results
    results = []
    not_found = []
    
    print("Aggregating results...")
    print("=" * 60)
    
    for network in networks:
        rf_file = os.path.join(BASE_DIR, network, "data", "ILS_low_dup_low", "rf_distance_results.txt")
        
        if os.path.exists(rf_file):
            raw_ad, percent_ad = extract_ad_from_file(rf_file)
            if raw_ad is not None:
                results.append({
                    'network': network,
                    'AD': percent_ad
                })
                print(f"? {network:30s} AD = {percent_ad:.2f}%")
            else:
                not_found.append(network)
                print(f"? {network:30s} Could not extract AD value")
        else:
            not_found.append(network)
            print(f"? {network:30s} Results file not found")
    
    # Save to CSV
    if results:
        # Calculate summary statistics
        ad_values = [r['AD'] for r in results]
        mean_ad = np.mean(ad_values)
        std_ad = np.std(ad_values, ddof=1)  # Using sample standard deviation
        
        with open(OUTPUT_CSV, 'w', newline='') as csvfile:
            fieldnames = ['network', 'AD']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            # Write individual network results
            for result in results:
                writer.writerow(result)
            
            # Add empty row for separation
            writer.writerow({'network': '', 'AD': ''})
            
            # Write summary statistics
            writer.writerow({'network': 'MEAN', 'AD': f'{mean_ad:.2f}'})
            writer.writerow({'network': 'STD', 'AD': f'{std_ad:.2f}'})
        
        print("\n" + "=" * 60)
        print(f"SUCCESS: Aggregated results saved to:")
        print(f"  {OUTPUT_CSV}")
        print(f"Total networks processed: {len(results)}/{len(networks)}")
        print(f"\nSummary Statistics:")
        print(f"  Mean AD: {mean_ad:.2f}%")
        print(f"  Std AD:  {std_ad:.2f}%")
        
        if not_found:
            print(f"\nNetworks with missing/incomplete results ({len(not_found)}):")
            for net in not_found:
                print(f"  - {net}")
    else:
        print("\nERROR: No results found!")
        exit(1)

if __name__ == "__main__":
    main()