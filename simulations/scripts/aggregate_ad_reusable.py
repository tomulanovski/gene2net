# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Aggregate AD (Average Distance) results from SimPhy simulations.
Usage: python aggregate_ad_simphy.py [config_name]
Example: python aggregate_ad_simphy.py conf_ils_low_10M
"""
import os
import re
import csv
import numpy as np
import sys

def extract_ad_from_file(filepath):
    """Extract AD value from rf_distance_results.txt file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            match = re.search(r'Average Distance \(AD\):\s+([\d.]+)\s+\(([\d.]+)%\)', content)
            if match:
                return float(match.group(1)), float(match.group(2))
    except Exception as e:
        print(f"  Error reading {filepath}: {e}")
    return None, None

def main():
    # Get configuration from command line or use default
    if len(sys.argv) > 1:
        SIMPHY_CONFIG = sys.argv[1]
    else:
        SIMPHY_CONFIG = "conf_ils_low_10M"
    
    BASE_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
    OUTPUT_CSV = f"/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis/aggregate_AD_{SIMPHY_CONFIG}.csv"
    
    NUM_REPLICATES = 5
    
    networks = [
        "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
        "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019", "Popp_2005", "Wu_2015",
        "Liu_2023", "Ren_2024", "Marcussen_2011", "Marcussen_2012", "Sessa_2012b", "Zhao_2021",
        "Hori_2014", "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014"
    ]
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    
    print("Aggregating results...")
    print("=" * 60)
    print(f"Configuration: {SIMPHY_CONFIG}")
    print("=" * 60)
    
    # Collect results per replicate and per network
    all_results = []  # Every single replicate
    network_summaries = []  # Aggregated per network
    not_found = []
    
    for network in networks:
        network_ads = []
        replicates_found = 0
        
        for rep in range(1, NUM_REPLICATES + 1):
            rf_file = os.path.join(BASE_DIR, network, "data", SIMPHY_CONFIG, 
                                   f"replicate_{rep}", "rf_distance_results.txt")
            
            if os.path.exists(rf_file):
                raw_ad, percent_ad = extract_ad_from_file(rf_file)
                if percent_ad is not None:
                    network_ads.append(percent_ad)
                    all_results.append({
                        'network': network,
                        'replicate': rep,
                        'AD': percent_ad
                    })
                    replicates_found += 1
        
        # Calculate network summary
        if network_ads:
            mean_ad = np.mean(network_ads)
            std_ad = np.std(network_ads, ddof=1) if len(network_ads) > 1 else 0
            network_summaries.append({
                'network': network,
                'replicates': replicates_found,
                'mean_AD': mean_ad,
                'std_AD': std_ad
            })
            print(f"? {network:30s} {replicates_found}/{NUM_REPLICATES} reps, Mean AD = {mean_ad:.2f}% +/- {std_ad:.2f}%")
        else:
            not_found.append(network)
            print(f"? {network:30s} No results found")
    
    # Save to CSV
    if network_summaries:
        # Calculate overall statistics
        all_mean_ads = [s['mean_AD'] for s in network_summaries]
        overall_mean = np.mean(all_mean_ads)
        overall_std = np.std(all_mean_ads, ddof=1)
        
        with open(OUTPUT_CSV, 'w', newline='') as csvfile:
            fieldnames = ['network', 'replicates', 'mean_AD', 'std_AD']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            # Write network summaries
            for summary in network_summaries:
                writer.writerow(summary)
            
            # Add separator
            writer.writerow({'network': '', 'replicates': '', 'mean_AD': '', 'std_AD': ''})
            
            # Write overall statistics
            writer.writerow({
                'network': 'OVERALL_MEAN',
                'replicates': len(network_summaries),
                'mean_AD': f'{overall_mean:.2f}',
                'std_AD': f'{overall_std:.2f}'
            })
        
        print("\n" + "=" * 60)
        print(f"SUCCESS: Results saved to:")
        print(f"  {OUTPUT_CSV}")
        print(f"\nNetworks processed: {len(network_summaries)}/{len(networks)}")
        print(f"Total replicates: {len(all_results)}")
        print(f"\nOverall Statistics:")
        print(f"  Mean AD: {overall_mean:.2f}% +/- {overall_std:.2f}%")
        print(f"  Range: {min(all_mean_ads):.2f}% - {max(all_mean_ads):.2f}%")
        
        if not_found:
            print(f"\nNetworks with missing results ({len(not_found)}):")
            for net in not_found:
                print(f"  - {net}")
        print("=" * 60)
    else:
        print("\nERROR: No results found!")
        exit(1)

if __name__ == "__main__":
    main()