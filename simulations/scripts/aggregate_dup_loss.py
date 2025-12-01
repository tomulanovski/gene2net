#!/usr/bin/env python3
"""
Aggregate duplication and loss statistics across all networks
(Handles multiple replicates per network)
"""

import sys
from pathlib import Path

def aggregate_summaries(base_dir, networks, data_subdir, analysis_dir, output_filename):
    """Aggregate dup/loss summaries across all networks and their replicates"""
    
    # Global stats (1 value per network, which is the average of its replicates)
    network_final_dups = []
    network_final_losses = []
    
    # To store details for the report
    network_report_stats = []
    
    for network in networks:
        network_path = Path(base_dir) / network / data_subdir
        
        # Find all replicate directories (replicate_1, replicate_2, etc.)
        replicate_dirs = sorted(list(network_path.glob("replicate_*")))
        
        if not replicate_dirs:
            print(f"WARNING: No replicate directories found for {network} in {network_path}")
            continue

        # Store stats for this specific network's replicates
        current_net_dups = []
        current_net_losses = []
        
        for rep_dir in replicate_dirs:
            summary_path = rep_dir / "dup_loss_summary.txt"
            
            if not summary_path.exists():
                # It's possible a specific replicate failed, just skip it
                continue
            
            # Parse the summary file
            try:
                avg_dup = None
                avg_loss = None
                
                with open(summary_path, 'r') as f:
                    for line in f:
                        if line.startswith("Average duplications:"):
                            avg_dup = float(line.split(':')[1].strip())
                        elif line.startswith("Average losses:"):
                            avg_loss = float(line.split(':')[1].strip())
                
                if avg_dup is not None and avg_loss is not None:
                    current_net_dups.append(avg_dup)
                    current_net_losses.append(avg_loss)
                    
            except Exception as e:
                print(f"ERROR parsing {summary_path}: {e}")
                continue

        # If we found valid data for this network
        if current_net_dups:
            # Calculate average across replicates for this network
            net_avg_dup = sum(current_net_dups) / len(current_net_dups)
            net_avg_loss = sum(current_net_losses) / len(current_net_losses)
            
            # Add to global lists
            network_final_dups.append(net_avg_dup)
            network_final_losses.append(net_avg_loss)
            
            # Store for writing to file
            network_report_stats.append((network, net_avg_dup, net_avg_loss, len(current_net_dups)))
            
            print(f"? {network}: dup={net_avg_dup:.4f}, loss={net_avg_loss:.4f} (from {len(current_net_dups)} reps)")
        else:
            print(f"WARNING: No valid summary files found for {network}")

    if not network_final_dups:
        print("ERROR: No valid data found for any network!")
        return False
    
    # Calculate overall averages (Average of averages)
    overall_avg_dup = sum(network_final_dups) / len(network_final_dups)
    overall_avg_loss = sum(network_final_losses) / len(network_final_losses)
    
    # Write aggregate summary
    Path(analysis_dir).mkdir(parents=True, exist_ok=True)
    output_path = Path(analysis_dir) / output_filename
    
    with open(output_path, 'w') as f:
        f.write(f"Aggregate Duplication and Loss Summary Across All Networks\n")
        f.write(f"{'='*60}\n\n")
        f.write(f"Dataset Config: {data_subdir}\n")
        f.write(f"Number of networks analyzed: {len(network_report_stats)}\n\n")
        f.write(f"Overall average duplications: {overall_avg_dup:.4f}\n")
        f.write(f"Overall average losses: {overall_avg_loss:.4f}\n\n")
        f.write(f"Note: Values below are averages across available replicates.\n")
        f.write(f"{'='*60}\n")
        f.write(f"Per-Network Statistics:\n")
        f.write(f"{'='*60}\n\n")
        
        for network, avg_dup, avg_loss, num_reps in network_report_stats:
            f.write(f"{network} ({num_reps} reps):\n")
            f.write(f"  Average duplications: {avg_dup:.4f}\n")
            f.write(f"  Average losses: {avg_loss:.4f}\n\n")
    
    print(f"\n{'='*60}")
    print(f"Aggregate summary written to: {output_path}")
    print(f"{'='*60}")
    print(f"Networks analyzed: {len(network_report_stats)}")
    print(f"Overall average duplications: {overall_avg_dup:.4f}")
    print(f"Overall average losses: {overall_avg_loss:.4f}")
    
    return True

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python aggregate_dup_loss.py <data_subdir> <output_filename>")
        sys.exit(1)
    
    data_subdir = sys.argv[1]
    output_filename = sys.argv[2]
    
    # Network list
    networks = [
        "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
        "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019", 
        "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011", 
        "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014", 
        "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014"
    ]
    
    base_dir = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
    analysis_dir = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis"
    
    success = aggregate_summaries(base_dir, networks, data_subdir, analysis_dir, output_filename)
    
    if not success:
        sys.exit(1)