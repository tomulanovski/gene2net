import pandas as pd
import re
from pathlib import Path

# Paths
txt_dir = Path("/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis/")
csv_file = Path("/groups/itay_mayrose/tomulanovski/gene2net/dup_loss_summary.csv")
output_file = Path("/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis/merged_dup_loss_summary.csv")

def parse_txt_file(file_path):
    """Parse a txt file and extract network statistics."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    networks = {}
    
    # Split by lines and look for the pattern: NetworkName:\n  Average duplications: X\n  Average losses: Y
    lines = content.split('\n')
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Check if this line ends with a colon (potential network name)
        if line and line.endswith(':') and not line.startswith('='):
            network_name = line[:-1]  # Remove the colon
            
            # Skip if it's a header line
            if network_name in ['Aggregate Duplication and Loss Summary Across All Networks', 
                                'Per-Network Statistics', 'Dataset']:
                i += 1
                continue
            
            # Look ahead for "Average duplications" and "Average losses"
            avg_dup = None
            avg_loss = None
            
            for j in range(i+1, min(i+5, len(lines))):
                if 'Average duplications:' in lines[j]:
                    match = re.search(r'Average duplications:\s*([\d.]+)', lines[j])
                    if match:
                        avg_dup = float(match.group(1))
                elif 'Average losses:' in lines[j]:
                    match = re.search(r'Average losses:\s*([\d.]+)', lines[j])
                    if match:
                        avg_loss = float(match.group(1))
            
            # If we found both values, add to networks
            if avg_dup is not None and avg_loss is not None:
                networks[network_name] = {'dup': avg_dup, 'loss': avg_loss}
        
        i += 1
    
    return networks

def extract_file_suffix(filename):
    """Extract the suffix from dup_loss_summary_*.txt files."""
    match = re.search(r'dup_loss_summary_(.+)\.txt', filename)
    if match:
        return match.group(1)
    return filename.replace('.txt', '')

# Step 1: Read the CSV file and build the backbone
print(f"Reading CSV file: {csv_file}")
df_csv = pd.read_csv(csv_file)
print(f"CSV shape: {df_csv.shape}")
print(f"CSV columns: {df_csv.columns.tolist()}")

# Filter out MEAN and STD rows from CSV
df_csv_filtered = df_csv[~df_csv['dataset'].isin(['MEAN', 'STD'])].copy()
print(f"CSV shape after filtering MEAN/STD: {df_csv_filtered.shape}")

# Create the backbone dataframe with network names from CSV
df_final = pd.DataFrame()
df_final['network'] = df_csv_filtered['dataset'].values
df_final['real_dup'] = df_csv_filtered['avg_duplications'].values
df_final['real_loss'] = df_csv_filtered['avg_losses'].values

print(f"\nBackbone created with {len(df_final)} networks from CSV")
print(f"Networks: {df_final['network'].tolist()}")

# Step 2: Process all txt files
txt_files = sorted(txt_dir.glob("dup_loss_summary_*.txt"))
print(f"\nFound {len(txt_files)} txt files to process")

# Collect all networks from all txt files first
all_txt_networks = set(df_final['network'].tolist())

for txt_file in txt_files:
    print(f"\nProcessing: {txt_file.name}")
    
    suffix = extract_file_suffix(txt_file.name)
    print(f"  Suffix: {suffix}")
    
    networks = parse_txt_file(txt_file)
    print(f"  Found {len(networks)} networks")
    print(f"  Networks: {list(networks.keys())}")
    
    # Add any new networks to our set
    all_txt_networks.update(networks.keys())
    
    # Create column names
    if suffix.endswith('_dup'):
        dup_col = suffix
        loss_col = suffix.replace('_dup', '_loss')
    else:
        dup_col = f'{suffix}_dup'
        loss_col = f'{suffix}_loss'
    
    # Initialize columns with NaN
    df_final[dup_col] = None
    df_final[loss_col] = None
    
    # Fill in values for networks that exist in this txt file
    for network_name, stats in networks.items():
        # Find the row for this network (or add if it doesn't exist)
        if network_name in df_final['network'].values:
            idx = df_final[df_final['network'] == network_name].index[0]
            df_final.at[idx, dup_col] = stats['dup']
            df_final.at[idx, loss_col] = stats['loss']
        else:
            # Network doesn't exist in dataframe yet, add a new row
            new_row = {'network': network_name, dup_col: stats['dup'], loss_col: stats['loss']}
            df_final = pd.concat([df_final, pd.DataFrame([new_row])], ignore_index=True)
            print(f"  Added new network: {network_name}")

# Step 3: Sort and clean up
df_final = df_final.sort_values('network').reset_index(drop=True)

# Convert numeric columns to float
for col in df_final.columns:
    if col != 'network':
        df_final[col] = pd.to_numeric(df_final[col], errors='coerce')

# Step 4: Calculate MEAN and STD
mean_row = {'network': 'MEAN'}
std_row = {'network': 'STD'}

for col in df_final.columns:
    if col != 'network':
        mean_row[col] = df_final[col].mean()
        std_row[col] = df_final[col].std()

df_mean = pd.DataFrame([mean_row])
df_std = pd.DataFrame([std_row])
df_final = pd.concat([df_final, df_mean, df_std], ignore_index=True)

# Step 5: Save and report
print(f"\n{'='*60}")
print(f"Final DataFrame shape: {df_final.shape}")
print(f"Networks (rows): {len(df_final) - 2} (+ MEAN and STD rows)")
print(f"Columns: {len(df_final.columns)}")
print(f"\nColumn names:")
for col in df_final.columns:
    print(f"  - {col}")

df_final.to_csv(output_file, index=False)
print(f"\n{'='*60}")
print(f"Merged CSV saved to: {output_file}")
print(f"\nFirst few rows:")
print(df_final.head(10))
print(f"\nLast few rows (including MEAN/STD):")
print(df_final.tail(5))
print(f"\nSummary of missing values:")
print(df_final.isnull().sum())