import pandas as pd
import os
from collections import defaultdict
import re
import unicodedata

# Define base paths
BASE_PATH = "/groups/itay_mayrose/tomulanovski/gene2net"
CSV_PATH = os.path.join(BASE_PATH, "data", "Genes&Reads.csv")
PAPERS_PATH = os.path.join(BASE_PATH, "papers")

# Method names
METHODS = ['grampa', 'padre', 'polyphest', 'mpallop', 'mpl', 'allopnet']

# Read the CSV file
df = pd.read_csv(CSV_PATH)

# Dictionary to track authors/years for lettering
author_year_count = defaultdict(int)

def normalize_name(name):
    """Remove accents and special characters from name"""
    # Normalize unicode characters
    normalized = unicodedata.normalize('NFKD', name)
    # Remove non-ASCII characters but keep basic punctuation
    normalized = ''.join(c for c in normalized if not unicodedata.combining(c))
    return normalized

def extract_author_year(paper_text):
    """Extract author and year from paper text."""
    # Handle papers that already have a/b suffix
    match = re.search(r'([A-Za-zÀ-ÿ\-]+)\s+et\s+al\s+(\d{4}[a-z]?)', paper_text)
    if match:
        author = match.group(1)
        year = match.group(2)
        # Normalize the author name
        author = normalize_name(author)
        # If year already has a letter, strip it for our counter
        base_year = year[:4]
        existing_letter = year[4:] if len(year) > 4 else ''
        return author, base_year, existing_letter
    return None, None, None

def create_paper_directories(row):
    # Extract author, year, and any existing letter
    author, year, existing_letter = extract_author_year(row['Ref'])
    if not author or not year:
        print(f"Warning: Could not extract author/year from: {row['Ref']}")
        return None
    
    # Create unique directory name
    key = f"{author}_{year}"
    
    # If there's an existing letter, use it
    if existing_letter:
        base_dir = f"{author}_{year}{existing_letter}"
    else:
        # Add letter if duplicate author/year
        base_dir = f"{author}_{year}"
        if author_year_count[key] > 0:
            base_dir += chr(97 + author_year_count[key])
    author_year_count[key] += 1
    
    # Create main paper directory and subdirectories
    paper_dir = os.path.join(PAPERS_PATH, base_dir)
    subdirs = ['genes', 'gene_trees', 'networks']
    
    # Create all directories
    os.makedirs(paper_dir, exist_ok=True)
    for subdir in subdirs:
        os.makedirs(os.path.join(paper_dir, subdir), exist_ok=True)
        
        # Create method subdirectories in networks
        if subdir == 'networks':
            for method in METHODS:
                method_dir = os.path.join(paper_dir, subdir, method)
                os.makedirs(method_dir, exist_ok=True)
    
    print(f"Created directory structure for: {base_dir}")
    return base_dir

# Create directories for each paper
print("Starting directory creation...")
for index, row in df.iterrows():
    dir_name = create_paper_directories(row)

print("Directory creation completed!")