#!/usr/bin/env python3
"""
Export all tables from a SQLite database to CSV files.
Usage: python "/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/extract_sql.py" path_to_database.db [output_dir(optional. default is database directory)]
"""

import sqlite3
import csv
import sys
import os
from pathlib import Path

def export_sqlite_to_csv(db_path, output_dir=None):
    """Export all tables from SQLite database to CSV files."""
    
    # Set output directory - default to same directory as the database
    if output_dir is None:
        output_dir = Path(db_path).parent
        if output_dir == Path('.'):  # Handle relative paths
            output_dir = Path.cwd()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Connect to database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Get all table names
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    
    print(f"Found {len(tables)} tables in {db_path}")
    print(f"Exporting to: {output_dir}/\n")
    
    # Export each table
    for table_name, in tables:
        print(f"Exporting table: {table_name}... ", end="")
        
        # Query all data from table
        cursor.execute(f"SELECT * FROM {table_name}")
        rows = cursor.fetchall()
        
        # Get column names
        column_names = [description[0] for description in cursor.description]
        
        # Write to CSV
        csv_path = os.path.join(output_dir, f"{table_name}.csv")
        with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(column_names)  # Write header
            writer.writerows(rows)  # Write data
        
        print(f"? ({len(rows)} rows)")
    
    conn.close()
    print(f"\nExport complete! Files saved to: {output_dir}/")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python export_sqlite.py database.db [output_dir]")
        sys.exit(1)
    
    db_path = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.exists(db_path):
        print(f"Error: Database file '{db_path}' not found")
        sys.exit(1)
    
    export_sqlite_to_csv(db_path, output_dir)