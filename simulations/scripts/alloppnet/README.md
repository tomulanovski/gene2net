# AlloppNET Scripts Directory

This directory contains scripts for running AlloppNET phylogenetic network inference on simulation data.

## Required Setup

**IMPORTANT:** Before running AlloppNET, you must copy the AlloppDT R scripts from the cluster:

```bash
# On the cluster:
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/alloppnet/

# Copy AlloppDT scripts
cp "/groups/itay_mayrose/tomulanovski/gene2net/jones_examples/2013-05-15-manual (2)/AlloppDT_5beastxml_toplevel.r" .
cp "/groups/itay_mayrose/tomulanovski/gene2net/jones_examples/2013-05-15-manual (2)/AlloppDT_6beastxml_bits.r" .
```

After copying, you should have:
- `AlloppDT_5beastxml_toplevel.r`
- `AlloppDT_6beastxml_bits.r`
- `generate_beast_xml.r` (wrapper script)
- `prepare_alloppnet_input.py` (preprocessing)
- `remove_copy_numbers.py` (post-processing)

## Scripts

### prepare_alloppnet_input.py
Converts PHY alignments to NEXUS, generates ploidy_level.json and taxa_table.txt

### generate_beast_xml.r
Generates BEAST XML using AlloppDT scripts

### remove_copy_numbers.py
Removes _0, _1 suffixes from consensus trees

## Usage

See `simulations/md_files/ALLOPPNET_GUIDE.md` for complete usage instructions.
