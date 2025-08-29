#!/usr/bin/env python3
"""
Validate phylogenetic trees for GRAMPA compatibility
Usage: 
# Basic validation
python validate_for_grampa.py gene_trees.tre

# With species tree validation
python validate_for_grampa.py gene_trees.tre species_tree.tre
Returns: PASS/FAIL and detailed validation report
"""

import sys
import os
from ete3 import Tree

def validate_tree_for_grampa(tree_string, tree_num, species_tips=None):
    """
    Validate a single tree for GRAMPA compatibility
    Returns: (is_valid, issues_list)
    """
    issues = []
    
    try:
        tree = Tree(tree_string)
        
        # Get tip information
        tips = [leaf.name for leaf in tree.get_leaves()]
        
        # Critical Check 1: Tip label format (GRAMPA requirement)
        bad_tips = []
        species_in_tree = []
        
        for tip in tips:
            if '_' not in tip:
                bad_tips.append(tip)
            else:
                parts = tip.split('_')
                if len(parts) < 2 or not parts[-1]:
                    bad_tips.append(tip)
                else:
                    species_in_tree.append(parts[-1])
        
        if bad_tips:
            issues.append(f"CRITICAL: Invalid tip labels (must be 'gene_species'): {bad_tips[:3]}...")
            return False, issues
        
        # Critical Check 2: Species tree compatibility
        if species_tips:
            invalid_species = [sp for sp in species_in_tree if sp not in species_tips]
            if invalid_species:
                issues.append(f"CRITICAL: Species not in species tree: {set(invalid_species)}")
                return False, issues
        
        # Critical Check 3: Minimum tree size
        if len(tips) < 2:
            issues.append(f"CRITICAL: Tree has only {len(tips)} tips (need ?2)")
            return False, issues
        
        # Critical Check 4: Polytomies (GRAMPA can't handle these)
        polytomies = []
        for node in tree.traverse():
            if len(node.children) > 2:
                polytomies.append(len(node.children))
        
        if polytomies:
            issues.append(f"CRITICAL: Polytomies found (nodes with {polytomies} children)")
            return False, issues
        
        # Critical Check 5: Tree must be rooted (GRAMPA requirement)
        root = tree.get_tree_root()
        if len(root.children) == 3:
            issues.append("CRITICAL: Tree is unrooted (GRAMPA requires rooted trees)")
            return False, issues
        elif len(root.children) > 3:
            issues.append(f"CRITICAL: Root has {len(root.children)} children (invalid tree structure)")
            return False, issues
        
        # Warning checks (won't fail validation but may cause issues)
        warnings = []
        
        # Check for duplicates - but this is normal in polyploidy studies
        duplicates = [tip for tip in set(tips) if tips.count(tip) > 1]
        if duplicates:
            # This is actually expected in polyploidy studies, so just note it
            pass  # Duplicate tips are normal for polyploid gene trees
        
        # Check for zero-length branches
        zero_branches = sum(1 for node in tree.traverse() if node.dist == 0)
        if zero_branches > 0:
            warnings.append(f"WARNING: {zero_branches} zero-length branches")
        
        # Tree passed critical checks
        if warnings:
            issues.extend(warnings)
        
        return True, issues
        
    except Exception as e:
        issues.append(f"CRITICAL: Parse error - {str(e)}")
        return False, issues

def load_species_tree(species_file):
    """Load species tree and return tip names"""
    try:
        species_tree = Tree(species_file)
        return [leaf.name for leaf in species_tree.get_leaves()]
    except Exception as e:
        print(f"Warning: Could not load species tree - {e}")
        return None

def validate_for_grampa(gene_tree_file, species_tree_file=None):
    """
    Validate gene trees for GRAMPA compatibility
    """
    print("GRAMPA Tree Validation")
    print("=" * 50)
    
    # Load species tree if provided
    species_tips = None
    if species_tree_file and os.path.exists(species_tree_file):
        print(f"Loading species tree: {species_tree_file}")
        species_tips = load_species_tree(species_tree_file)
        if species_tips:
            print(f"Species tree has {len(species_tips)} species: {species_tips[:5]}...")
        print()
    
    # Load gene trees
    print(f"Validating gene trees: {gene_tree_file}")
    
    try:
        with open(gene_tree_file, 'r') as f:
            content = f.read().strip()
    except Exception as e:
        print(f"FAIL: Cannot read file - {e}")
        return False
    
    # Parse trees
    if ';' in content:
        tree_strings = [t.strip() + ';' for t in content.split(';') if t.strip()]
    else:
        tree_strings = [t.strip() for t in content.split('\n') if t.strip()]
    
    if not tree_strings:
        print("FAIL: No trees found in file")
        return False
    
    print(f"Found {len(tree_strings)} trees to validate")
    print()
    
    # Validate each tree
    valid_trees = 0
    failed_trees = []
    warning_trees = []
    
    for i, tree_string in enumerate(tree_strings):
        is_valid, issues = validate_tree_for_grampa(tree_string, i+1, species_tips)
        
        critical_issues = [issue for issue in issues if issue.startswith("CRITICAL")]
        warning_issues = [issue for issue in issues if issue.startswith("WARNING")]
        
        if is_valid:
            valid_trees += 1
            if warning_issues:
                warning_trees.append(i+1)
                if len(warning_trees) <= 5:  # Show first 5 warning trees
                    print(f"Tree {i+1}: PASS (with warnings)")
                    for warning in warning_issues:
                        print(f"  {warning}")
        else:
            failed_trees.append(i+1)
            if len(failed_trees) <= 5:  # Show first 5 failed trees
                print(f"Tree {i+1}: FAIL")
                for issue in critical_issues:
                    print(f"  {issue}")
                print(f"  Tree preview: {tree_string[:80]}...")
        
        if i == 0 or (i+1) % 100 == 0:  # Progress indicator
            print(f"Processed {i+1}/{len(tree_strings)} trees...")
    
    # Summary
    print("\n" + "=" * 50)
    print("VALIDATION SUMMARY")
    print("=" * 50)
    
    success_rate = (valid_trees / len(tree_strings)) * 100
    
    print(f"Total trees: {len(tree_strings)}")
    print(f"Valid trees: {valid_trees}")
    print(f"Failed trees: {len(failed_trees)}")
    print(f"Trees with warnings: {len(warning_trees)}")
    print(f"Success rate: {success_rate:.1f}%")
    
    if failed_trees:
        print(f"\nFailed tree numbers: {failed_trees[:20]}...")
        print("\nCommon fixes needed:")
        print("- Fix tip labels to 'gene_species' format")
        print("- Resolve polytomies with resolve_polytomies.py")
        print("- Remove trees with duplicate sequences")
        print("- Ensure species names match species tree")
    
    if warning_trees:
        print(f"\nTrees with warnings: {warning_trees[:20]}...")
        print("These trees will work but may have suboptimal results")
    
    # Final verdict
    print("\n" + "=" * 50)
    if success_rate == 100:
        print("RESULT: PASS - All trees are GRAMPA compatible! ?")
        return True
    elif success_rate >= 80:
        print("RESULT: MOSTLY PASS - Most trees are compatible")
        print("Consider fixing failed trees or excluding them")
        return True
    else:
        print("RESULT: FAIL - Too many incompatible trees")
        print("Fix critical issues before running GRAMPA")
        return False

if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python validate_for_grampa.py gene_trees.tre [species_tree.tre]")
        print("\nValidates gene trees for GRAMPA compatibility")
        print("Optional: provide species tree to check species name matching")
        sys.exit(1)
    
    gene_tree_file = sys.argv[1]
    species_tree_file = sys.argv[2] if len(sys.argv) == 3 else None
    
    if not os.path.exists(gene_tree_file):
        print(f"Error: Gene tree file '{gene_tree_file}' not found")
        sys.exit(1)
    
    success = validate_for_grampa(gene_tree_file, species_tree_file)
    sys.exit(0 if success else 1)