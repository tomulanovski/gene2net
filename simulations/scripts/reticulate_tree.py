import webbrowser
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path
from ete3 import Tree, TreeNode
from pyvis.network import Network
from collections import defaultdict, Counter
import os

# The per-reticulation "malformed inference" warning below floods scoring output
# (thousands of lines on over-reticulate reconstructions) and is not actionable
# during benchmarking. Off by default; set RETICULATE_WARN_MALFORMED=1 to re-enable.
_WARN_MALFORMED_RETIC = os.environ.get("RETICULATE_WARN_MALFORMED") == "1"

# Convension for this class:
# _func() for private methods, may or may not be static
# func() for public methods, therefore we should try to make them static or classmethods [not to modify the instance state]

class ReticulateTree:

    def __init__(self, input_data, is_multree=False, threshold=None, normalize=None,
                 refold=False):
        self.tree = None
        self.tree_str = None
        self.dag = None
        self.retnodes = None
        self.threshold = threshold
        self.normalize = normalize
        self.refold = refold

        input_type = self._determine_input_type(input_data)

        if input_type == 'str:newick':
            self.tree_str = input_data
            self.tree = Tree(input_data, format=1)
            self._init_from_tree(is_multree)

        elif input_type == 'ete3tree':
            self.tree = input_data
            self.tree_str = input_data.write(format=1)
            self._init_from_tree(is_multree)

        elif input_type == 'str:enewick':
            self._init_from_dag(
                self.__class__.enewick_to_dag(input_data)
            )

        elif input_type.endswith('nexml'):
            raise NotImplementedError('NeXML input is not yet supported for ReticulateTree initialization.')
            #dag = self._parse_nexml_to_dag(input_data, input_type.split(':')[0])
            #self._init_from_dag(dag)

        elif input_type == 'nxdag':
            self._init_from_dag(input_data)

        # XOR None check to ensure correct param tracking
        if (self.threshold is None) ^ (self.normalize is None):
            print('Warning: Under the current settings, MUL-tree folding parameters are not used.')
            self.threshold = None
            self.normalize = None
        
        self.retnodes = self.get_reticulation_nodes()

    def _init_from_tree(self, is_multree):
        if is_multree or self.check_duplicated():
            if self.threshold is not None or self.normalize is not None:
                # Polyphest's multree folding defaults
                if self.threshold is None:
                    self.threshold = 0.2
                if self.normalize is None:
                    self.normalize = True
                print("Note: Polyphest's multree folding does NOT guarantee DAG back-translatability to the original (MUL-)tree!")
                self.dag = self.__class__.polyphest_multree_to_dag(self.tree, threshold=self.threshold, normalize=self.normalize)
            else:
                self.dag = self.__class__.multree_to_dag(self.tree, strict=True)
        else:
            self.dag = self.tree_to_dag(self.tree)

    def _init_from_dag(self, dag):
        '''
        Take a network that arrives already folded, such as extended Newick.

        With refold=True the network is unfolded to the MUL-tree it exhibits
        and folded again with Holm, which maps any binary resolution back onto
        the canonical N(T). Huber et al. 2006 prove a minimal binary network
        exhibiting a MUL-tree is a resolution of N(T), so this compares such a
        network at the same resolution level as every method whose output is a
        MUL-tree, instead of penalising it for committing to an event ordering
        that no MUL-tree determines.

        Off by default. get_modified_mu_representation rebuilds a
        ReticulateTree around a suppressed graph, and refolding there would
        silently change the metric.
        '''
        self.tree = self.__class__.dag_to_multree(dag, match_dag=True)
        self.tree_str = self.tree.write(format=1)
        self.dag = (self.__class__.multree_to_dag(self.tree, strict=True)
                    if self.refold else dag)
        
    def check_duplicated(self):
        ''' Check if the tree has duplicated labels (i.e., multiple leaves with the same name). '''
        return any(v > 1 for v in self.get_leaf_counts().values())
  
    @staticmethod
    def _determine_input_type(input_data):
        '''
        Determine the type of input data provided to ReticulateTree.
        Returns a string indicating the type:
        - 'str:newick' for Newick strings
        - 'str:enewick' for enriched Newick strings with #H markers
        - 'str:nexml' for NeXML strings
        - 'ete3tree' for ETE3 Tree objects
        - 'nxdag' for NetworkX directed acyclic graphs (DAGs)
        - 'file:nexml' for file-like objects (assumed to be NeXML)
        '''
        input_type = None
        if isinstance(input_data, str):
            input_str = input_data.strip()
            # Can be Newick, enriched Newick, or NeXML - check format
            if input_str.endswith('nexml>'):
                input_type = 'str:nexml'
                print('Warning: NeXML input is only supported for networks!')
            else:
                if not input_str.endswith(';'):
                    raise ValueError('String format unrecognized. Must end with ";" for Newick / enriched Newick, or "nexml>" for NeXML.')
                else:
                    if '#' in input_str or '//' in input_str:
                        input_type = 'str:enewick'
                    else:
                        input_type = 'str:newick'
        elif isinstance(input_data, Tree):
            input_type = 'ete3tree'
        elif isinstance(input_data, nx.DiGraph):
            input_type = 'nxdag'
        elif hasattr(input_data, 'read'):
            # Assume file-like object
            input_type = 'file:nexml'  # Default to NeXML for file-like objects
            print('Warning: NeXML input is only supported for networks!')
        else:
            raise ValueError(f'Unsupported input type: {type(input_data)}')
        return input_type

    @staticmethod
    def get_dag_root(G):
        ''' Get the root of a directed acyclic graph (DAG). '''
        roots = [n for n, d in G.in_degree() if d == 0]
        if len(roots) != 1:
            raise ValueError('Expected exactly one root.')
        return roots[0]

    @staticmethod
    def compute_depths(G):
        ''' Compute DAG depths to mimic layered layout. '''
        root = ReticulateTree.get_dag_root(G)
        depths = {node: float('-inf') for node in G.nodes}
        depths[root] = 0
        # Topological sort guarantees that all parents are visited before children
        for node in nx.topological_sort(G):
            for child in G.successors(node):
                proposed_depth = depths[node] + 1
                if proposed_depth > depths[child]:
                    depths[child] = proposed_depth
        return depths

    @staticmethod
    def is_reticulation_node(G, node):
        ''' Check if a DAG node is a reticulation node (including autopolyploidy nodes). '''
        return G.in_degree(node) > 1 or G.nodes[node].get('autopolyploidy', False)

    # ───── MODULE: From simple tree to network ─────

    @staticmethod
    def tree_to_dag(t):
        '''
        Convert ete3 singly-labeled tree to networkx DAG (no reticulations)
        '''
        G = nx.DiGraph()
        for node in t.traverse('postorder'):
            nid = id(node)
            G.add_node(
                nid,
                ete=node,         # original ete3 node
                label=node.name,  # display label # node.name if node.is_leaf() else str(nid)
            )
            if node.up:
                if hasattr(node, 'dist'):
                    G.add_edge(id(node.up), nid, length=node.dist)
                else:
                    G.add_edge(id(node.up), nid)
        return G

    # ───── MODULE: From MUL-tree to network ─────

    @staticmethod
    def _compute_height_map(root):
        height_map = defaultdict(list)
        heights = {}  # only used inside this loop to compute parent's height

        for node in root.traverse('postorder'):
            if node.is_leaf():
                h = 0
            else:
                h = max(heights[child] for child in node.children) + 1
            heights[node] = h
            height_map[h].append(node)

        return height_map
    
    @staticmethod
    def _compute_sorted_multisets(root):
        '''
        Returns a dict: {node: tuple(sorted multiset items)} for each node.
        '''
        sorted_multisets = {}
        multisets = {}

        for node in root.traverse('postorder'):
            if node.is_leaf():
                mset = Counter([node.name])
            else:
                mset = Counter()
                for child in node.children:
                    mset += multisets[child]

            # Store the Counter for children use, but also store sorted tuple for output
            multisets[node] = mset
            sorted_multisets[node] = tuple(sorted(mset.items()))

        return sorted_multisets

    @staticmethod
    def _compute_canonical_forms(root):
        canonical_map = {}

        def recurse(node):
            '''
            Recursively compute a canonical string representation of a subtree rooted at `node`,
            such that isomorphic subtrees (with possibly duplicated leaves) get the same string.
            '''
            if node.is_leaf():
                canon = node.name
            else:
                child_encodings = sorted([recurse(c) for c in node.children])
                canon = '(' + ','.join(child_encodings) + ')'
            canonical_map[node] = canon
            return canon

        recurse(root)
        return canonical_map
    
    @staticmethod
    def _insert_reticulation(G, src, tgt):
        u = f'retic_{src}_{tgt}'
        G.add_node(u, label='Ret.')
        # Keep the length attribute if present
        old_edge = G[src][tgt]
        G.remove_edge(src, tgt)
        if 'length' in old_edge:
            G.add_edge(src, u, length=old_edge['length'])
        else:
            G.add_edge(src, u)
        G.add_edge(u, tgt)
        return u

    @staticmethod
    def _simplify_redundant_reticulations(G):
        to_remove = []
        for node in list(G.nodes):
            preds = list(G.predecessors(node))
            succs = list(G.successors(node))

            # Skip autopolyploidy nodes - they represent WGD events and should be preserved
            if G.nodes[node].get('autopolyploidy', False):
                continue

            # Here we can't use `is_reticulation_node` because it checks in-degree > 1,
            # which is still not the case for these temporary reticulation nodes (to be removed).
            if 'ete' not in G.nodes[node] and len(preds) == 1 and len(succs) == 1:
                parent, child = preds[0], succs[0]
                if not G.has_edge(parent, child):
                    G.add_edge(parent, child)
                to_remove.append(node)

        G.remove_nodes_from(to_remove)

    @staticmethod
    def _get_single_child(G, retnode):
        children = list(G.successors(retnode))
        x = len(children)
        if x != 1:
            raise ValueError(f'Reticulation node {retnode} should have exactly one child, found {x}')
        return children[0]

    @classmethod
    def multree_to_dag(cls, tree, strict=True):
        '''
        An implementation of the HOLM2006 algorithm to convert an ete3 multree to a networkx DAG phylogenetic network.
        This method should NOT modify the original tree.

        Optimizations:
        - Group nodes at each height by canonical form (or multiset) instead of doing O(k^2) pairwise comparisons.
        - Cache sorted multiset representations used for ordering.
        '''
        # Compute heights and multisets (same as before)
        height_map = cls._compute_height_map(tree)

        # Create the initial networkx graph
        G = cls.tree_to_dag(tree)

        # Prepare canonical map or use multisets, depending on strict flag
        if strict:
            # Expect compute_canonical_forms to return a dict: {node: canon_str}
            canonical_map = cls._compute_canonical_forms(tree)
            def key_fn(n):
                return canonical_map[n]
        else:
            # Precompute a stable, hashable representation for each node's multiset
            sorted_multisets = cls._compute_sorted_multisets(tree)
            def key_fn(n):
                return sorted_multisets[n]

        # Process heights from root down to leaves (same order as before)
        for h in sorted(height_map.keys(), reverse=True):
            nodes_at_height = height_map[h]
            # Sort nodes in a stable way using the same key as original code used for ordering.
            # (Original sorted by sorted(multisets[n].items()); for strict True canonical_map ordering is used implicitly.)
            nodes_at_height = sorted(nodes_at_height, key=lambda n: key_fn(n))

            # Group nodes by their key (canonical string or sorted multiset tuple)
            groups = defaultdict(list)
            for n in nodes_at_height:
                groups[key_fn(n)].append(n)

            # For each group with >1 node, pick the first as the retained node and treat the rest as isomorphs
            for group_nodes in groups.values():
                if len(group_nodes) <= 1:
                    continue

                leader = group_nodes[0]  # the node we keep and attach the reticulation above
                # If leader has no parent, we cannot insert a reticulation above it (mirrors previous behavior)
                parent_leader = leader.up
                if not parent_leader:
                    # This mirrors previous code that only inserts reticulation when parent exists.
                    # Skip group if there's no valid parent for leader.
                    continue

                pid = id(parent_leader)
                nid_leader = id(leader)

                # Count distinct parents among alive nodes for event classification
                alive_parents = [pid]  # leader's parent (leader is always alive at this point)

                # Insert reticulation node above leader in G
                retic_node = cls._insert_reticulation(G, pid, nid_leader)

                # For every other isomorphic node in the group, reroute its parent's edge(s) to the retic node
                for iso in group_nodes[1:]:
                    iso_parent = iso.up
                    iso_pid = id(iso_parent) if iso_parent else None

                    # Track parent only if this iso is still alive in the DAG
                    if iso_pid and id(iso) in G:
                        alive_parents.append(iso_pid)

                    # Detect autopolyploidy: identical subtrees under the same parent
                    if iso_pid == pid:
                        G.nodes[retic_node]['autopolyploidy'] = True

                    # If the DAG contains the edge we expect, replace it so the parent's edge now goes to the retic node
                    if iso_pid and G.has_edge(iso_pid, id(iso)):
                        G.remove_edge(iso_pid, id(iso))
                        G.add_edge(iso_pid, retic_node)

                    # Remove the isomorphic subtree nodes from the DAG (exactly as before)
                    to_delete = [id(n) for n in iso.traverse()]
                    G.remove_nodes_from(to_delete)

                # Store event counts: auto = copies under same parent, allo = connections between parents
                k = len(alive_parents)
                d = len(set(alive_parents))
                G.nodes[retic_node]['auto_events'] = k - d
                G.nodes[retic_node]['allo_events'] = d - 1

        # Remove reticulation nodes that have exactly one predecessor and one successor
        cls._simplify_redundant_reticulations(G)
        return G

    @classmethod
    def polyphest_multree_to_dag(cls, tree, threshold=0.2, normalize=True):
        '''
        Polyphest's implementation of the HOLM2006 algorithm:

        Converts an ETE3 MUL-tree to a DAG by merging near-isomorphic subtrees
        based on graph edit distance <= threshold.

        This follows the structure of polyphest's folding.convert_multree_to_network(),
        but uses ETE3 + NetworkX instead of PhyNetPy + Polyphest's helper functions.
        '''
        # Height map {height: [nodes at that height]}
        height_map = cls._compute_height_map(tree)

        # DAG copy of the tree (node IDs are id(ete3_node))
        G = cls.tree_to_dag(tree)

        # Helper: get set of leaf names under node
        def get_leaf_set(n):
            return frozenset(leaf.name for leaf in n.iter_leaves())

        # Near-isomorphic check (polyphest-style)
        def is_near_isomorphic(n1, n2):
            # Quick leaf-set filter
            if get_leaf_set(n1) != get_leaf_set(n2):
                return False

            H1 = cls.tree_to_dag(n1)
            H2 = cls.tree_to_dag(n2)
            dist = nx.graph_edit_distance(
                H1, H2,
                node_match=lambda u, v: u.get('label') == v.get('label')
            )
            if normalize:
                dist /= (H1.number_of_nodes() + H2.number_of_nodes())
            return dist <= threshold

        # Process from top (root) to bottom (leaves)
        for h in sorted(height_map.keys(), reverse=True):
            nodes_to_check = list(height_map[h])
            while nodes_to_check:
                node = nodes_to_check.pop(0)
                isomorphic_trees = [node]
                remaining = []

                for other in nodes_to_check:
                    if is_near_isomorphic(node, other):
                        isomorphic_trees.append(other)
                    else:
                        remaining.append(other)

                nodes_to_check = remaining

                if len(isomorphic_trees) > 1:
                    leader = isomorphic_trees[0]
                    parent_leader = leader.up
                    if not parent_leader:
                        continue  # Can't insert reticulation above root

                    pid = id(parent_leader)
                    nid_leader = id(leader)

                    # Count distinct parents among alive nodes for event classification
                    alive_parents = [pid]  # leader's parent

                    retic_node = cls._insert_reticulation(G, pid, nid_leader)

                    for iso in isomorphic_trees[1:]:
                        iso_parent = iso.up
                        iso_pid = id(iso_parent) if iso_parent else None

                        # Track parent only if this iso is still alive in the DAG
                        if iso_pid and id(iso) in G:
                            alive_parents.append(iso_pid)

                        # Detect autopolyploidy: identical subtrees under the same parent
                        if iso_pid == pid:
                            G.nodes[retic_node]['autopolyploidy'] = True

                        if iso_parent and G.has_edge(iso_pid, id(iso)):
                            G.remove_edge(iso_pid, id(iso))
                            G.add_edge(iso_pid, retic_node)
                        G.remove_nodes_from([id(n) for n in iso.traverse()])

                    # Store event counts
                    k = len(alive_parents)
                    d = len(set(alive_parents))
                    G.nodes[retic_node]['auto_events'] = k - d
                    G.nodes[retic_node]['allo_events'] = d - 1

        cls._simplify_redundant_reticulations(G)
        return G

    # ───── MODULE: From exact network to (MUL-)tree ─────

    @classmethod
    def dag_to_multree(cls, G, match_dag=False):
        '''
        Convert a DAG (from multree or similar) back into an ete3 Tree.
        This handles multi-labeled trees and reticulations.
        '''
        def clone_node(orig_node_data):
            node = TreeNode()
            label = orig_node_data.get('name') or orig_node_data.get('label') or None
            if label:
                node.name = label
            return node

        # Identify reticulation nodes by having multiple parents
        retic_nodes = {n for n in G.nodes if cls.is_reticulation_node(G, n)}

        # Identify root (node with in-degree 0)
        root_id = cls.get_dag_root(G)
        # Root must not be a reticulation node
        if root_id in retic_nodes:
            raise ValueError('Root node cannot be a reticulation node.')

        def build_tree(node_id):
            node_data = G.nodes[node_id]

            # Skip reticulation nodes and attach their single child directly
            if node_id in retic_nodes:
                return build_tree( ReticulateTree._get_single_child(G, node_id) )

            # Normal node
            curr_node = clone_node(node_data)

            # Attach to DAG node if requested
            if match_dag:
                G.nodes[node_id]['ete'] = curr_node
                G.nodes[node_id]['label'] = curr_node.name

            for child_id in G.successors(node_id):
                child_subtree = build_tree(child_id)
                curr_node.add_child(child_subtree)

            return curr_node

        # Post-process: tag reticulation nodes (only if match_dag)
        if match_dag:
            for node in retic_nodes:
                G.nodes[node].pop('ete', None)
                G.nodes[node]['label'] = 'Ret.'

        return build_tree(root_id)

    # ───── MODULE: From H-enriched newick to network ─────

    @staticmethod
    def _find_reticulation_defs(G):
        defs = {}
        for n, data in G.nodes(data=True):
            name = data.get('label', '')
            if name.startswith('#') and G.out_degree(n) > 0:
                hid = name[1:]  # strip the '#' prefix
                if hid in defs:
                    raise ValueError(f"Duplicate definition for reticulation id '{hid}' found.")
                defs[hid] = n
        return defs

    @staticmethod
    def _link_and_track_references(G, defs):
        renamed_targets = set()
        for n, data in list(G.nodes(data=True)):
            name = data.get('label', '')
            if name.startswith('#') and G.out_degree(n) == 0:  # leaf reference
                hid = name[1:]  # strip the '#' prefix
                if hid not in defs:
                    continue

                parent = next(G.predecessors(n), None)
                if parent is None:
                    continue

                edata = dict(G[parent][n])  # preserve edge data (branch length)
                G.remove_node(n)  # remove leaf reference

                target = defs[hid]
                if G.out_degree(target) > 1:  # cladistic case
                    G.add_edge(parent, target, **edata)
                    new_node = hid + '_ref'
                    G.add_node(new_node)
                    for child in list(G.successors(target)):
                        G.add_edge(new_node, child, **G[target][child])
                        G.remove_edge(target, child)
                    G.add_edge(target, new_node)  # structural only
                else:  # direct case
                    G.add_edge(parent, target, **edata)

                renamed_targets.add(hid)

        return renamed_targets

    @classmethod
    def enewick_to_dag(cls, enewick_str):
        '''
        Convert enriched Newick with #H reticulation markers into a DAG.
        - Branch lengths from #H reference leaves are preserved on the edge to the reticulation node.
        - Reticulation -> new internal node edges have no branch length (structural only).
        '''
        # Step 1: parse as-is into DAG
        tree = Tree(enewick_str, format=1) # format for Newick TBD !!!
        G = cls.tree_to_dag(tree)

        # Step 2: find definition nodes (internal nodes ending in #H)
        defs = cls._find_reticulation_defs(G)

        # Step 3: link references and track retnodes
        renamed_targets = cls._link_and_track_references(G, defs)

        # Step 4: rename all reticulation targets
        for hid in renamed_targets:
            target = defs[hid]
            new_tgt_key = 'retic_' + hid
            nx.relabel_nodes(G, {target: new_tgt_key}, copy=False)
            G.nodes[new_tgt_key]['label'] = 'Ret.'  # mark as reticulation

        return G

    def to_enewick(self, filename=None):
        G = self.dag
        retnodes = set(self.retnodes)  # reticulation nodes (no numbering)
        
        # Assign new #H IDs to reticulation nodes
        hid_map = {node: str(i) for i, node in enumerate(sorted(retnodes))}
        
        # Find root
        root = self.get_dag_root(G)

        defined = set()  # keep track of reticulations already defined

        def build(node, parent=None):
            # If reticulation node
            if node in retnodes:
                hid = hid_map[node]
                if node not in defined and (parent is None or parent not in retnodes):
                    # First encounter from main parent → definition
                    defined.add(node)
                    children = list(G.successors(node))
                    if not children:
                        # Reticulation leaf alone
                        return f"#H{hid}"
                    inner = ",".join(build(c, node) for c in children)
                    return f"({inner})#H{hid}"
                else:
                    # Reference form (from other parents)
                    return f"#H{hid}"

            # Non-reticulation node
            children = list(G.successors(node))
            if not children:
                # leaf: output label or node id
                label = G.nodes[node].get("label", str(node))
                return label
            
            label = G.nodes[node].get("label", "")
            inner = ",".join(build(c, node) for c in children)
            return f"({inner}){label}"

        newick_str = build(root) + ";"

        if filename:
            with open(filename, "w") as f:
                f.write(newick_str)
        return newick_str


    # ───── MODULE: Visualization ─────

    def __str__(self):
        return self.tree_str
    
    def visualize(self, filename=None, uid_labeled=False):
        G = self.dag
        depths = self.compute_depths(G)

        # Annotate graph with depth-based 'layer' for multipartite layout
        for node, depth in depths.items():
            G.nodes[node]['layer'] = depth

        pos = nx.multipartite_layout(G, subset_key='layer')

        # Color and shape encoding
        node_colors = []
        node_shapes = {}
        for node in G.nodes():
            if self.is_reticulation_node(G, node):
                node_shapes[node] = 'D'
                node_colors.append('orange')
            elif G.out_degree(node) == 0:
                node_shapes[node] = 'o'
                node_colors.append('lightgreen')
            else:
                node_shapes[node] = 'o'
                node_colors.append('lightblue')

        for shape in set(node_shapes.values()):
            shaped_nodes = [n for n in G.nodes if node_shapes[n] == shape]
            nx.draw_networkx_nodes(
                G, pos, nodelist=shaped_nodes,
                node_shape=shape,
                node_color=[node_colors[list(G.nodes).index(n)] for n in shaped_nodes],
                node_size=1000
            )

        nx.draw_networkx_edges(G, pos, arrows=True)

        # Labels (with uid fallback)
        labels = {}
        for n, d in G.nodes(data=True):
            label = d.get('label', None)
            labels[n] = str(n) if uid_labeled and label == '' else label
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=10)

        plt.axis('off')
        plt.tight_layout()
        if filename:
            plt.savefig(filename, dpi=690)
            print(f'Graph saved to {filename}')
        plt.show()

    def interact(self, filename, uid_labeled=False, launch=True):
        '''
        Interactive DAG visualization using pyvis.
        Emulates the structure and appearance of the matplotlib-based layout.
        '''
        filename = Path(filename)
        net = Network(height='800px', width='100%', directed=True, layout=True)
        net.toggle_physics(False)  # Static layout for hierarchy-like appearance

        G = self.dag
        depths = self.compute_depths(G)

        # Add nodes
        for node, data in G.nodes(data=True):
            label = data.get('label', '')
            if not uid_labeled and label == '':
                label = ' '

            node_type = 'internal'
            shape = 'ellipse'
            color = 'lightblue'

            if self.is_reticulation_node(G, node):
                node_type = 'reticulation'
                shape = 'diamond'
                color = 'orange'
            elif G.out_degree(node) == 0:
                node_type = 'leaf'
                color = 'lightgreen'

            level = depths.get(node, 0)

            net.add_node(
                str(node),
                label=label,
                title=f'{node_type} node: {node}',
                level=level,
                shape=shape,
                color=color,
                size=20
            )

        # Add edges
        for u, v in G.edges():
            net.add_edge(str(u), str(v), width=1)

        # Save and optionally open
        net.save_graph(str(filename))
        print(f'Interactive graph saved to {filename}.')
        if launch:
            try:
                net.show(str(filename))
            except Exception as e:
                print(f'Failed to display HTML. Reason: {e}\nTrying to open in browser...')
                try:
                    webbrowser.open(str(filename))
                    # Wait for the user to come back
                    input('Press Enter to continue...')
                except:
                    print('Failed to open in browser. Please open the file manually.')
                    pass

    # ───── MODULE: Comparative features ─────

    def get_reticulation_nodes(self):
        '''
        Get all reticulation nodes in the DAG.
        '''
        return [n for n in self.dag.nodes() if self.is_reticulation_node(self.dag, n)]
    
    def get_reticulation_count(self):
        '''
        Count total reticulation events (auto + allo).
        Uses event counts stored during folding when available,
        falls back to in_degree - 1 for extended Newick input.
        '''
        count = 0
        for n in self.retnodes:
            if 'auto_events' in self.dag.nodes[n]:
                count += self.dag.nodes[n]['auto_events'] + self.dag.nodes[n]['allo_events']
            else:
                # Extended Newick input: no folding, use in-degree
                count += max(self.dag.in_degree(n) - 1, 1)
        return count

    def get_auto_event_count(self):
        '''Count autopolyploidization events (copies under same parent).'''
        return sum(self.dag.nodes[n].get('auto_events', 0) for n in self.retnodes)

    def get_allo_event_count(self):
        '''Count allopolyploidization events (copies under different parents).'''
        count = 0
        for n in self.retnodes:
            if 'allo_events' in self.dag.nodes[n]:
                count += self.dag.nodes[n]['allo_events']
            else:
                # Extended Newick input: all events assumed allo
                count += max(self.dag.in_degree(n) - 1, 1)
        return count
    
    def get_reticulation_leaves(self):
        '''
        Get all leaves that are descendants of reticulation nodes.
        '''
        leaves = {}
        for ret in self.retnodes:
            desc = nx.descendants(self.dag, ret)
            leaf_labels = [
                self.dag.nodes[n].get('label', str(n))
                for n in desc
                if self.dag.out_degree(n) == 0
            ]
            leaves[ret] = leaf_labels
        return leaves
    
    def get_reticulation_sisters(self):
        '''
        For each reticulation node, return its two sister clades (one from each parent).
        Each clade is a list of ete3 TreeNodes (from the 'ete' attribute in the DAG).
        
        Returns:
            dict: {ret_node_id: [clade1_nodes, clade2_nodes]}
        '''
        def collect_subtree_nodes(node):
            '''Return node as ete3 TreeNodes under the given DAG node.'''
            data = self.dag.nodes[node]
            # Check if n is a ret node, if so, take its only child

            if self.is_reticulation_node(self.dag, node):
                data = self.dag.nodes[ self._get_single_child(self.dag, node) ]
            try:
                return data['ete']
            except KeyError:
                raise ValueError(f'No ete3 node found for DAG node {node}')
            
        sisters = {}
        for ret in self.retnodes:
            parents = list(self.dag.predecessors(ret))
            if len(parents) != 2:
                if self.dag.nodes[ret].get('autopolyploidy', False):
                    # Autopolyploidy node has 1 parent - this is expected
                    pass
                else:
                    if _WARN_MALFORMED_RETIC:
                        print(f'Warning: Reticulation node {ret} does not have exactly two parents: {parents}')
                        print('This indicates that a Strict MUL-to-Net mode was used on a malformed inference.')

            clades = []
            for parent in parents:
                # Identify the specific child of the parent that is the reticulation node
                # We want the 'sibling' subtrees, i.e., other children of the parent
                sibling_children = [child for child in self.dag.successors(parent) if child != ret]
                # Collect all ete3 nodes in the sibling clades
                clade_nodes = []
                for child in sibling_children:
                    clade_nodes.extend(collect_subtree_nodes(child))
                clades.append(clade_nodes)

            sisters[ret] = clades

        # Convert to a more user-friendly format: list of sets of names
        sisters = { ret:
            [set([n.name for n in side]) for side in sister_pair]
            for ret, sister_pair in sisters.items()
        }

        return sisters
    
    def get_leaf_counts(self):
        return dict(Counter(leaf.name for leaf in self.tree.iter_leaves()))

    def measure(self, printout=False):
        measurements = {
            'reticulation_count': self.get_reticulation_count(),
            'reticulation_leaves': self.get_reticulation_leaves(),
            'reticulation_sisters': self.get_reticulation_sisters(),
            'leaf_counts': self.get_leaf_counts(),
        }
        if printout:
            for key, value in measurements.items():
                k_str = key.replace('_', ' ').title()
                print(f'{k_str}: {value}')
        return measurements

    def get_edit_distance(self, other: 'ReticulateTree', normalize=True) -> float:
        '''
        Compute the graph edit distance between two ReticulateTree instances.
        This is a placeholder for future implementation.
        '''
        distance = next(nx.optimize_graph_edit_distance(
            self.dag, other.dag,
            node_match=lambda u, v: u.get('label') == v.get('label')
        ))
        if normalize:
            normalization = max(
                len(self.dag.nodes) + len(self.dag.edges),
                len(other.dag.nodes) + len(other.dag.edges)
            )
            distance /= normalization
        return distance
    
    def get_edit_distance_multree(self, other: 'ReticulateTree', normalize=True, timeout=300,
                                  canonical=True) -> float:
        '''
        Compute graph edit distance on MUL-trees (before folding to networks).
        This compares the tree structures directly without network folding.

        Args:
            other: Another ReticulateTree instance
            normalize: If True, normalize by max tree size
            timeout: Max seconds per comparison (default 300s = 5 min).
                     Returns NaN if exceeded. GED is NP-hard so large trees
                     that don't finish in 5 min won't finish in 24 hours either.
            canonical: If True (default), children are visited in canonical order so the
                     result cannot depend on the child ordering of the input Newick.
                     Only set False to reproduce pre-2026-08-05 numbers.

        Returns:
            Edit distance between the two MUL-trees, or NaN if timed out
        '''
        import signal

        def _timeout_handler(signum, frame):
            raise TimeoutError("GED computation timed out")

        # Convert trees to graphs (without folding - just tree structure)
        def tree_to_simple_graph(tree_obj):
            """Convert ete3 tree to NetworkX graph preserving tree structure.

            optimize_graph_edit_distance below returns the FIRST complete edit path its
            greedy search finds, not the minimum, and that search enumerates candidates
            in node INSERTION order. So without canonical ordering the score depends on
            the child order of the input Newick, which is phylogenetically meaningless
            and differs between the programs that wrote each tree: two identical trees
            can score up to ~0.6 instead of 0.

            Same idea as _compute_canonical_forms above and _iso_leaf_map in
            apply_diploidization.py, order both sides canonically so they align.
            Does not modify the tree.
            """
            G = nx.DiGraph()

            if not canonical:
                # Kept verbatim so canonical=False reproduces pre-fix numbers exactly
                # (traverse() defaults to levelorder, not the preorder walk below)
                for node in tree_obj.tree.traverse():
                    node_id = id(node)
                    # Leaf nodes get their species name as label
                    if node.is_leaf():
                        G.add_node(node_id, label=node.name)
                    else:
                        # Internal nodes don't get labels (or get None)
                        G.add_node(node_id, label=None)

                    # Add edge from parent to this node
                    if not node.is_root():
                        parent_id = id(node.up)
                        G.add_edge(parent_id, node_id)

                return G

            canon = {}

            def canonical_form(node):
                if node.is_leaf():
                    form = node.name
                else:
                    form = '(' + ','.join(sorted(canonical_form(c) for c in node.children)) + ')'
                canon[node] = form
                return form

            canonical_form(tree_obj.tree)

            def add_subtree(node):
                node_id = id(node)
                G.add_node(node_id, label=node.name if node.is_leaf() else None)

                if not node.is_root():
                    G.add_edge(id(node.up), node_id)

                # Ties (identical sibling subtrees, i.e. autopolyploidy) keep their input
                # order, which is harmless because the subtrees are identical
                for child in sorted(node.children, key=lambda c: canon[c]):
                    add_subtree(child)

            add_subtree(tree_obj.tree)
            return G

        # Get graphs for both trees
        graph1 = tree_to_simple_graph(self)
        graph2 = tree_to_simple_graph(other)

        # GED is NP-hard with exponential memory — skip for very large trees
        MAX_NODES_FOR_GED = 500
        n1, n2 = len(graph1.nodes), len(graph2.nodes)
        if max(n1, n2) > MAX_NODES_FOR_GED:
            print(f"  Skipping MUL-tree GED: trees too large ({n1}, {n2} nodes > {MAX_NODES_FOR_GED} limit)")
            return float('nan')

        # Compute edit distance with timeout. SIGALRM is POSIX-only, so on Windows
        # we run without one (fine for the small trees used off-cluster).
        use_alarm = bool(timeout) and hasattr(signal, 'SIGALRM')
        if use_alarm:
            old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
            signal.alarm(timeout)
        try:
            distance = next(nx.optimize_graph_edit_distance(
                graph1, graph2,
                node_match=lambda u, v: u.get('label') == v.get('label')
            ))
        except TimeoutError:
            n1, n2 = len(graph1.nodes), len(graph2.nodes)
            print(f"  WARNING: MUL-tree GED timed out after {timeout}s "
                  f"(tree sizes: {n1}, {n2} nodes). Returning NaN.")
            return float('nan')
        finally:
            if use_alarm:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)

        if normalize:
            normalization = max(
                len(graph1.nodes) + len(graph1.edges),
                len(graph2.nodes) + len(graph2.edges)
            )
            if normalization > 0:
                distance /= normalization

        return distance

    def _mu_taxa(self, G=None):
        '''
        Sorted taxon labels of the folded network's leaves.

        The mu-representation indexes its vectors by taxon, so the network must
        carry exactly one leaf per taxon. Folding collapses the MUL-tree copies,
        so this holds for a correctly folded network and a violation means the
        network is not a valid comparison target.
        '''
        G = self.dag if G is None else G
        labels = [G.nodes[n].get('label') for n, d in G.out_degree() if d == 0]
        if any(label is None for label in labels):
            raise ValueError('Network has an unlabeled leaf; cannot index mu-vectors by taxon.')
        if len(labels) != len(set(labels)):
            duplicated = sorted({x for x in labels if labels.count(x) > 1})
            raise ValueError(
                f'Network has repeated leaf labels {duplicated}; the mu-representation '
                'requires one leaf per taxon. Was this folded?'
            )
        return sorted(labels)

    @staticmethod
    def _suppress_elision_nodes(G):
        '''
        Drop nodes with in-degree 1 and out-degree 1.

        Folding can leave these behind. They carry no phylogenetic information,
        and semi-binary requires every internal tree node to have out-degree 2,
        so they must go before the class preconditions or the representation
        are evaluated.
        '''
        H = G.copy()
        changed = True
        while changed:
            changed = False
            for node in list(H.nodes()):
                if H.in_degree(node) == 1 and H.out_degree(node) == 1:
                    parent = next(H.predecessors(node))
                    child = next(H.successors(node))
                    H.remove_node(node)
                    if parent != child:
                        H.add_edge(parent, child)
                    changed = True
                    break
        return H

    def _mu_network(self):
        '''The folded network in the form the mu-representation is defined on.'''
        return self._suppress_elision_nodes(self.dag)

    def _check_mu_preconditions(self, G, label='network'):
        '''
        Enforce the hypotheses of Theorem 1 of arXiv:2412.05107.

        Raises rather than returning a number the theorem does not licence. A
        network outside the class can share its representation with a
        non-isomorphic one, so a distance of 0 would not mean what it appears
        to mean.
        '''
        reticulations = {u for u, d in G.in_degree() if d >= 2}

        stacked = [(p, c) for p, c in G.edges()
                   if p in reticulations and c in reticulations]
        if stacked:
            raise ValueError(
                f'The {label} is not stack-free: {len(stacked)} reticulation(s) are '
                'children of reticulations. Theorem 1 of the modified '
                'mu-representation does not apply, so no distance is reported.'
            )

        offending = []
        for node in G.nodes():
            in_degree, out_degree = G.in_degree(node), G.out_degree(node)
            if out_degree == 0:
                if in_degree != 1:
                    offending.append((in_degree, out_degree))
            elif node in reticulations:
                if out_degree != 1:
                    offending.append((in_degree, out_degree))
            elif out_degree != 2:
                offending.append((in_degree, out_degree))
        if offending:
            raise ValueError(
                f'The {label} is not semi-binary: {len(offending)} node(s) violate the '
                f'degree constraints, e.g. (in,out) = {offending[:5]}. Tree nodes must '
                'have out-degree 2 and reticulations out-degree 1.'
            )

    @staticmethod
    def _suppress_node(G, node):
        '''Bypass a node left with one parent and one child, or a stalled root.'''
        if node not in G:
            return
        if G.in_degree(node) == 1 and G.out_degree(node) == 1:
            parent = next(G.predecessors(node))
            child = next(G.successors(node))
            G.remove_node(node)
            if parent != child:
                G.add_edge(parent, child)
        elif G.in_degree(node) == 0 and G.out_degree(node) == 1:
            G.remove_node(node)

    @classmethod
    def _pick_reducible_pair(cls, G):
        '''
        Pick one cherry or reticulated cherry from G, reducing it in place.

        A cherry is two leaves sharing a parent. A reticulated cherry is a leaf
        whose parent is a reticulation that is itself a child of the parent of
        another leaf. Returns the pair of leaf labels, or None if neither
        exists, in which case the reduction is stuck.
        '''
        leaves = [n for n, d in G.out_degree() if d == 0]

        by_parent = defaultdict(list)
        for leaf in leaves:
            if G.in_degree(leaf) == 1:
                by_parent[next(G.predecessors(leaf))].append(leaf)
        for parent, siblings in by_parent.items():
            if len(siblings) >= 2:
                removed, kept = siblings[0], siblings[1]
                pair = (G.nodes[removed].get('label'), G.nodes[kept].get('label'))
                G.remove_node(removed)
                cls._suppress_node(G, parent)
                return pair

        reticulations = {n for n, d in G.in_degree() if d >= 2}
        for leaf in leaves:
            if G.in_degree(leaf) != 1:
                continue
            reticulation = next(G.predecessors(leaf))
            if reticulation not in reticulations:
                continue
            for parent in list(G.predecessors(reticulation)):
                partners = [c for c in G.successors(parent)
                            if c != reticulation and G.out_degree(c) == 0]
                if not partners:
                    continue
                pair = (G.nodes[leaf].get('label'), G.nodes[partners[0]].get('label'))
                G.remove_edge(parent, reticulation)
                cls._suppress_node(G, reticulation)
                cls._suppress_node(G, parent)
                return pair

        return None

    def orchard_certificate(self):
        '''
        Try to reduce the network to a single leaf by cherry picking.

        Returns (True, sequence) when it succeeds, where sequence is the list of
        picked pairs and therefore a certificate that the network is orchard.
        Returns (False, None) when the reduction gets stuck.

        The search is greedy. A success is a proof because the sequence is
        exhibited. A failure is only a proof if cherry picking is confluent for
        this class, so treat (False, None) as "no certificate found" rather than
        as established non-orchard.
        '''
        G = self._mu_network()
        sequence = []
        while True:
            leaves = [n for n, d in G.out_degree() if d == 0]
            if len(leaves) <= 1:
                return True, sequence
            pair = self.__class__._pick_reducible_pair(G)
            if pair is None:
                return False, None
            sequence.append(pair)

    def is_orchard(self):
        '''True when a cherry picking sequence reducing the network was found.'''
        return self.orchard_certificate()[0]

    def get_modified_mu_representation(self):
        '''
        Modified mu-representation, Definition 1 of arXiv:2412.05107.

            mu-bar(v) = ( in_degree(v), mu_1(v), ..., mu_n(v) )

        where mu_i(v) counts directed paths from v to leaf i. The multiset runs
        over ALL nodes, reticulations included.

        This is not the extended mu-representation of Cardona et al. 2024. That
        one prepends the number of paths to a reticulation instead of the
        in-degree, excludes reticulation nodes, and encodes only BINARY orchard
        networks. Figure 1 of arXiv:2412.05107 exhibits two non-isomorphic
        semi-binary stack-free networks sharing an extended mu-representation,
        and MUL-tree folding produces semi-binary networks, so the extended
        form is not sound here.

        The mu-vectors themselves come from the authors' `phylonetwork`
        package, so the path counts are the published implementation. Making
        them semi-binary ready needs only the in-degree prepended, which the
        package already exposes. Ours is the suppression of degree-2 nodes that
        semi-binary requires, the class preconditions, and the distance.
        '''
        G = self._mu_network()
        self._check_mu_preconditions(G)

        try:
            import phylonetwork
        except ImportError as exc:
            raise ImportError(
                'The mu-vectors are computed by the phylonetwork reference '
                'implementation (Cardona et al.), which is not installed. '
                'Install it with: pip install phylonetwork==2.2.2'
            ) from exc

        network = phylonetwork.PhylogeneticNetwork(
            eNewick=self.__class__(G).to_enewick()
        )
        return [tuple([network.in_degree(u)] + [int(x) for x in network.mu_dict[u]])
                for u in network.nodes()]

    def _mu_path_counts(self, G):
        '''
        Standard mu-vector of every node: mu_i(v) counts paths from v to leaf i.

        Nodes are visited in reverse topological order so every child is
        complete before its parent.
        '''
        taxa = self._mu_taxa(G)
        position = {taxon: i for i, taxon in enumerate(taxa)}

        vectors = {}
        for node in reversed(list(nx.topological_sort(G))):
            vector = [0] * len(taxa)
            for child in G.successors(node):
                child_vector = vectors[child]
                for i in range(len(taxa)):
                    vector[i] += child_vector[i]
            if G.out_degree(node) == 0:
                vector[position[G.nodes[node]['label']]] += 1
            vectors[node] = vector

        return vectors, taxa

    def get_mu_root_vector(self):
        '''
        Taxa and the root's per-taxon path counts.

        mu_i(root) counts root-to-leaf-i paths, which equals the number of
        copies of taxon i in the MUL-tree, so this is the copy-number (ploidy)
        vector of the network.
        '''
        G = self._mu_network()
        vectors, taxa = self._mu_path_counts(G)
        root = self.__class__.get_dag_root(G)
        return taxa, tuple(vectors[root])

    def get_mu_distance(self, other: 'ReticulateTree', normalize=True) -> float:
        '''
        Extended mu-distance between two folded networks.

        d(N1, N2) = |mu-bar(N1) symmetric-difference mu-bar(N2)|, Theorem 1 of
        arXiv:2412.05107: for N1 semi-binary stack-free orchard and N2
        semi-binary stack-free, this is 0 exactly when N1 and N2 are isomorphic.

        Uses the MODIFIED mu-representation, which prepends each node's
        in-degree to its standard mu-vector and ranges over all nodes. The
        extended mu-representation of Cardona et al. 2024, which the
        phylonetwork package implements, encodes only BINARY orchard networks
        and is not sound on the semi-binary networks MUL-tree folding produces.
        Our standard mu-vectors are cross-checked against that package.

        Networks outside the semi-binary stack-free class, and pairs on
        different taxon sets, raise rather than returning a number the theorem
        does not licence.

        With normalize=True the count is divided by the total number of
        mu-vectors, giving a value in [0, 1] that can be averaged across
        networks of different sizes.
        '''
        network_self = self._mu_network()
        network_other = other._mu_network()
        self._check_mu_preconditions(network_self, 'first network')
        other._check_mu_preconditions(network_other, 'second network')

        taxa_self = self._mu_taxa(network_self)
        taxa_other = other._mu_taxa(network_other)
        if taxa_self != taxa_other:
            only_self = sorted(set(taxa_self) - set(taxa_other))
            only_other = sorted(set(taxa_other) - set(taxa_self))
            raise ValueError(
                'Cannot compute mu-distance between networks on different taxa. '
                f'Only in first: {only_self}. Only in second: {only_other}.'
            )

        mu_self = self.get_modified_mu_representation()
        mu_other = other.get_modified_mu_representation()
        set_self, set_other = set(mu_self), set(mu_other)
        distance = float(sum(1 for vector in mu_self + mu_other
                             if vector not in set_self or vector not in set_other))

        if normalize:
            normalization = len(mu_self) + len(mu_other)
            if normalization == 0:
                raise ValueError('Both networks have empty mu-representations.')
            distance /= normalization

        return distance

    def __sub__(self, other: 'ReticulateTree') -> float:
        '''
        Overriding the minus operator to compute the edit distance between two ReticulateTree instances.
        '''
        if not isinstance(other, ReticulateTree):
            raise TypeError('Subtraction is only supported between ReticulateTree instances.')
        return self.get_edit_distance(other)

