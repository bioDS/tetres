"""
Utility functions for reading and writing phylogenetic tree data
WIP: This should probably contain other functions such as write_nexus, etc.
"""
import logging
from pathlib import Path

from tetres.trees.time_trees import TimeTreeSet

logger = logging.getLogger(__name__)


def write_clustered_treesets(clustered_trees: dict[int, list],
                             map_data, output_dir: Path, name: str, k: int,
                             _overwrite: bool = False):
    logger.warning(f"Splitting tree set {name} into {k} subsets...")
    for k_cluster in range(k):
        cur_treeset = TimeTreeSet()
        cur_treeset.map = map_data
        cur_treeset.trees = clustered_trees[k_cluster]

        output_tree_file = output_dir / f"split_k-{k}-c-{k_cluster}-{name}.trees"
        if output_tree_file.exists():
            if _overwrite:
                logger.warning(f"Overwriting existing subset tree file: {output_tree_file}")
                output_tree_file.unlink()
            else:
                logger.warning(f"Skipping existing file: {output_tree_file} "
                               f"(pass _overwrite=True to overwrite)")
                continue

        cur_treeset.write_nexus(file_name=output_tree_file)
