"""
Utility functions for reading and writing phylogenetic tree data
WIP: This should probably end up in the trees module and contain other functions such as
write_nexus, etc.
"""
import logging
from pathlib import Path

from tetres.trees.time_trees import TimeTreeSet

logger = logging.getLogger(__name__)


def write_clustered_treesets(clustered_trees: dict[int, list],
                             map_data,
                             working_dir: Path,
                             name: str,
                             k: int,
                             _overwrite: bool = False):
    logger.warning(f"Splitting tree set {name} into {k} subsets...")
    for k_cluster in range(k):
        cur_treeset = TimeTreeSet()
        cur_treeset.map = map_data
        cur_treeset.trees = clustered_trees[k_cluster]

        output_path = working_dir / f"trees_k-{k}-c-{k_cluster}-{name}.trees"
        if output_path.exists():
            if _overwrite:
                logger.warning(f"Overwriting existing subset tree file: {output_path}")
                output_path.unlink()
            else:
                logger.warning(f"Skipping existing file: {output_path} "
                               f"(pass _overwrite=True to overwrite)")
                continue

        cur_treeset.write_nexus(file_name=output_path)
