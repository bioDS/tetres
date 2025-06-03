"""
Utility for handling IO related things.
Currently only splitting and writing log files according to a clustering.
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from numpy._typing import NDArray

logger = logging.getLogger(__name__)


def write_clustered_logfiles(df: pd.DataFrame, clustering: NDArray[np.int_], output_dir: Path,
                             name: str, k: int, step: int = 1000, _overwrite: bool = False):
    """
    Splits a DataFrame based on clustering labels and writes each subset to a tsv .log file.

    :param step: For fixing the sample column in the log files, increment size between samples
    :param df: DataFrame to split
    :param clustering: List/array/Series of cluster labels (same length as df)
    :param output_dir: Path to output directory
    :param name: Name for output files
    :param k: Number of clusters
    :param _overwrite: If True, existing files will be overwritten
    """
    logger.warning("Writing clustered logfiles...")
    clustering_series = pd.Series(clustering, index=df.index)

    for cluster_id, group in df.groupby(clustering_series):
        output_log_file = output_dir / f"split_k-{k}-c-{cluster_id}-{name}.log"
        if output_log_file.exists() and not _overwrite:
            logger.warning(f"Skipping existing file: {output_log_file}")
            continue

        group = group.copy()
        group["Sample"] = range(0, step * len(group), step)

        group.to_csv(output_log_file, index=False, sep="\t")
        logger.warning(f"Wrote: {output_log_file}")
