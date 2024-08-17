#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-08-13
# @Filename: create_confSummaryS.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import logging
import multiprocessing
import os
import pathlib
from functools import partial

from typing import TYPE_CHECKING

import numpy
import polars
from astropy.table import Table

from sdsstools._vendor.yanny import write_ndarray_to_yanny, yanny
from sdsstools.logger import get_logger
from target_selection.skies import is_valid_sky


if TYPE_CHECKING:
    pass


ROOT_DIR = pathlib.Path(__file__).parents[1]

# List of directories to skip. The path is relative to the root directory.
# Add paths here once all the files in that directory have been processed to
# speed up the process.
SKIP_DIRECTORIES: list[str] = [
    "apo/summary_files/000XXX",
    "apo/summary_files/001XXX",
    "apo/summary_files/002XXX",
    "apo/summary_files/003XXX",
    "apo/summary_files/004XXX",
    "apo/summary_files/005XXX",
    "apo/summary_files/006XXX",
    "apo/summary_files/007XXX",
    "apo/summary_files/008XXX",
    "apo/summary_files/009XXX",
    "apo/summary_files/010XXX",
    "apo/summary_files/011XXX",
    "apo/summary_files/012XXX",
    "apo/summary_files/013XXX",
    "apo/summary_files/014XXX",
    "lco/summary_files/10000XXX"
    "lco/summary_files/10001XXX"
    "lco/summary_files/10002XXX"
    "lco/summary_files/10003XXX"
    "lco/summary_files/10004XXX"
    "lco/summary_files/10005XXX"
    "lco/summary_files/10006XXX"
    "lco/summary_files/10007XXX"
    "lco/summary_files/10008XXX",
]

# List of catalogues against which to validate the sky coordinates.
CATALOGUES: list[str] = ["gaia_dr3_source", "twomass_psc", "tycho2", "twomass_xsc"]

# Name of the fake carton to assign to valid skies.
CARTON: str = "ops_sky_apogee_best_unassigned"

# Connection URI.
DATABASE_URI: str = "postgresql://sdss_user@operations.sdss.org/sdss5db"


log = get_logger("sdsscore.create_confsummarys", use_rich_handler=True)


def read_confSummary(path: str | pathlib.Path) -> tuple[dict, polars.DataFrame]:
    """Reads a configuration summary file and returns the header and data frame."""

    if not os.path.exists(str(path)):
        raise FileNotFoundError(f"{path} does not exist.")

    y = yanny(str(path))
    header = dict(y)

    fibermap = header.pop("FIBERMAP")

    df = polars.DataFrame(fibermap)
    df = df.with_columns(polars.selectors.binary().cast(polars.String()))

    for key, value in header.items():
        try:
            header[key] = int(value)
        except ValueError:
            try:
                header[key] = float(value)
            except ValueError:
                pass

    return header, df.sort(["positionerId", "fiberType"])


def get_confSummaryS_file(file_: pathlib.Path):
    """Creates the ``confSummaryS`` path for a ``confSummary`` file."""

    root, rest = file_.name.split("-")

    return file_.with_name(f"{root}S-{rest}")


def get_new_files():
    """Returns a list of ``confSummary(F)`` files that need processing."""

    # TODO: there may be a faster way to do this.

    new_files: list[pathlib.Path] = []

    # List XXX directories.
    xxx_dirs = ROOT_DIR.glob("**/*XXX")

    # Create the list of directories to skip.
    skip_dirs = [ROOT_DIR / path for path in SKIP_DIRECTORIES]

    for xxx_dir in xxx_dirs:
        if xxx_dir in skip_dirs:
            continue

        files = list(xxx_dir.glob("**/confSummary*.par"))

        for file_ in files:
            root, _ = file_.name.split("-")

            # Skip confSummaryS files
            if root.endswith("S"):
                continue

            # Check if the associated confSummary(F)S file exists
            fileS = get_confSummaryS_file(file_)
            if fileS not in files:
                new_files.append(file_)

    return new_files


def process_file(
    path: pathlib.Path | str,
    database_uri: str = DATABASE_URI,
    overwrite: bool = False,
    create_parquet: bool = True,
):
    """Processes a ``confSummary`` to create the corresponding ``confSummaryS`` file."""

    path = pathlib.Path(path)
    pathS = get_confSummaryS_file(path)

    if pathS.exists() and not overwrite:
        raise ValueError(f"{pathS!s} already exists. Use overwrite=True.")

    pathS.unlink(missing_ok=True)

    log.debug(f"Processing file {path.name} ...")

    header, df = read_confSummary(path)

    # Select only APOGEE fibres that are not on target.
    unassigned = df.filter(
        polars.col.assigned == 0,
        polars.col.on_target == 0,
        polars.col.fiberType == "APOGEE",
        polars.col.fiberId > 0,
        polars.col.ra >= 0,
        polars.col.dec > -999,
    )
    coords = unassigned.select(["ra", "dec"]).to_numpy()

    # Get the skies mask.
    if coords.size > 0:
        try:
            mask = is_valid_sky(
                coords,
                database_uri,
                catalogues=CATALOGUES,
                epoch=header["epoch"],
            )
        except Exception as ee:
            log.error(f"Error validating skies in {path.name}: {ee}")
            log.warning(f"A copy of {path.name} will be used as confSummaryS.")
            mask = numpy.zeros(coords.shape[0], dtype=bool)
    else:
        log.warning(f"No unassigned APOGEE fibers in {path.name}.")
        mask = numpy.array([])

    # Limit unassigned to the columns we care about and replace the program
    # and category columns for valid skies.
    unassigned = unassigned.select(
        [
            "positionerId",
            "fiberType",
            "firstcarton",
            "program",
            "category",
        ]
    ).with_columns(
        firstcarton=polars.Series(numpy.where(mask, CARTON, unassigned["firstcarton"])),
        program=polars.Series(numpy.where(mask, "ops_sky", unassigned["program"])),
        category=polars.Series(numpy.where(mask, "sky_apogee", unassigned["category"])),
    )

    # Join back to the full fibermap df.
    df = df.join(unassigned, on=["positionerId", "fiberType"], how="left")

    df = df.with_columns(
        firstcarton=polars.coalesce(["firstcarton_right", "firstcarton"]),
        program=polars.coalesce(["program_right", "program"]),
        category=polars.coalesce(["category_right", "category"]),
    ).drop(["firstcarton_right", "program_right", "category_right"])

    # Write the new confSummaryS file.
    fibermap = Table.from_pandas(df.to_pandas())
    fibermap["mag"] = fibermap["mag"].astype(numpy.dtype(("f4", 5)))  # type: ignore

    write_ndarray_to_yanny(
        str(pathS),
        [fibermap.filled()],
        structnames=["FIBERMAP"],
        hdr=header,
        enums={"fiberType": ("FIBERTYPE", ("BOSS", "APOGEE", "METROLOGY", "NONE"))},
    )

    if create_parquet:
        parquet_file = pathS.with_suffix(".parquet")

        df_full = df.with_columns(**{kk: polars.lit(vv) for kk, vv in header.items()})
        df_full.write_parquet(parquet_file)


def _process_file_wrapper(path: pathlib.Path, **kwargs):
    """Helper function to wrap `.process_file` in try-except."""

    try:
        return process_file(path, **kwargs)
    except Exception as ee:
        log.error(f"Error processing {path}: {ee}")


def process_summaries(
    database_uri: str = DATABASE_URI,
    n_cpus: int = 4,
    verbose: bool = False,
    create_parquet: bool = True,
):
    """Processes all the ``confSummary`` files in ``sdsscore``.

    Produces a list of ``confSummary`` and ``confSummaryF`` files without
    a corresponding ``confSummaryS`` file. For each of these files runs
    `.process_file` to generate the ``confSummaryS``. The process is parallelised
    using ``n_cpus`` CPUs.

    Parameters
    ----------
    database_uri
        The database URI to use for the sky queries.
    n_cpus
        The number of CPUs to use for parallel processing.
    verbose
        If `True`, will print progress information.
    create_parquet
        If `True`, will create a parquet file for each ``confSummaryS`` file.

    """

    log.set_level(logging.DEBUG if verbose else logging.INFO)
    log.sh.setLevel(logging.DEBUG if verbose else logging.INFO)

    log.info("Retrieving list of new confSummary files to process.")
    new_files = sorted(get_new_files())

    n_new = len(new_files)

    log.info(f"Found {n_new} new confSummary files to process.")

    if verbose:
        if n_new > 1000:
            log.debug("Too many new files to list.")
        else:
            log.debug("New files:")
            for file_ in new_files:
                log.debug(str(file_.name))

    log.info("Processing new files ...")

    partial_func = partial(
        _process_file_wrapper,
        database_uri=database_uri,
        create_parquet=create_parquet,
    )

    with multiprocessing.Pool(n_cpus) as pool:
        pool.map(partial_func, new_files)


if __name__ == "__main__":
    process_summaries(verbose=True)
