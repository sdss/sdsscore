#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2025-08-13
# @Filename: __main__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

from typing import Annotated

import typer

from sdsscore.confsummary import process_summaries


app = typer.Typer()

confsummary_app = typer.Typer()
app.add_typer(
    confsummary_app,
    name="confsummary",
    help="Handle configuration summary files",
)


@confsummary_app.command()
def update(
    quiet: Annotated[
        bool,
        typer.Option("--quiet", "-q", help="Suppress output messages"),
    ] = False,
):
    """Generate confSummaryS files and Parquet versions of other summary files.

    This command should only be run at Utah.

    """

    process_summaries(verbose=not quiet)
