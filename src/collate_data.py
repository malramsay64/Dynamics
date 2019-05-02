#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Collate a collection of gsd files into a single file."""

import logging
from pathlib import Path
from typing import Tuple

import click
import pandas

logger = logging.getLogger(__name__)
logging.basicConfig(level="WARNING")


@click.command()
@click.argument("output", type=click.Path(file_okay=True, dir_okay=False))
@click.argument(
    "infiles", nargs=-1, type=click.Path(exists=True, file_okay=True, dir_okay=False)
)
def collate_files(output: Path, infiles: Tuple[Path, ...]) -> None:
    with pandas.HDFStore(output, "w") as dst:
        for file in infiles:
            with pandas.HDFStore(file) as src:
                for key in ["dynamics", "molecular_relaxations"]:
                    try:
                        dst.append(key, src.get(key))
                    except KeyError:
                        logger.warning("%s doesn't contain key %s", file, key)


if __name__ == "__main__":
    collate_files()
