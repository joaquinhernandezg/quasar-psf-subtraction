import numpy as np
import glob
import argparse
from astropy.table import Table
from collections import defaultdict
from quasar_psf_subtraction.utils.hsc_footprints import filter_catalog_by_footprint


def main():
    parser = argparse.ArgumentParser(description="A script in my Python module.")
    parser.add_argument('input_catalog', type=str, help="Path to the configuration file")
    parser.add_argument('hsc_release', type=str, help="Release of the HSC data")
    parser.add_argument('ra_column', type=str, help="Ra column name")
    parser.add_argument('dec_column', type=str, help="Dec column name")
    # optional catalog filename
    parser.add_argument('--output_catalog', type=str, help="Output catalog filename")
    parser.add_argument('--catalog_format', type=str, help="Catalog format")
    parser.add_argument('--radec_units', type=str, default="deg", help="Units of the RA and Dec columns")

    args = parser.parse_args()


    catalog = args.input_catalog
    ra_column = args.ra_column
    dec_column = args.dec_column
    release = args.hsc_release
    output_catalog = args.output_catalog

    if args.catalog_format is not None:
        catalog_format = args.catalog_format
        table = Table.read(catalog, format=catalog_format)
    else:
        table = Table.read(catalog)
    print("Initial number of objects: ", len(table))

    from os.path import dirname, abspath

    module_path = dirname(dirname(dirname(abspath(__file__))))

    if release in ["pdr2_wide", "pdr2_dud", "pdr3_wide", "pdr3_dud"]:
        tract_dir = f"{module_path}/data/hsc_tracts/{release}"
    else:
        raise ValueError("Release not supported")

    table = filter_catalog_by_footprint(table, ra_column, dec_column, tract_dir, radec_units=args.radec_units)

    if output_catalog is None:
        extension = catalog.split(".")[-1]
        output_catalog = catalog.replace(f".{extension}", f"_in_{release}_footprint.{extension}")
    
    if catalog_format is not None:
        table.write(output_catalog, format=catalog_format, overwrite=True)
    else:
        table.write(output_catalog, overwrite=True)
    
    print("Final number of objects: ", len(table))
