import numpy as np
import glob
import argparse
from astropy.table import Table
from collections import defaultdict
from astropy.coordinates import SkyCoord

def parse_tract_corners(file_path):
    """
    Parse the tract corner information from the file.
    """
    tracts = defaultdict(lambda: {"corners": []})
    with open(file_path, 'r') as file:
        lines = file.readlines()
        current_tract = None

        for line in lines:
            if "Tract:" in line and "Corner" in line and "Patch" in line:
                parts = line.split()
                current_tract = parts[1]
                current_patch = parts[3]
                current_tract = f"{current_tract}-{current_patch}"
                ra = float(parts[-3].strip("(").strip(","))
                dec = float(parts[-1].strip(")"))
                tracts[current_tract]["corners"].append((ra, dec))
    return tracts

def read_all_tracts(file_dir):
    """
    Read all text files in the directory and extract tract information.
    """
    all_tracts = {}
    txt_files = glob.glob(f"{file_dir}/*.txt")
    for file_path in txt_files:
        tracts = parse_tract_corners(file_path)
        all_tracts.update(tracts)
    return all_tracts

def is_point_in_tract(ra_array, dec_array, tract_corners):
    """
    Vectorized version to check if points (RA, Dec) are within a tract defined by corners.
    Assumes corners are in order (clockwise or counterclockwise).
    """
    def points_in_polygon(x, y, polygon):
        """
        Vectorized ray-casting algorithm to check if points (x, y) are inside the polygon.
        """
        polygon = np.array(polygon)
        x_coords, y_coords = polygon[:, 0], polygon[:, 1]
        n = len(polygon)

        # Initialize an array to store the "inside" status for all points
        inside = np.zeros(len(x), dtype=bool)

        for i in range(n):
            x1, y1 = x_coords[i], y_coords[i]
            x2, y2 = x_coords[(i + 1) % n], y_coords[(i + 1) % n]

            # Check if the y-coordinates of points fall between the edge's y-coordinates
            cond = (y1 <= y) & (y < y2) | (y2 <= y) & (y < y1)

            # Compute intersection points of the edge with horizontal lines at y
            xinters = (y - y1) * (x2 - x1) / (y2 - y1) + x1

            # Check if the x-coordinate of points is to the left of the intersection
            cross = x < xinters

            # Toggle the inside state for points meeting the condition
            inside ^= cond & cross

        return inside

    return points_in_polygon(ra_array, dec_array, tract_corners)

def create_mask(ra_list, dec_list, all_tracts):
    """
    Vectorized function to create a mask for a list of RA and Dec coordinates.
    """
    ra_array = np.array(ra_list)
    dec_array = np.array(dec_list)
    mask = np.zeros(len(ra_array), dtype=bool)

    for tract_id, tract_info in all_tracts.items():
        corners = tract_info["corners"]
        # Vectorized check for all points within the current tract
        mask |= is_point_in_tract(ra_array, dec_array, corners)

    return mask


def filter_catalog_by_footprint(table, ra_column, dec_column, tract_dir, radec_units="deg"):
    """
    Filter the input catalog by the HSC PDR2 Wide footprint.
    """
    print("Initial number of objects: ", len(table))
    ra, dec = table[ra_column], table[dec_column]
    if radec_units == "deg":
        pass
    elif radec_units == "hourangle":
        coord = SkyCoord(ra=ra, dec=dec, unit=("hourangle", "deg"))
        ra, dec = coord.ra.deg, coord.dec.deg

    all_tracts = read_all_tracts(tract_dir)
    mask = create_mask(ra, dec, all_tracts)
    table = table[mask]
    print("Final number of objects: ", len(table))
    return table

