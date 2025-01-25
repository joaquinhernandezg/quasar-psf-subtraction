import numpy as np
import glob

def parse_tract_corners(file_path):
    """
    Parse the tract corner information from the file.
    """
    tracts = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        current_tract = None

        for line in lines:
            if "Tract:" in line and "Corner" in line and "Patch" not in line:
                parts = line.split()
                current_tract = int(parts[1])
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

def is_point_in_tract(ra, dec, tract_corners):
    """
    Check if a point (ra, dec) is within the tract defined by corners.
    Assumes corners are in order (clockwise or counterclockwise).
    """
    def point_in_polygon(x, y, polygon):
        # Ray-casting algorithm to check if point is inside polygon
        n = len(polygon)
        inside = False
        px, py = x, y
        for i in range(n):
            x1, y1 = polygon[i]
            x2, y2 = polygon[(i + 1) % n]
            if min(y1, y2) < py <= max(y1, y2) and px <= max(x1, x2):
                if y1 != y2:
                    xinters = (py - y1) * (x2 - x1) / (y2 - y1) + x1
                if x1 == x2 or px <= xinters:
                    inside = not inside
        return inside

    return point_in_polygon(ra, dec, tract_corners)

def create_mask(ra_list, dec_list, all_tracts):
    """
    Create a mask for a list of RA and Dec coordinates.
    """
    mask = np.zeros(len(ra_list), dtype=bool)
    for i, (ra, dec) in enumerate(zip(ra_list, dec_list)):
        for tract_id, tract_info in all_tracts.items():
            if is_point_in_tract(ra, dec, tract_info["corners"]):
                mask[i] = True
                break
    return mask