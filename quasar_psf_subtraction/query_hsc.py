import os

import astropy.units as u
from astropy.coordinates import SkyCoord

from unagi.task import hsc_cutout
from unagi.task import hsc_psf

from astropy.coordinates import SkyCoord
import astropy.units as u

from unagi import task



def query_hsc(targname, ra_hh_mm_ss, dec_hh_mm_ss, archive, filters=None,
              cutout_size=3*u.arcsec, output_basename="hsc_data"):

    #filters = 'gri'
    coord = SkyCoord(ra_hh_mm_ss, dec_hh_mm_ss, unit=(u.hourangle, u.deg))

    if filters is None:
        coverage = task.hsc_check_coverage(coord, archive=archive, verbose=True)
        filters = set(coverage["filter01"])
        filters = list(map(lambda x: x.split("-")[-1].lower(), filters))
    if not len(filters):
        return None, None
    # Output dir
    output_dir = os.path.join(output_basename, targname)
    os.makedirs(output_dir, exist_ok=True)
    s_ang = cutout_size

    cutouts_hduls = []
    psf_hduls = []

    for filter in filters:
        try:
            cutouts = hsc_cutout(coord, cutout_size=s_ang, filters=filter, archive=archive,
                                    use_saved=True, output_dir=output_dir, verbose=True,
                                    save_output=True)
            psfs = hsc_psf(coord, filters=filter, archive=archive,
                           output_dir=output_dir, use_saved=True, verbose=True)
        except Exception as e:
            print(f"Error with {filter}: {e}")
            continue
        cutouts_hduls.append(cutouts)
        psf_hduls.append(psfs)

    return cutouts_hduls, psf_hduls
