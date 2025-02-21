import astropy.units as u
from astropy.coordinates import SkyCoord

from .query_hsc import query_hsc


class Catalog:

    def __init__(self, catalog_file, catalog_format="fits"):
        self.catalog_file = catalog_file
        self.catalog = self._read_catalog()

    def _read_catalog(self):
        pass

    


class Target:

    def __init__(self, ra, dec, ra_unit="deg", dec_unit="deg", targetid=None, target_prefix="SDSS"):
        self.ra = ra * u.Unit(ra_unit)
        self.dec = dec * u.Unit(dec_unit)
        self.coord = SkyCoord(ra=self.ra, dec=self.dec)
        self.targetid = self.construct_target_id() if targetid is None else targetid

        self.cutouts = None
        self.psfs = None
        

    def get_cutouts(self, archive, filters=None, cutout_size=None, output_download_basename=None):
        
        if self.cutouts is None:
            self.cutouts, self.psfs = query_hsc(self.targetid, self.coord, archive,
                              filters=filters, cutout_size=cutout_size,
                              output_basename=output_download_basename)
        return self.cutouts, self.psfs


    def construct_target_id(self):
        ra_hhmmss = self.coord.ra.to_string(unit="hour", sep=":", precision=2)
        dec_ddmmss = self.coord.dec.to_string(unit="deg", sep=":", precision=2)
        return f"{self.target_prefix} J{ra_hhmmss}{dec_ddmmss}"
    
    def cds_url(self):
        ra_hh_mm_ss = self.coord.ra.to_string(unit="hour", sep="", precision=2, pad=True)
        dec_dd_mm_ss = self.coord.dec.to_string(unit="deg", sep="", precision=1, pad=True, alwayssign=True)
        url = f"http://cdsportal.u-strasbg.fr/?target=SDSS%20J{ra_hh_mm_ss}%20{dec_dd_mm_ss}"
        return url

