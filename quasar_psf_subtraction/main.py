from .make_galfit import make_run_galfit
from .query_hsc import query_hsc
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

import numpy as np
import os
import astropy.units as u
from astropy.io import fits
from .utils.measurements import measure_background, make_central_mask, make_mask_sources
from .utils.plot import plot_target_results, plot_all_targets
from astropy.table import Table
from unagi import hsc
from astropy.coordinates import SkyCoord
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")

def fit_source(targname,
                data_path,
                ra,
                dec,
                archive,
                filters=None,
                cutout_size=3*u.arcsec,
                substract_sky=True,
                r_mask_central_source=10,
                n_sigma_mask_sources=3,
                only_central_mask=False,
                output_download_basename="hsc_data",
                run_single_band=True,
                run_multi_band=True,
                pixscale=None,
                plot=True):

    # get the data and PSFs
    coord = SkyCoord(ra, dec, unit=(u.deg, u.deg))

    cutouts, psfs = query_hsc(targname, coord, archive,
                              filters=filters, cutout_size=cutout_size,
                              output_basename=output_download_basename)
    if cutouts is None:
        print("no data")
        return

    filters = list([hdu[0].header["FILTER"] for hdu in cutouts])
    output_dir = os.path.join(data_path, targname)
    os.makedirs(output_dir, exist_ok=True)
    if not len(cutouts):
        os.rmdir(output_dir)
        return


    filenames_data = []
    filenames_psf = []
    filenames_mask = []
    filters_list = []
    for i, image_hdul, psf_hdul in zip(range(len(cutouts)), cutouts, psfs):
        filter = image_hdul[0].header["FILTER"]
        filename_data = f'{output_dir}/cutout_{filter}.fits'
        filename_data = os.path.join(os.getcwd(), filename_data)
        filename_psf = f'{output_dir}/psf_{filter}.fits'
        filename_psf = os.path.join(os.getcwd(), filename_psf)
        filename_mask = f'{output_dir}/mask_{filter}.fits'
        filename_mask = os.path.join(os.getcwd(), filename_mask)

        from astropy.stats import sigma_clipped_stats
        data = image_hdul[1].data

        background = measure_background(data, box_size=20, filter_size=3, n_sigma_clip=3)

        if substract_sky:
            data -= background

        # 1 where to fit 0 where is masked
        mask_central = make_central_mask(data, radius=r_mask_central_source)
        mask_sources = make_mask_sources(data, n_sigma=n_sigma_mask_sources)

        if only_central_mask:
            mask = mask_central
        else:
            mask = np.logical_or(mask_central, mask_sources)

        mask = np.logical_not(mask)

        fits.writeto(filename_data, image_hdul[1].data, image_hdul[1].header, overwrite=True)
        fits.writeto(filename_psf, psf_hdul[0].data, psf_hdul[0].header, overwrite=True)
        fits.writeto(filename_mask, mask.astype(int), overwrite=True)
        filenames_data.append(filename_data)
        filenames_psf.append(filename_psf)
        filenames_mask.append(filename_mask)
        filters_list.append(filter)

    make_run_galfit(output_dir, filenames_data, filenames_psf,
                    filenames_mask, filters_list, cutouts[0][1].data.shape,
                    single_band=run_single_band, multi_band=run_multi_band)


    if plot:
        fig, axs = plot_target_results(targname, data_path, ra=ra, dec=dec, pixscale=pixscale)
        fig.savefig(f"{output_dir}/results.pdf", bbox_inches='tight', dpi=300)


def from_config_file(config_file):
    # table in csv format targname ra hh:mm:ss, dec hh:mm:ss
    import yaml

    # Load the YAML configuration file
    with open(config_file, "r") as file:
        config = yaml.safe_load(file)

    n_processes = config["n_processes"]
    fit_first_n = config["fit_first"]
    if fit_first_n is not None:
        fit_first_n = int(fit_first_n)

    catalog = Table.read(config["catalog_filename"])
    catalog = catalog[:fit_first_n]

    ra_column, dec_column, targname_column = config["ra_column"], config["dec_column"], config["targname_column"]
    ra_unit, dec_unit = config["ra_unit"], config["dec_unit"]
    ra_list, dec_list, targname_list = catalog[ra_column], catalog[dec_column], catalog[targname_column]
    output_dir = config["output_dir"]

    data_release = config["hsc_query"]["data_release"]
    output_hsc = config["hsc_query"]["output_hsc"]
    rerun = config["hsc_query"]["rerun"]
    config_file_credentials = config["hsc_query"]["config_file_credentials"]
    cutout_size = float(config["hsc_query"]["cutout_size"])*u.arcsec
    pix_scale = float(config["hsc_query"]["pixel_scale"])

    archive = hsc.Hsc(dr=data_release, rerun=rerun, config_file=config_file_credentials)

    substract_sky = config["galfit"]["substract_sky"]
    r_mask_central_source = config["galfit"]["r_mask_central_source"]
    only_mask_central_source = config["galfit"]["only_mask_central_source"]
    n_sigma_mask_sources = config["galfit"]["n_sigma_mask_sources"]
    run_single_band = config["galfit"]["run_single_band"]
    run_multi_band = config["galfit"]["run_multi_band"]
    plot = config["galfit"]["plot"]


    def run_mp(i):
        ra, dec, targname = ra_list[i], dec_list[i], targname_list[i]
        coord = SkyCoord(ra, dec, unit=(u.Unit(ra_unit), u.Unit(dec_unit)))
        ra_deg, dec_deg = coord.ra.deg, coord.dec.deg
        fit_source(targname, output_dir, ra_deg, dec_deg, archive,
               filters=None, cutout_size=cutout_size,
               output_download_basename=output_hsc, plot=plot,
               pixscale=pix_scale,
                substract_sky=substract_sky,
                r_mask_central_source=r_mask_central_source,
                n_sigma_mask_sources=n_sigma_mask_sources,
                only_central_mask=only_mask_central_source,
                run_single_band=run_single_band,
                run_multi_band=run_multi_band,
                )



    from multiprocessing.dummy import Pool as ThreadPool

    # copy config file inside output path
    os.makedirs(output_dir, exist_ok=True)
    import shutil
    shutil.copy(config_file, os.path.join(output_dir, os.path.basename(config_file)))

    p = ThreadPool(n_processes)
    p.map(run_mp, range(len(ra_list)))
    p.close()
    p.join()

    plot_all_targets(output_dir, f"{output_dir}/fits_results.pdf")