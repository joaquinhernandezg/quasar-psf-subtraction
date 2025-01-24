from .make_galfit import make_run_galfit
from .query_hsc import query_hsc
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

import numpy as np
import os
import astropy.units as u
from astropy.io import fits
from .utils.measurements import measure_background, make_central_mask, make_mask_sources
from .utils.plot import plot_target_results
from astropy.table import Table

def fit_source(targname,
                data_path,
                ra_hh_mm_ss,
                dec_hh_mm_ss,
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
                plot=True):

    # get the data and PSFs
    cutouts, psfs = query_hsc(targname, ra_hh_mm_ss, dec_hh_mm_ss, archive,
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
        fig, axs = plot_target_results(targname, data_path)
        fig.savefig(f"{output_dir}/results.pdf", bbox_inches='tight', dpi=300)


def main(config_file):
    # table in csv format targname ra hh:mm:ss, dec hh:mm:ss



    import yaml

    # Load the YAML configuration file
    with open(config_file, "r") as file:
        config = yaml.safe_load(file)
