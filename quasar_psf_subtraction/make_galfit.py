import os
from string import Template
from .utils.measurements import get_guess_mags, measure_background, make_central_mask, make_mask_sources
from .utils.plot import plot_target_results
from .query_hsc import query_hsc
import numpy as np
from astropy.io import fits
import astropy.units as u

WAVES = {"u": 0.35*1000, "g": 0.48*1000, "r": 0.62*1000, "i": 0.75*1000, "z": 0.89*1000, "i2": 0.75*1000,
             "r2": 0.62*1000, "z2": 0.89*1000, "g2": 0.48*1000, "y": 1.0*1000}


def fill_template(template_file, params):
    template_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), "templates", template_file)
    with open(template_file, "r") as f:
        template = Template(f.read())
        content = template.substitute(params)
    return content


def make_run_galfit(output_dir, filenames_data, filenames_psf, filenames_mask, bands, shape,
                    single_band=True, multi_band=True, substract_sky=True):
    #full path
    output_dir = os.path.join(os.getcwd(), output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    if multi_band:
        multi_band_dir = os.path.join(output_dir, 'galfit_multiband')

        galfit_filename = make_galfit_multiband_file(multi_band_dir, filenames_data, filenames_psf, filenames_mask, bands, shape)
        os.system(f"galfitm {galfit_filename}")

    if single_band:
        single_band_dir = os.path.join(output_dir, 'galfit_singleband')
        for filename_data, filename_psf, filename_mask, band in zip(filenames_data, filenames_psf, filenames_mask, bands):
            config_filename = f"config_{band}.galfit"
            output_filename = f"output_{band}.fits"
            galfit_filename = make_galfit_multiband_file(single_band_dir, [filename_data], [filename_psf], [filename_mask], [band], shape,
                                                         config_filename=config_filename, output_filename=output_filename)
            os.system(f"galfitm {galfit_filename}")
            print(f"Done with {band}")



def make_galfit_multiband_file(output_dir,
                              filenames_data,
                              filenames_psf,
                              filenames_mask,
                              bands,
                              shape,
                              config_filename="config.galfit",
                              output_filename="output.fits"):
    #full path
    #output_dir = os.path.join(os.getcwd(), output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #sort band by wavelength
    #print(bands)
    #bands = sorted(bands, key=lambda x: WAVES[x])

    waves = [str(WAVES[band]) for band in bands]

    output_filename = os.path.join(output_dir, output_filename)

    galfit_filename = os.path.join(output_dir, config_filename)

    guess_mags = get_guess_mags(filenames_data)

    pixscale = 0.262
    shape_y, shape_x = shape
    center_x, center_y = shape_x / 2, shape_y / 2
    size_fit = 10
    sigma_image = "none"
    psf_fine_sampling_factor = 1
    constraints_file = "none"
    convolution_size = 20
    zeropoint = 0


    galfit_params = dict(
        input_image = ','.join(filenames_data),
        bands_names = ','.join(iter(bands)),
        bands_waves = ','.join(waves),
        output_block = output_filename,
        sigma_image = sigma_image,
        psf_image = ','.join(filenames_psf),
        psf_fine_sampling_factor = psf_fine_sampling_factor,
        bad_pixel_mask = ','.join(filenames_mask),
        constraints_file = constraints_file,
        x_min = 0,
        x_max = shape_x,
        y_min = 0,
        y_max = shape_y,
        conv_size_x = convolution_size,
        conv_size_y = convolution_size,
        phot_zeropoint = zeropoint,
        plate_scale_x = pixscale,
        plate_scale_y = pixscale,
    )

    psf_params = dict(
        pos_x = center_x,
        pos_y = center_y,
        pos_x_free = 1,
        pos_y_free = 1,
        magnitude =  ",".join(list(map(str, guess_mags))),
        magnitude_free = ",".join("1" for i in bands),
        skip_model_output="0"
    )

    sky_params = dict(
        sky_level = ",".join(["0" for i in bands]),
        sky_level_free = ",".join(["1" for i in bands]),
        slope_x_sky = ",".join(["0" for i in bands]),
        slope_y_sky = ",".join(["0" for i in bands]),
        slope_x_sky_free = ",".join(["1" for i in bands]),
        slope_y_sky_free = ",".join(["1" for i in bands]),
        skip_model_output="0",
    )

    galfit_content = fill_template("galfit_multiband.template", galfit_params)
    psf_content = fill_template("psf.template", psf_params)
    sky_content = fill_template("sky.template", sky_params)

    with open(galfit_filename, "w") as f:
        f.write(galfit_content)
        f.write(psf_content)
        f.write(sky_content)

    return galfit_filename



if __name__ == "__main__":
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from unagi import hsc

    ra_targ = "01:25:55.11"
    dec_targ ="-01:29:25.00"
    coord = SkyCoord(ra_targ, dec_targ, unit=(u.hourangle, u.deg))

    targname = "Test"
    data_path = "Test"
    archive_pdr2_wide = hsc.Hsc(dr='pdr2', rerun='pdr2_wide',config_file="/data/estudiantes/jhernandez/PasquierPSF/config_file.txt")

    from unagi import task

    #coverage_pdr2_wide = task.hsc_check_coverage(coord, archive=archive_pdr2_wide, verbose=True)
    #filters = set(coverage_pdr2_wide["filter01"])
    #filters = list(map(lambda x: x.split("-")[-1].lower(), filters))
    #query_hsc(targname, ra_targ, dec_targ, archive_pdr2_wide, filters=filters, cutout_size=3*u.arcsec, output_basename="hsc_data")

    #coverage_pdr3_wide = task.hsc_check_coverage(coord, archive=archive_pdr3_wide, verbose=True)
    #print(coverage_pdr2_wide)
    #print(coverage_pdr2_wide)

    ra_targ = "01:25:55.11"
    dec_targ ="-01:29:25.00"
    coord = SkyCoord(ra_targ, dec_targ, unit=(u.hourangle, u.deg))

    #coverage_pdr2_wide = task.hsc_check_coverage(coord, archive=archive_pdr2_wide, verbose=True)
    #filters = set(coverage_pdr2_wide["filter01"])



    #coverage_pdr3_wide = task.hsc_check_coverage(coord, archive=archive_pdr3_wide, verbose=True)
    #print(coverage_pdr2_wide)
    data_path = "/data/estudiantes/jhernandez/PasquierPSF/Test"
    fit_source(targname, data_path, ra_targ, dec_targ, archive_pdr2_wide,
               filters=None, cutout_size=10*u.arcsec,
               output_download_basename="hsc_data", plot=True,
               )