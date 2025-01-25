import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.io import fits
from astropy.visualization import simple_norm
import glob
from scipy.ndimage import gaussian_filter
from matplotlib.backends.backend_pdf import PdfPages
from astropy.coordinates import SkyCoord


def plot_target_results(targname, data_path, ra=None, dec=None, pixscale=None):
    target_path = os.path.join(data_path, targname)

    output_galfit_multiband = os.path.join(target_path, "galfit_multiband/output.fits")
    output_galfit_single_band  = glob.glob(os.path.join(target_path, "galfit_singleband/output_*.fits"))
    n_bands = len(output_galfit_single_band )
    try:
        output_multiband_hdul = fits.open(output_galfit_multiband)
    except Exception as e:
        print(f"Error with {targname}")
        return None, None

    n_cols = 3
    n_rows = n_bands
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(3*n_cols, 3*n_rows))

    for i in range(n_bands):
        band_name = output_galfit_single_band[i].split("/")[-1].split("_")[-1].split(".")[0]
        try:
            output_single_band_hdul = fits.open(output_galfit_single_band[i])
        except Exception as e:
            continue

        data = output_single_band_hdul[1].data
        model_single_band = output_single_band_hdul[2].data
        res_single_band = output_single_band_hdul[3].data

        model_multi_band = output_multiband_hdul[(1+n_bands*1+i)].data
        residual_multi_band = output_multiband_hdul[(1+n_bands*2+i)].data

        norm = simple_norm(residual_multi_band, min_percent=1, max_percent=99)


        center_y, center_x = data.shape[0]//2, data.shape[1]//2



        if n_bands == 1:
            axs_row = axs
        else:
            axs_row = axs[i]
        axs_row[0].imshow(data, cmap="viridis", origin="lower", norm=norm)
        # text in 0.2, 0.8
        axs_row[0].text(0.2, 0.8, band_name, transform=axs_row[0].transAxes, color="w", fontsize=20)


        #axs_row[1].imshow(model_single_band, cmap="viridis", origin="lower")

        axs_row[1].imshow(res_single_band, cmap="viridis", origin="lower", norm=norm)
        data_ = gaussian_filter(model_single_band, 0.2)
        xx, yy = np.meshgrid(np.arange(data_.shape[1]), np.arange(data_.shape[0]))
        max_val = np.max(data_)
        median_val = np.median(data_)
        percents = []
        levels = [median_val + (max_val - median_val) * p / 100 for p in percents]
        axs_row[1].contour(xx, yy, data_, levels=levels, colors="k", alpha=0.5)
        axs_row[0].contour(xx, yy, data_, levels=levels, colors="k", alpha=0.5)

        # draw a ruler indicating 5 arcseconds if pixscale is not none
        if pixscale is not None:
            scale = 5 / pixscale
            # plot it in frac position x 0.5 y 0.8 with transaxes
            shape_y, shape_x = data.shape
            x0, x1 = shape_x//2, shape_x//2 + scale
            y0, y1 = shape_y*0.8, shape_y*0.8
            axs_row[1].plot([x0, x1], [y0, y1], "r-", color="w", lw=2)
            axs_row[1].text(0.5, 0.85, "5 arcsec", transform=axs_row[1].transAxes, color="w", fontsize=12)



        axs_row[2].imshow(residual_multi_band, cmap="viridis", origin="lower", norm=norm)
        data_ = gaussian_filter(model_multi_band, 0.2)
        xx, yy = np.meshgrid(np.arange(data_.shape[1]), np.arange(data_.shape[0]))
        max_val = np.max(data_)
        median_val = np.median(data_)
        percents = []
        levels = [median_val + (max_val - median_val) * p / 100 for p in percents]
        axs_row[2].contour(xx, yy, data_, levels=levels, colors="r", alpha=0.5)
        axs_row[0].contour(xx, yy, data_, levels=levels, colors="r", alpha=0.5)


        for ax in axs_row[:3]:
            ax.set_axis_off()
        if i==0:
            axs_row[0].set_title("Data")
            axs_row[1].set_title("Residual single band")
            axs_row[2].set_title("Residual multi band")

    title = targname

    if ra is not None and dec is not None:
        coord = SkyCoord(ra=ra, dec=dec, unit="deg")
        ra_hhmmss = coord.ra.to_string(unit="hour", sep=":", precision=2)
        dec_ddmmss = coord.dec.to_string(unit="deg", sep=":", precision=2)
        title += f"\nJ{ra_hhmmss} {dec_ddmmss}"

    fig.suptitle(title, fontsize=12)
    return fig, axs

def plot_all_targets(data_path, pdf_name):
    pp = PdfPages(pdf_name)
    for targname in os.listdir(data_path):
        try:
            fig, axs = plot_target_results(targname, data_path)
            pp.savefig(fig)
        except Exception as e:
            print(f"Error with {targname}")


    pp.close()