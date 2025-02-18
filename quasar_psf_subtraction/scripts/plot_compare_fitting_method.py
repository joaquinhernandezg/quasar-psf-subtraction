import argparse
import matplotlib
import os
from collections import defaultdict
import glob
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
matplotlib.use("Agg")


def main():
    parser = argparse.ArgumentParser(description="Plot the comparison of different fitting methods")
    parser.add_argument("--pattern", type=str, default="*",
                        help="Pattern to search for the directories")
    args = parser.parse_args()
    pattern = args.pattern


    current_dir = os.getcwd()
    # read all directories in current directory, not files matching pattern
    # use glob
    directories = glob.glob(pattern)

    subdirs = {}
    all_targets = []

    for directory in directories:
        # list subdirs
        subdirs[directory] = [d for d in os.listdir(os.path.abspath(directory)) if os.path.isdir(os.path.abspath(os.path.join(directory, d)))]
        all_targets.extend(subdirs[directory])

    all_targets = list(set(all_targets))
    all_dirs = list(set(subdirs.keys()))


    for target in all_targets:

        available_bands = []
        for directory in all_dirs:
            cutout_files = glob.glob(os.path.join(directory, target, "cutout*.fits"))
            for cutout_file in cutout_files:
                band = cutout_file.split("_")[-1].split(".")[0]
                available_bands.append(band)
        available_bands = sorted(list(set(available_bands)))

        n_bands, n_methods = len(available_bands), len(all_dirs)
        fig, axes = plt.subplots(n_bands, n_methods+1, figsize=((n_methods+1)*5, n_bands*5))

        for i in range(n_bands):
            band = available_bands[i]
            axes[i, 0].set_title(band)

            for j in range(1, n_methods+1):
                method = all_dirs[j-1]
                title = method
                every = 20
                title = '\n'.join(title[i:i+every] for i in range(0, len(title), every))
                axes[i, j].set_title(title)
                print(target, band, method)

                cutout_file = os.path.join(method, target, "cutout_{}.fits".format(band))


                residual_file = os.path.join(method, target, "galfit_singleband", "output_{}.fits".format(band))

                try:
                    cutout_data = fits.getdata(cutout_file, ext=0)
                    mean, median, std = sigma_clipped_stats(cutout_data)
                    x_center, y_center = cutout_data.shape[1]//2, cutout_data.shape[0]//2
                    size = 40#pix


                    axes[i, 0].imshow(cutout_data, origin="lower", cmap="viridis", vmin=mean-5*std, vmax=mean+5*std)
                    axes[i, 0].axis("off")
                    axes[i, 0].set_xlim(x_center-size, x_center+size)
                    axes[i, 0].set_ylim(y_center-size, y_center+size)

                except Exception as e:
                    print(e)
                    continue

                try:
                    residual_data = fits.getdata(residual_file, ext=3)
                    # plot data in first column plot residual in column j

                    axes[i, j].imshow(residual_data, origin="lower", cmap="coolwarm", vmin=mean-5*std, vmax=mean+5*std)
                    axes[i, j].axis("off")
                    axes[i, j].set_xlim(x_center-size, x_center+size)
                    axes[i, j].set_ylim(y_center-size, y_center+size)
                except Exception as e:
                    print(e)
                    continue

        plt.savefig(os.path.join(current_dir, "plots", "{}_comparison.png".format(target)))










