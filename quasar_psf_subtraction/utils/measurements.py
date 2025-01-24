from astropy.io import fits
import numpy as np
import os
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.background import Background2D, MedianBackground


def get_guess_mag(data, x_center=None, y_center=None, radius=5):
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    data -= median
    if x_center is None:
        x_center = data.shape[1] // 2
    if y_center is None:
        y_center = data.shape[0] // 2

    xx, yy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
    r = np.sqrt((xx - x_center)**2 + (yy - y_center)**2)
    mask = r < radius
    data[~mask] = np.nan
    guess_mag = -2.5 * np.log10(np.abs(np.nansum(data))+1e-40)
    return guess_mag


def get_guess_mags(filenames, x_center=None, y_center=None, radius=5):
    guess_mags = []
    for filename in filenames:
        data = fits.open(filename)[0].data
        guess_mag = get_guess_mag(data, x_center=x_center, y_center=y_center, radius=radius)
        guess_mags.append(guess_mag)
    return guess_mags

def measure_background(data, box_size=20, filter_size=3, n_sigma_clip=3):
    sigma_clip = SigmaClip(sigma=n_sigma_clip)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (box_size, box_size), filter_size=(filter_size, filter_size),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    return bkg.background


def make_central_mask(data, x_center=None, y_center=None, radius=5):
    if x_center is None:
        x_center = data.shape[1] // 2
    if y_center is None:
        y_center = data.shape[0] // 2

    xx, yy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
    r = np.sqrt((xx - x_center)**2 + (yy - y_center)**2)
    mask = r < radius
    return mask

def make_mask_sources(data, n_sigma=5):
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    mask = data < median + n_sigma * std
    return mask


