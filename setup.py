from setuptools import setup, find_packages

setup(
    name="quasar_psf_subtraction",                  # Package name
    version="0.1.0",                    # Package version
    description="A package to fit and substract PSF from Subaru Data",  # Short description
    long_description=open("README.md").read(),  # Long description (README file)
    long_description_content_type="text/markdown",
    author="Joaquin A. Hernández",
    author_email="joaquin.hernandezg@uc.cl",
    url="https://github.com/joaquinhernandezg/quasar-psf-subtraction",  # Project URL
    license="MIT",                      # License
    packages=find_packages(),           # Automatically find all sub-packages
    install_requires=[
        # List dependencies here
        "numpy",
        "scipy",
        "astropy",
        "matplotlib",
        "unagi",

    ],

    python_requires=">=3.10.16",            # Supported Python versions
)