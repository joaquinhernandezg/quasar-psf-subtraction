import argparse
from quasar_psf_subtraction.main import from_config_file
import matplotlib
matplotlib.use("Agg")

def main():
    parser = argparse.ArgumentParser(description="A script in my Python module.")
    parser.add_argument('config_file', type=str, help="Path to the configuration file")

    args = parser.parse_args()

    config_file = args.config_file
    from_config_file(config_file)

if __name__ == "__main__":
    main()
