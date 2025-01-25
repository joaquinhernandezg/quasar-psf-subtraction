import argparse
from quasar_psf_subtraction.main import from_config_file
import os
from string import Template

def main():
    template_file = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "templates", "config.yaml")
    # copy file to current path
    with open(template_file, "r") as f:
        template = f.read()
    with open("config.yaml", "w") as f:
        f.write(template)
    print("config.yaml file created")



if __name__ == "__main__":
    main()
