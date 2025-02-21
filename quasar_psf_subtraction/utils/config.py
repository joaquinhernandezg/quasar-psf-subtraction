

class Config:

    def __init__(self, config_file):
        self.config_file = config_file
        self.config = self.read_config()

    @property
    def ra_deg(self):
        ra_column = self.config['ra_column']
        
        return self.config[ra_column]

    