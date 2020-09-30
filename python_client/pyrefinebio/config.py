# WIP
import os
from pathlib import Path

import yaml


class Config:

    _instance = None

    def __new__(cls, config_file=None):
        if not cls._instance:
            cls._instance = super(Config, cls).__new__(cls)

            cls.config_file = (
                os.getenv("CONFIG_FILE") or config_file or str(Path.home()) + "/.refinebio.yaml"
            )

        return cls._instance
