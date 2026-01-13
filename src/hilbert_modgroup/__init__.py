import logging

from .version import version as __version__

logging.basicConfig(level=logging.ERROR, format='%(message)s')
logging.captureWarnings(True)
