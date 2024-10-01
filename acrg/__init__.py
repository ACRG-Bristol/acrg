import logging

from rich.logging import RichHandler as _RichHandler

# Start module level logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logging.captureWarnings(capture=True)

# Create console handler - set to WARNING (lower level)
consoleHandler = _RichHandler()
consoleFormatter = logging.Formatter("%(levelname)s:%(name)s:%(message)s", datefmt="%Y-%m-%dT%H:%M:%S%z")
consoleHandler.setFormatter(consoleFormatter)
consoleHandler.setLevel(logging.INFO)
logger.addHandler(consoleHandler)
