import logging
# Setting up a logger with a tree structure for the MPL module
logger = logging.getLogger(__name__)

# If we are loading the module without the rest of madgraph, setup a basic logger
if len(logger.handlers) == 0:
    logging.basicConfig()  # Really basic setup here

# The MPL module loggers all inherit from this one and will respect this level
logger.setLevel(logging.INFO)

from mpl_dispatch import MPL_dispatcher
from mpl_implementation import mpl_implementations

G = MPL_dispatcher(mpl_implementations)
