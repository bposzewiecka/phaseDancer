# pylint: disable=no-member

import logging


def set_logger():

    logger = logging.getLogger("phaseDancer")
    logger.setLevel(logging.DEBUG)

    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)

    logger.addHandler(handler)

    return logger


pLogger = set_logger()
