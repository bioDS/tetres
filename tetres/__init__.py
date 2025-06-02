import logging

if not logging.getLogger().hasHandlers():
    logging.basicConfig(
        format='%(asctime)s [%(name)s:%(lineno)d] %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S',
        level=logging.INFO,
    )
