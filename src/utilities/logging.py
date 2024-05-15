import logging
logger = logging.getLogger('dysfunction_logger')
logger.setLevel(logging.DEBUG)
hd = logging.FileHandler('debuglog.log')
hd.setLevel(logging.DEBUG)
cons = logging.StreamHandler()
cons.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
hd.setFormatter(formatter)
cons.setFormatter(formatter)
logger.addHandler(hd)
logger.addHandler(cons)

