
import copy

def info_as_debug(logger):
	quiet = copy.copy(logger)
	quiet.info = quiet.debug
	return quiet
