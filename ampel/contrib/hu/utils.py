import copy


def info_as_debug(logger):
    quiet = copy.copy(logger)
    quiet.info = quiet.debug
    return quiet


def _flatten(*args):
    for l in args:
        if isinstance(l, list | tuple):
            for item in l:
                yield from _flatten(item)
        else:
            yield l


def flatten(*args):
    "Flatten a list; similar to nltk.util.flatten"
    return list(_flatten(*args))
