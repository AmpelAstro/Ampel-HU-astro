# copied from https://github.com/nltk/nltk/blob/8c233dc585b91c7a0c58f96a9d99244a379740d5/nltk/util.py
#
# Natural Language Toolkit: Utility functions
#
# SPDX-SnippetCopyrightText: 2001-2023 NLTK Project
# SPDX-Author: Steven Bird <stevenbird1@gmail.com>
# SPDX-Author: Eric Kafe <kafe.eric@gmail.com> (acyclic closures)

# ruff: noqa


def flatten(*args):
    """
    Flatten a list.

        >>> from nltk.util import flatten
        >>> flatten(1, 2, ['b', 'a' , ['c', 'd']], 3)
        [1, 2, 'b', 'a', 'c', 'd', 3]

    :param args: items and lists to be combined into a single list
    :rtype: list
    """

    x = []
    for l in args:
        if not isinstance(l, (list, tuple)):
            l = [l]
        for item in l:
            if isinstance(item, (list, tuple)):
                x.extend(flatten(item))
            else:
                x.append(item)
    return x
