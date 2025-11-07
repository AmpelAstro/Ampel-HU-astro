#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/aiotns
# License:             BSD-3-Clause
# Author:              Jakob van Santen <jakob.van.santen@desy.de>
# Date:                04.11.2019
# Last Modified Date:  04.11.2019
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import math


class TNSName:
    """
    TNS objects are numbered with a variable-length, 1-indexed letter code,
    e.g. 2019a represents (2019,1), while 2019aa represents (2019,27)
    """

    __slots__ = ("number", "year")
    BASE = 26

    def __init__(self, year: int, number: int):
        self.year = year
        self.number = number

    @classmethod
    def from_str(cls, name):
        assert 0 < len(name) - 4 < 11
        year = int(name[:4])
        number = 0
        for i, c in enumerate(reversed(name[4:])):
            number += (ord(c.lower()) - ord("a") + 1) * (cls.BASE**i)
        return cls(year, number)

    @classmethod
    def from_index(cls, packed):
        return cls((packed >> 56) + 1970, packed & ((1 << 56) - 1))

    def __int__(self):
        """
        Pack year into top 8 bits, intra-year number into bottom 56
        """
        return ((self.year - 1970) << 56) + self.number

    def __str__(self):
        name = str(self.year)
        if self.number <= 0:
            return name
        # find number of digits to emit
        digits = math.ceil(math.log(self.number) / math.log(self.BASE))
        for i in range(digits - 1, -1, -1):
            unit = self.BASE**i
            name += chr(ord("a") + (self.number % (unit * self.BASE)) // unit - 1)
        return name

    def __repr__(self):
        return f"TNSName({self.year},{self.number})"
