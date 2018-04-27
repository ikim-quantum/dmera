# This is part of the DMERA project(https://github.com/ikim-quantum/dmera).
# Copyright (C) 2018 Isaac H. Kim.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from qubit import Qubit


class Coordinate(tuple):
    def __init__(self, *args):
        self = (*args)

    def __mul__(self, other):
        if isinstance(other, int):
            return ([n * other for n in self])
        elif isinstance(other, tuple):
            return ([n * factor for n, factor in zip(self, other)])

    def __imul__(self, other):
        if isinstance(other, int):
            self = self * other
        elif isinstance(other, tuple):
            self = self * other

    def __add__(self, other):
        return ([i + j for i, j in zip(self, other)])

    def __iadd__(self, other):
        self = self + other

    def __sub__(self, other):
        return ([i - j for i, j in zip(self, other)])

    def __isub__(self, other):
        self = self - other


class Lattice():

    qubits = []
    d = 0
    size = 0

    def __init__(self, size=1):
        """
        Generates a lattice of size l_1 x l_2 x ... for size = (l_1, l_2,
        ...). Each lattice site contains a qubit.

        Args:
            size(tuple of int): length of the lattice in each directions
        """
        self.coordinates = [Coordinate(loc) for loc in np.ndindex(size)]
        self.qubits = {v: Qubit(v) for v in self.coordinates}
        self.d = dim_spatial(size)
        self.size = size

    def fine_grain(self, blowup_factor):
        """
        Fine-grain the existing lattice into a larger lattice
        """
        new_size = self.size * blowup_factor
        old_qubits = {
            v * blowup_factor: self.qubits[v] for v in self.coordinates}
        pts = [Coordinate(loc) for loc in np.ndindex(new_size)]
        pts_old = [v * blowup_factor for v in self.coordinates]
        pts_new = [loc for loc in pts if loc not in pts_old]

        for v in pts_old:
            old_qubits[v].label_circuit = v

        for v in pts_new:
            old_qubits[v] = Qubit(v)


def dim_spatial(size):
    if isinstance(size, int):
        return 1
    elif isinstance(size, tuple):
        for n in size:
            if not isinstance(n, int):
                raise TypeError("The size should be specified as integers.")
            elif n < 1:
                return ValueError("Integer is not positive.")
        return len(size)
    else:
        raise TypeError("Input should be an integer or a tuple of integers.")
