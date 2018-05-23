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
import matplotlib.pyplot as plt
from qubit import Qubit


class Point(tuple):
    def __init__(self, *v):
        """
        Points on a square lattice. Inherited from the tuple class.
        """
        self = tuple(v)

    def __mul__(self, other):
        """
        Elementwise multiplication.

        Args:
            self (Point): self.
            other (tuple): Other Point/tuple instance.

        Returns:
            (Point): Elementwise multiplication of self and other.

        Example:
            Suppose we have a Point instance called pt=Point((2,3,4)).

            >>> pt * (1, 2, 3)
            (2, 6, 12)
        """
        if isinstance(other, int):
            return Point([n * other for n in self])
        elif isinstance(other, tuple):
            return Point([n * factor for n, factor in zip(self, other)])

    __rmul__ = __mul__

    def __imul__(self, other):
        """
        Elementwise multiplication.

        Args:
            self (Point): self.
            other (tuple): Other Point/tuple instance.

        Example:
            Suppose we have a Point instance called pt=Point((2,3,4)).

            >>> pt *= (1, 2, 3)
            >>> print(pt)
            (2, 6, 12)
        """
        self = self * other
        return self

    def __mod__(self, other):
        """
        Elementwise %.

        Args:
            self (Point): self.
            other (tuple): Other Point/tuple instance.

        Example:
            Suppose we have a Point instance called pt=Point((2,3,4)).

            >>> pt % (1, 2, 3)
            (0, 1, 1)
        """
        return Point([i % j for i, j in zip(self, other)])

    def __imod__(self, other):
        """
        Elementwise %.

        Args:
            self (Point): self.
            other (tuple): Other Point/tuple instance.

        Example:
            Suppose we have a Point instance called pt=Point((2,3,4)).

            >>> pt %= (1, 2, 3)
            >>> print(pt)
            (0, 1, 1)
        """
        self = self % other
        return self

    def __add__(self, other):
        """
        Elementwise addition.

        Args:
            self (Point): self.
            other (tuple): Other Point/tuple instance.

        Example:
            Suppose we have a Point instance called pt=Point((2,3,4)).

            >>> pt + (1, 2, 3)
            (3, 5, 7)
        """
        return Point([i + j for i, j in zip(self, other)])

    def __iadd__(self, other):
        """
        Elementwise addition.

        Args:
            self (Point): self.
            other (tuple): Other Point/tuple instance.

        Example:
            Suppose we have a Point instance called pt=Point((2,3,4)).

            >>> pt += (1, 2, 3)
            >>> print(pt)
            (3, 5, 7)
        """
        self = self + other
        return self

    def __sub__(self, other):
        """
        Elementwise subtraction.

        Args:
            self (Point): self.
            other (tuple): Other Point/tuple instance.

        Example:
            Suppose we have a Point instance called pt=Point((2,3,4)).

            >>> pt - (1, 2, 3)
            (1, 1, 1)
        """
        return Point([i - j for i, j in zip(self, other)])

    def __isub__(self, other):
        """
        Elementwise subtraction.

        Args:
            self (Point): self.
            other (tuple): Other Point/tuple instance.

        Example:
            Suppose we have a Point instance called pt=Point((2,3,4)).

            >>> pt -= (1, 2, 3)
            >>> print(pt)
            (1, 1, 1)
        """
        self = self - other
        return self


class Lattice():
    def __init__(self, *args):
        """
        Generates a lattice of size l_1 x l_2 x ... for size = (l_1, l_2,
        ...). Each lattice site contains a qubit.

        Args:
            size(list of int): length of the lattice in each directions

        Examples:
            The following code will create an instance my_lattice, whose size
            is 2 x 3.

            >>> my_lattice = Lattice(2, 3)
            >>> print(my_lattice.pts)
            [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2)]
        """
        if args:
            self.pts = [Point(loc) for loc in np.ndindex(args)]
            self.qubits = {v: Qubit(v) for v in self.pts}
        else:
            self.pts = []
            self.qubits = {}

    def __iter__(self):
        """
        Iterate over the points.

        Example:
            One can iterate over the points in the lattice in a following way.

            >>> my_lattice = Lattice(2,2)
            >>> for pt in my_lattice:
                    print(pt)
            (0, 0)
            (0, 1)
            (1, 0)
            (1, 1)
        """
        return np.ndindex(self.size)

    def __add__(self, shift):
        """
        Shift the qubits by shift. Periodic boundary condition is assumed.

        Args:
            shift(list of int): Lattice vector for the shift.

        Returns:
            Lattice : New lattice instance. The qubit at pt is the qubit at pt
                      + shift in the original instance.

        Examples:
            Let's say we create a Lattice instance with dimension 3 x 2, and a
            new lattice that is shifted by a vector (1, 1).

            >>> my_lattice = Lattice(3, 2)
            >>> new_lattice = my_lattice + (1, 1)


            These two lattices contain the same lattice points. However, the
            assignment of the qubit is different.

            >>> my_lattice.qubits[0,0] == new_lattice.qubits[0,0]
            False


            This is because the new_lattice is shifted by (1,1). Once we
            account for this shift, one can see that the qubits can be
            identified with each other.

            >>> my_lattice.qubits[1,1] == new_lattice.qubits[0,0]
            True
        """
        out = Lattice()
        out.qubits = {(v + Point(shift)) %
                      self.size: self.qubits[v] for v in self.pts}
        out.pts = [pt for pt in self.pts]
        return out

    def __iadd__(self, shift):
        """
        Shift the qubits by shift. Periodic boundary condition is assumed. See
        __add__ for the details/examples.

        Args:
            shift(list of int): Lattice vector for the shift.
        """
        self.qubits = {(v + Point(shift)) %
                       self.size: self.qubits[v] for v in self.pts}
        return self

    def __sub__(self, shift):
        """
        Shift the qubits by -shift. Periodic boundary condition is assumed.

        Args:
            shift(list of int): Lattice vector for the shift.

        Returns:
            Lattice : New lattice instance. The qubit at pt is the qubit at pt
                      + shift in the original instance.

        Examples:
            Let's say we create a Lattice instance with dimension 3 x 2, and a
            new lattice that is shifted by a vector (-1, -1).

            >>> my_lattice = Lattice(3, 2)
            >>> new_lattice = my_lattice - (1, 1)


            These two lattices contain the same lattice points. However, the
            assignment of the qubit is different.

            >>> my_lattice.qubits[0,0] == new_lattice.qubits[0,0]
            False


            This is because the new_lattice is shifted by (1,1). Once we
            account for this shift, one can see that the qubits can be
            identified with each other.

            >>> my_lattice.qubits[0,0] == new_lattice.qubits[1,1]
            True
        """
        out = Lattice()
        out.qubits = {(v - Point(shift)) %
                      self.size: self.qubits[v] for v in self.pts}
        out.pts = [pt for pt in self.pts]
        return out

    def __isub__(self, shift):
        """
        Shift the qubits by -shift. Periodic boundary condition is assumed. See
        __sub__ for the details/examples.

        Args:
            shift(list of int): Lattice vector for the shift.
        """
        self.qubits = {(v - Point(shift)) %
                       self.size: self.qubits[v] for v in self.pts}
        return self

    @property
    def size(self):
        """
        Returns the size of the lattice.

        Returns:
            tuple: length of the lattice in each directions.

        Example:
            >>> my_lattice = Lattice(2,3,4)
            >>> my_lattice.size
            (2, 3, 4)
        """
        return Point([max([pt[i] for pt in self.pts]) + 1 for
                      i in range(self.d)])

    @property
    def d(self):
        """
        Returns the dimension of the lattice.

        Returns:
            int: Dimension of the lattice.

        Examples:
            >>> my_lattice = Lattice(3,6,5)
            >>> my_lattice.d
            3
            >>> my_lattice = Lattice(4,4)
            >>> my_lattice.d
            2
        """
        return dim_spatial(self.pts[0])

    def draw(self, *v_sublattice):
        """
        Draws lattice. Only works for 2D lattice.

        Args:
            v_sublattice(tuple of int): The lattice vector that defines a
                                        sublattice. The first component is
                                        x component, and the second component
                                        is the y component.

        Examples:
            The following command draws a 4 x 4 lattice. Sites with even
            x coordinate are colored as blue. Sites with odd x coordinate are
            colored as red.

            >>> my_lattice = Lattice(4,4)
            >>> my_lattice.draw(2,1)

            The coloring pattern changes depending on the input. In the
            following example, sites with even y coordinate are colored as
            blue. Sites with odd y coordinate are colored as red.
        """
        try:
            if self.d == 1:
                for site in self.restrict(*v_sublattice):
                    nbhr = site + v_sublattice - [1]
                    plt.plot([site, nbhr], [0, 0], color='k')
                    plt.plot(site, 0, "o", color='blue')
                    plt.plot(nbhr, 0, "o", color='red')

            else:
                for site in self.restrict(*v_sublattice):
                    nbhr = (site + v_sublattice) + (-1, -1)
                    plt.plot([site[0], nbhr[0]], [site[1], nbhr[1]], color='k')
                    plt.plot(site[0], site[1], "o", color='blue')
                    plt.plot(nbhr[0], nbhr[1], "o", color='red')

        except ValueError:
            print("Only 2D is supported.")

        plt.axis('off')
        plt.show()

    def restrict(self, *v):
        """
        Restrict to lattice sites which are integer multiples of v.
        Returns a Lattice instance which only contains these lattice points.
        The qubits of the new instance correspond to the qubits of the original
        instance.

        Args:
            v (tuple): A lattice vector that defines a sublattice.

        Returns:
            Lattice: Returns a restricted lattice, whose points are integer
                     multiples of v. The qubits of this instance is identified
                     with the qubits of the original instance.
        """
        slattice = Lattice()
        slattice.pts = [Point(pt) for pt in self.pts if is_zero(pt % v)]
        slattice.qubits = {v: self.qubits[v] for v in slattice.pts}
        return slattice

    def expand(self, *blowup_factor):
        """
        Expand the existing lattice to a larger lattice.
        Returns a Lattice instance which are expanded by the blowup_factor.
        The qubits of the original instance is used as some of the qubits
        of the new instance. New qubits are created as well.

        Args:
            blowup_factor (tuple): A factor by which each direction is blown
            up.

        Returns:
            Lattice: Returns an expanded lattice. The points in the original
                     instance is blown up by the blowup_factor. Missing sites
                     between these sites is filled up. The qubits on the blown
                     up sites are identified with the qubits of the original
                     instance. For the new sites, new qubits are created.
        """
        if isinstance(blowup_factor, int):
            factor = tuple([blowup_factor for i in range(self.d)])
        else:
            factor = tuple(blowup_factor)
        # set of new points
        pts = [Point(loc) for loc in np.ndindex(self.size * factor)]
        pts_inherited = [v * factor for v in self.pts]
        pts_new = [loc for loc in pts if loc not in pts_inherited]

        # embed the old qubits to a new lattice
        elattice = Lattice()
        elattice.pts = pts
        elattice.qubits = {v * factor: self.qubits[v] for v in self.pts}
        # change the label_circuit for the old qubits
        for v in pts_inherited:
            elattice.qubits[v].label_circuit = v
        # introduce the new qubits
        for v in pts_new:
            elattice.qubits[v] = Qubit(v)
        return elattice

    def nnpairs(self):
        """
        Returns pairs of nearest neighbors.

        Returns:
            list: The element of the list is another list,
                  which is a pair of lattice instances. The first
                  lattice in the pair is a sublattice, and the 
                  second lattice is a sublattice shifted by a unit
                  vector.
        """
        pairs = []
        for direction in range(self.d):
            vec = unit_vector(direction, self.d)
            ones = one_vector(self.d)
            mymod = vec + ones
            first = self.restrict(*mymod)
            shifted = self - vec
            second = shifted.restrict(*mymod)
            shift2 = self + vec
            third = shift2.restrict(*mymod)
            pairs.append([first, second])
            pairs.append([first, third])
        return pairs

def one_vector(dim):
    """
    Returns a Point instance that has 1s.

    Args:
        dim (int): Number of spatial dimensions.

    Returns:
        Point: Unit vector of length (dim), filled with 1s.
    """
    return Point([1] * dim)
    
def unit_vector(direction, dim):
    """
    Returns a Point instance that represents a unit vector.

    Args:
        direction (int): The direction of the unit vector.
        dim (int): Number of spatial dimensions.

    Returns:
        Point: Unit vector in the direction.

    Example:
        >>> unit_vector(3, 5)
        (0, 0, 0, 1, 0)
    """
    out = [0] * dim
    out[direction] = 1
    return Point(out)

def dim_spatial(size):
    """
    Returns a number of spatial dimensions.

    Args:
        size (tuple): Integer-valued tuple which specifies the size.

    Returns:
        int: Number of spatial dimensions.
    """
    if isinstance(size, int):
        return 1
    elif isinstance(size, tuple):
        for n in size:
            if not isinstance(n, int):
                raise TypeError("The size should be specified as integers.")
            elif n < 0:
                return ValueError("Integer is negative.")
        return len(size)
    else:
        raise TypeError("Input should be an integer or a tuple of integers.")


def is_zero(coordinate):
    """
    Checks if the coordinate(specified in terms of a tuple) is at the origin.

    Args:
        coordinate (tuple): Input coordinate.

    Returns:
        bool: True if coordinate is at the origin, False otherwise.
    """
    for n in range(len(coordinate)):
        if coordinate[n] != 0:
            return False
    return True
