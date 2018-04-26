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


class Qubit:
    """
    Qubit class. Even though this class contains no method, this class
    will become handy in recursively generating DMERA ansatz at different
    scales.
    """
    def __init__(self, label_circuit="", label_physical=""):
        """
        Initialize the qubit.

        Args:
            label_circuit: The label of the qubit in the circuit diagram.
            label_physical: Assignment of this qubit in the hardware.
        """
        self.label_circuit = label_circuit
        self.label_physical = label_physical

    def __str__(self):
        """
        Returns:
            str: Returns a string with circuit and physical qubit.
        """
        return "[Circuit={}, Physical={}]".format(self.label_circuit,
                                                  self.label_physical)
