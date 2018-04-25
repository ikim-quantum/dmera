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


class gate:
    qubits = []
    unitary = np.eye(0)

    def __init__(self, qubits=None, unitary=None):
        """
        Initialize the gate in terms of the qubits that it acts on
        and the unitary implemented on these qubits.

        Args:
            qubits: The set of qubits on which the gate is acted on.
            unitary: The matrix that specifies the gate.
        """

        if qubits is None:
            # if qubits are not specified, create an empty gate
            self.qubits = []
            self.unitary = []
        elif not isinstance(qubits, list):
            # if qubits are not specified as list, return a TypeError
            raise TypeError("Qubits must be specified as a list.")
        else:
            # if qubits are specified in terms of a list, store it
            self.qubits = qubits
            if unitary is None:
                # if unitary is not specified, initialize it to identity
                self.unitary = np.eye(len(qubits) ** 2)
            elif not isinstance(unitary, np.ndarray):
                # if unitary is notndarray, return a TypeError
                raise TypeError("Unitary must be a numpy ndarray.")
            elif unitary.ndim != 2:
                # if unitary is not a matrix, return a TypeError
                raise TypeError("Unitary must be a matrix.")
            elif unitary.shape[0] != unitary.shape[1]:
                # if unitary is not a square matrix, return a TypeError
                raise TypeError("Unitary must be a square matrix.")
            elif unitary.shape[0] != len(qubits) ** 2:
                # if the unitary dimension is incorrect, return ValueError
                raise ValueError("Dimension must match the qubit list.")
            else:
                # if no error occurs, store the unitary
                self.unitary = unitary

    def has_qubits(self, qubits):
        """
        Checks if the gate acts on qubits.

        Args:
            qubits: Qubits of interest

        Rerturns:
            bool: True if the gate acts on some qubits, False otherwise.
        """
        if not set(qubits).isdisjoint(self.qubits):
            return True
        else:
            return False


class measure:
    qubits = []

    def __init__(self, qubits=None):
        """
        Specify qubits to measure. The measurement is done in the
        computational basis.

        Args:
            qubits: Qubits of interest
        """

        if qubits is None:
            # if qubits are not specified, create an empty list
            self.qubits = []
        elif not isinstance(qubits, list):
            # if qubits are not specified as list, return a TypeError
            print("TypeError : Qubits must be specified as a list.")
            raise TypeError
        else:
            self.qubits = qubits

    def has_qubits(self, qubits):
        """
        Checks if the measurement is applied to the qubits.

        Args:
            qubits: Qubits of interest

        Returns:
            bool: True if at least one qubit is measured, False otherwise.
        """
        if not set(qubits).isdisjoint(self.qubits):
            return True
        else:
            return False


class prepare:
    qubits = []

    def __init__(self, qubits=None):
        """
        Specify qubits to prepare. The preparation is done in the
        computational basis.

        Args:
            qubits: Qubits of interest
        """

        if qubits is None:
            # if qubits are not specified, create an empty list
            self.qubits = []
        elif not isinstance(qubits, list):
            # if qubits are not specified as list, return a TypeError
            raise TypeError("Qubits must be specified as a list.")
        else:
            self.qubits = qubits

    def has_qubits(self, qubits):
        """
        Checks if qubits were prepared.

        Args:
            qubits: Qubits of interest

        Returns:
            True if at least one of the qubits were prepared, False otherwise.
        """
        if not set(qubits).isdisjoint(self.qubits):
            return True
        else:
            return False


def is_a_gate(new_gate):
    """
    Checks if new_gate is a proper gate.

    Args:
        new_gate: Gate of interest

    Returns:
        True if the input is a proper gate, False otherwise.
    """
    if isinstance(new_gate, gate):
        return True
    elif isinstance(new_gate, measure):
        return True
    elif isinstance(new_gate, prepare):
        return True
    else:
        return False
