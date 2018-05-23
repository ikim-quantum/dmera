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

from circuit import Circuit
from lattice import Lattice
import gate as gt

class DMERA(Circuit):
    def __init__(self, depth, scale, *unit):
        super().__init__()
        self.D = depth
        self.n = scale
        self.lattice = Lattice(unit)

        for s in self.n:
            self.expand(2)
            for d in self.D:
                self.entangle_nn()

    @property
    def idle_qubits(self):
        """
        Returns a list of qubits which are in the lattice but not in the circuit.
        """
        qubits_lattice = self.lattice.qubits.values()
        qubits_circuit = self.qubits
        return [qubit for qubit in qubits_lattice if qubit not in qubits_circuit]

    def initialize_idle_qubits(self):
        """
        Initialize idle qubits.
        """
        self.extend([gt.Prepare(qubit) for qubit in self.idle_qubits])

    def expand(self, factor):
        """
        Expand the lattice and intialize the newly introduced qubits.
        """
        self.lattice.expand(factor)
        self.initialize_idle_qubits()

    def entangle_nn(self, unitary=None):
        """
        Entangle nearest neighbors.
        """
        dim = self.lattice.d
        # Create unit vectors in each directions.

        # Iterate over:
        # Create sublattice in each directions.

        # Make a list of qubits on the sublattice and the shifted sublattice.

        # Add entangling gates.
        self.extend(Circuit.doubles(family1, family2, unitary))
