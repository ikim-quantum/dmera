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

import gate as gt


class circuit(list):
    """
    The circuit class is inherited from Python's original list class.
    Circuit class overrides the list methods in a way that forces the
    list element to be a gate object.
    """

    def __init__(self, input_circuit=[]):
        """
        Initiazlies the circuit. Default gate input is an empty list.

        Args:
            input_circuit (list): Initial specification of the circuit.
        """
        self.extend(input_circuit)

    def __pow__(self, exponent):
        """
        Returns the circuit applied multiple times.

        Args:
            exponent (int): Repeat the circuit this many times.

        Returns:
            circuit: Returns the circuit applied exponent times.
        """
        out_circuit = circuit()
        for n in range(exponent):
            out_circuit.extend(self)
        return out_circuit

    def __ipow__(self, exponent):
        """
        Update the circuit to the same circuit applied multiple times.

        Args:
            exponent (int): Repeat the circuit this many times.

        Returns:
            circuit: Returns the circuit exponent times.
        """
        out_circuit = circuit() 
        for n in range(exponent):
            out_circuit.extend(self)
        self = out_circuit
        return out_circuit

    def append(self, new_gate):
        """
        Appends new_gate. For measure/prepare, qubit list is broken down
        to operations on individual qubits.

        Args:
            new_gate (gate): Append new_gate to the circuit.
        """
        if gt.is_a_gate(new_gate):
            if isinstance(new_gate, gt.gate):
                super().append(new_gate)
            elif isinstance(new_gate, gt.prepare):
                for qubit in new_gate.qubits:
                    super().append(gt.prepare([qubit]))
            elif isinstance(new_gate, gt.measure):
                for qubit in new_gate.qubits:
                    super().append(gt.measure([qubit]))
        else:
            raise TypeError("Input is not a proper gate.")

    def extend(self, other_circuit):
        """
        Extend the existing circuit by other_circuit. For measure/prepare,
        qubit list is broken down to operations on individual qubits.

        Args:
            other_circuit (circuit): Extend other_circuit to the circuit.
        """
        for gate in other_circuit:
            if not gt.is_a_gate(gate):
                raise TypeError("Input must be a circuit.")
        else:
            for gate in other_circuit:
                self.append(gate)

    def index(self, gate_or_qubit):
        """
        Returns the index of the gate_or_qubit

        Args:
            gate_or_qubit: It can be either gate datatype or anything else.

        Returns:
            int: If gate_or_qubit is (gate), return the first index.
                 Else, return the index of the first gate acting gate_or_qubit.
        """
        if gt.is_a_gate(gate_or_qubit):
            return super().index(gate_or_qubit)
        else:
            n = 0
            for gate in self:
                n += 1
                if gate_or_qubit in gate.qubits:
                    return n - 1
            raise ValueError(
                "Qubit {} not in the circuit".format(gate_or_qubit))

    def is_parallel(self):
        """
        Checks if the circuit can be run in parallel.

        Returns:
            bool: True if circuit can be run in parallel, False otherwise.
        """
        qubits = []
        for gate in self:
            for qubit in gate.qubits:
                if qubit in qubits:
                    return False
                else:
                    qubits.append(qubit)
        else:
            return True

    def prepared(self, qubit):
        """
        Checks if the qubit was prepared properly. If the qubit is present
        in the circuit and if it is prepared before other gates are enacted
        on this qubit, the qubit is properly prepared. Otherwise, the qubit
        is not properly prepared.

        Args:
            qubit: This method checks if the qubit is prepared properly.

        Returns:
            bool: True if the qubit is prepared properly, False otherwise.
        """
        try:
            if isinstance(self[self.index(qubit)], gt.prepare):
                return True
            else:
                return False
        except ValueError:
            return False

    def all_prepared(self):
        """
        Check if every qubit in the circuit is prepared properly.

        Returns:
            bool: True if every qubit is prepared properly, False otherwise.
        """
        for qubit in self.qubits():
            if not self.prepared(qubit):
                return False
        return True

    def find_gates(self, qubits):
        """
        Returns all gates that act on the qubits.

        Args:
            qubits: Qubits of interest

        Returns:
            list: List of gates that act on the qubits.
        """
        my_gates = []
        for gate in self:
            for qubit in qubits:
                if qubit in gate.qubits:
                    my_gates.append(gate)
        return my_gates

    def qubits(self):
        """
        Returns a list of qubits that appeared in the circuit.

        Returns:
            list: List of qubits in the circuit.
        """
        qubits = []
        for gate in self:
            for qubit in gate.qubits:
                if qubit not in qubits:
                    qubits.append(qubit)
        return qubits

    def has_qubits(self, qubits):
        """
        Checks if all the qubits are present in the circuit.

        Args:
            qubits: Qubits of interest

        Returns:
            bool: True if all the qubits appear, False otherwise.
        """
        circuit_qubits = self.qubits()
        for qubit in qubits:
            if qubit not in circuit_qubits:
                return False
        return True

    def pcc(self, qubits):
        """
        Returns the past causal cone.
        * Past causal cone refers to an ordered list of gates that can
        influence the reduced density matrix of the observables of the
        qubits. Every circuit element outside the past causal cone has
        absolutely no influence on such observables.

        Args:
            qubits: Qubits of interest

        Returns:
            circuit: Past causal cone of qubits.
        """
        tracked_qubits = qubits
        my_pcc = circuit()
        for gate in reversed(self):
            if gate.has_qubits(qubits):
                tracked_qubits += gate.qubits
                my_pcc.append(gate)
        return my_pcc[::-1]
