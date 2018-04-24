# Authors : Isaac Kim
# Circuit is a sequence of gates
import gate as gt


class circuit(list):
    def __init__(self, gates=[]):
        for gate in gates:
            self.append(gate)

    def append(self, new_gate):
        # append a gate
        if gt.is_a_gate(new_gate):
            super().append(new_gate)
        else:
            print("TypeError : input is not a proper gate.")
            raise TypeError

    def is_parallel(self):
        # returns True if and only if parallelizable
        qubits = []
        for gate in super():
            for qubit in gate.qubits:
                if qubit in qubits:
                    return False
                else:
                    qubits += qubit
        else:
            return True

    def index(self, gate_or_qubit):
        if gt.is_a_gate(gate_or_qubit):
            return super().index(gate_or_qubit)
        # return the index of the first gate that contains qubit
        else:
            n = 0
            for gate in super():
                n += 1
                if gate_or_qubit in gate.qubits:
                    return n-1
                print("Qubit {} was not used.".format(gate_or_qubit))
            return False

    def prepared(self, qubit):
        # return true if qubit was prepared before the other gates
        if gt.gate_type(super()[super().index(qubit)]) == "prepare":
            return True
        else:
            return False

    def all_prepared(self):
        # return true if every qubit was prepared before the other gates
        for gate in super():
            if not super().prepared(gate):
                return False
        return True

    def find_gates(self, qubits):
        # return a list of gates that act on the qubits
        my_gates = circuit()
        for qubit in qubits:
            for gate in super():
                if qubit in gate.qubits:
                    my_gates.append(gate)
        return my_gates
