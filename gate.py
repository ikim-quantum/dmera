# Author : Isaac Kim
# A class that describes a gate
import numpy as np


class gate:
    qubits = []
    unitary = np.eye(0)

    def __init__(self, qubits=None, unitary=None):
        """
        Initialize the gate in terms of the qubits that it acts on
        and the unitary implemented on these qubits.
        """

        if qubits is None:
            # if qubits are not specified, create an empty gate
            self.qubits = []
            self.unitary = []
        elif not isinstance(qubits, list):
            # if qubits are not specified as list, return a TypeError
            print("TypeError : Qubits must be specified as a list.")
            raise TypeError
        else:
            # if qubits are specified in terms of a list, store it
            self.qubits = qubits
            if unitary is None:
                # if unitary is not specified, initialize it to identity
                self.unitary = np.eye(len(qubits) ** 2)
            elif not isinstance(unitary, np.ndarray):
                # if unitary is notndarray, return a TypeError
                print("TypeError : Unitary must be a numpy ndarray")
                raise TypeError
            elif unitary.ndim != 2:
                # if unitary is not a matrix, return a TypeError
                print("TypeError : Unitary must be a matrix.")
                raise TypeError
            elif unitary.shape[0] != unitary.shape[1]:
                # if unitary is not a square matrix, return a TypeError
                print("TypeError : Unitary must be a square matrix.")
                raise TypeError
            elif unitary.shape[0] != len(qubits) ** 2:
                # if the unitary dimension is incorrect, return ValueError
                print("ValueError : Dimension must match the qubit list.")
                print("This gate uses {} qubits but the dimension of the"
                      " unitary is {}".format(len(qubits), unitary.shape[0]))
                raise ValueError
            else:
                # if no error occurs, store the unitary
                self.unitary = unitary


class measure:
    qubits = []

    def __init__(self, qubits=None):
        """
        Specify qubits to measure. The measurement is done in the
        computational basis.
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


class prepare:
    qubits = []

    def __init__(self, qubits=None):
        """
        Specify qubits to prepare. The preparation is done in the
        computational basis.
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


def is_a_gate(new_gate):
    if isinstance(new_gate, gate):
        return True
    elif isinstance(new_gate, measure):
        return True
    elif isinstance(new_gate, prepare):
        return True
    else:
        return False


def gate_type(new_gate):
    if isinstance(new_gate, gate):
        return "gate"
    elif isinstance(new_gate, measure):
        return "measure"
    elif isinstance(new_gate, prepare):
        return "prepare"
    else:
        return False
