########################### 
# Tests for Lattice class #
# 9/15/2020               #
# Isaac H. Kim            #
###########################
import numpy as np
from lattice.lattice import Lattice

n_trials = 100

# Shift tests
for i in range(n_trials):
    d = np.random.randint(3)+1
    size = tuple(np.random.randint(10, size=d) + 1)
    l = Lattice(*size)
    shift = np.random.randint(10, size=d)

    k = np.random.randint(len(l.pts))
    a = l.qubits[tuple(l.pts[k])]
    l += shift
    l -= shift
    b = l.qubits[tuple(l.pts[k])]

    assert a==b
    

