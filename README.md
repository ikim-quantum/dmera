# What is DMERA?

DMERA is an acronym for deep multi-scale entanglement renormalization ansatz. 

# Immediate Goals
## Higher-dimensional generalizations
The existing code outputs a random instance of DMERA circuit for 1D DMERA. I expect that, with a modest effort, this code be extended to higher-dimensional generalizations of DMERA.

## Data structure for DMERA
It will be desirable to develop a subroutine that, given a DMERA circuit and a local observable, outputs a sequence of gates and measurements that estimates the expectation value of this local observable. (Currently the code simply generates a random DMERA circuit and does this.)

## Variational Contraction
The energy obtained from the existing method may not be variational. There is a contraction scheme that can circumvent this problem, but it hasn't been implemented yet.

## Outputting into other quantum programming languages
Right now only QASM is supported. 