![Alt Text](https://raw.githubusercontent.com/ikim-quantum/dmera/master/dmera.png)
# What is DMERA?

DMERA is an acronym for deep multi-scale entanglement renormalization ansatz. This ansatz was proposed in [this preprint](https://arxiv.org/abs/1711.07500). DMERA is a tensor network that was specifically designed to be contracted on a quantum computer. One can certainly simulate this contraction on a classical computer, the cost will be huge. The cost scales exponentially with the number of variational parameters on a classical computer, but only linearly on a classical computer.

# Why DMERA?

If you want to do something meaningful in the [noisy intermediate scale quantum era](https://arxiv.org/abs/1801.00862), DMERA will be one compelling tool. Why?

1) DMERA solves a problem of practical interest. Simulation of strongly interacting quantum many-body systems is one of the outstanding problems in physics. On many [theoretical grounds](https://arxiv.org/abs/1711.07500), we expect the ansatz to be an accurate approximation of the ground state of such systems.

2) DMERA is deep. Unlike most near-term quantum algorithms, the depth of DMERA scales with the system size. This gives a larger expressability then circuits that are equipped with smaller depth.

3) DMERA is noise-resilient. Even though DMERA is deep, the circuit is resilient to noise, even without performing error correction! This means that for DMERA calculation, one can circumvent the high overhead cost associated to quantum error correction.

# Immediate Goals
## Higher-dimensional generalizations
The existing code outputs a random instance of DMERA circuit for 1D DMERA. I expect that, with a modest effort, this code be extended to higher-dimensional generalizations of DMERA.

## Data structure for DMERA
It will be desirable to develop a subroutine that, given a DMERA circuit and a local observable, outputs a sequence of gates and measurements that estimates the expectation value of this local observable. (Currently the code simply generates a random DMERA circuit and does this.)

## Variational Contraction
The energy obtained from the existing method may not be variational. There is a contraction scheme that can circumvent this problem, but it hasn't been implemented yet.

## Outputting into other quantum programming languages
Right now only QASM is supported. 