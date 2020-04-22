# TreeRing
A framework for exact moment propagation for nonlinear systems.

## Installation and Dependencies
To install the package, you can simply run `pip3 install -e .` in the root of this repo.

## Example
`example/example_uncontrolled_agent.py` provides an example of using TreeRing on the Dubin's system we use both with the reduced and un-reduced moment update form (MUF). As it is configured right now, running it will run TreeRing to search for a complete moment basis w.r.t. the reduced MUF and print out the moments in the moment basis along with the dynamics of the moments.

## Variable Name Syntax
`StateVariable` and `DisturbanceVariable` should be letters, no special symbols.

When we create instances of `BasisVariable` as the product of multiple state variables, the composing variables have their names concatenated:
`statevarpow1_statevarpow2`