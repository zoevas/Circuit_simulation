# Circuit_simulation
A linear system is solved for Modified Nodal Analysis in order to compute the circuit node voltages using conjugate gradient method and Biconjugate gradient stabilized method.

## Files
- main.c This file runs the MNA method upon a given circuit (netlist)
- mna.c Contains the functions that are used in the mna method
- mna.h The interfaces of the functions used in mna.c
- parser.c The implementation of the parser of the circuit
- parser.h The interface of the functions used in parser.c
- netlist The circuit

## How to run it
Compile: ```
gcc parser.c mna.c main.c -o main -lm```


Run:```
 ./main <input_file_name>```


where input_file_name is the netlist
