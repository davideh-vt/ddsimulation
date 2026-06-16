# ddsimulation
Contains code used for the simulations for the linearized hyperbolic model of delay discounting paper.

The file non-lin-ln-k-iterative-3-pt-approx.R contains a function to estimate ln(k) values for individual subjects under the nonlinear hyperbolic model.

The file t1e-power-simulations.R contains a script to simulate a dataset under our linearized hyperbolic model, analyze it using the linearized model, and nonlinear hierarchical and nonlinear individual procedures.
It saves information useful for analyzing the mean square errors for the estimated subject and group mean parameters using each method, and the power of the hypothesis test for group mean equality.

The file individual-ln-k-sim.R contains a script to simulate a dataset under our linearized hyperbolic model and save the individual subject ln(k) values, along with the estimated subject ln(k) by each estimator considered in the paper.

The file plot-code.R gets useful summary information from the simulation data, and plots it. It specifically contains the code for reproducing plots 3, 4, and 5 from the paper.
