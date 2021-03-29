Basic Mathworks Matlab implementations for the cardiomyocyte experiments
"How to steer and control ERK and the ERK signalling cascade exemplified by looking at cardiac insufficiency" 
by Tim Breitenbach published under Creative Commons license CC BY-NC-SA 

main_comp_therapies.m	Template main file in which all the relevant 
			settings can be made and the model with its external stimuli are given. Execution
			with Matlab where the other functions have to be in the same folder.
			Here other models can be inserted in the field for the right hand-side f.
			Set to compare different treatment strategies.

main_effective_treatment.m	Template main file for the experiments in which all the 
				relevant settings can be made and the model with its external stimuli 
				are given. Execution with Matlab where the other functions have to be in 
				the same folder. Set to determine the effective treatment

forward.m		Solves for given time curves of external stimuli u the model equations for
			the  state varialbe x to obtain the corresponding time curve

get_xd.m 		Sets the desired steady state

combinatorial_method.m	Is the implementation of Algorithm 1 of the supplementary material; 
			a heuristic method to select external
			stimuli that have a lower target functional value than the unperturbed time 
			evolution of the model

projected_gradient_method.m	Implementatioin of the projected gradient method

SQH_method.m		Implementation of the sequential quadratic Hamiltonian method

createJacobian.m	Creates the Jacobian of the right hand side of the model; Derivatives of f(x,u)
			with respect to x and u

projection.m		Projects each component of u into [0,1]

setControls.m		Sets the value of the external stimuli to 1 for a certain period of time, 
			else zero

backward.m		Solves the adjoint equation for the projected gradient methdod

backward_SQH.m		Solves the adjoint equation for the sequential quadratic Hamilonian method

get_J.m			Calculates the target functional value for the projected gradient method

get_J_SQH.m		Calculates the traget functional value for the sequential quadratic Hamiltonian 
			method

get_gradient.m		Assembles the gradient for the projected gradient method

drawStimuli.m		Draws the time curves of the resulting external stimuli u

drawStates.m		Draws the time curves of the resulting states x

In order to execute the Matlabfile, one needs a Matlab version with a symbolic math toolbox. 
Additionally the parallel computing toolbox is recommended. If this toolbox is not available, 
then just put "for" instead of "parfor" in the function createJacobian.m.
