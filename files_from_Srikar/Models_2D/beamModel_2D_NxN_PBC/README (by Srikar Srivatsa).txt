Beam_2D_NxN_PBC
**refer to MATLAB file for full code breakdown

PATHS TO CHANGE FOR EACH USER:
-in the parent MATLAB file: 
	- the path for the ANSYS executable (used in the system call) (line 70)
	- the ANSYS license type (currently aa_t_i in default code) (line 82)
-in the APDL text file (Truss2D_NxN): 
	- (8x) the path for the input text file, in eight separate read statements in the input section
	- (1x) the path for the output text file, in the write statement at the end