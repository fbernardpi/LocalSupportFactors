# GENERAL
This is an implementation of the method described in [1] for obtaining
linear shape deformation factors that have local support. When using this 
code, please cite [1]. 

Note that this code does not implement the factor splitting procedure. As 
such, there may be multiple active regions in a single factor, see [1] for 
details.

[1] F. Bernard, P. Gemmar, F. Hertel, J. Goncalves, and J. Thunberg. 
"Linear Shape Deformation Models with Local Support Using Graph-based 
Structured Matrix Factorisation". CVPR, 2016.

For feedback, suggestions or other comments please contact Florian Bernard 
per email at f.bernardpi[at]gmail[dot]com .

The FactorTVL1L2 code is provided with kind permission by Ben Haeffele.

# USAGE 
The file demo.m contains a simple toy example that shows how to use the 
method. The main functionality can be accessed through the function
computeLocalSupportDeformationModel().










