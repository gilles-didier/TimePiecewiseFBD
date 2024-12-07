# TimePiecewiseFBD
Maximum likelihood estimation of the skyline FBD model

Software includes

 - 'getML'
	computes the maximum likelihood and the corresponding parameters from a model, a single tree and the fossil stratigraphic intervals

type

	> make getML (or make all)
	in a console opened on the src directory for compiling the software.


Directory "src" contains the C sources of the software

Directory "data" contains the dataset studied in "The rise and fall of Varanopidae† (Amniota, Synapsida), the last ectothermic stem-mammals":


 - 'Simulated_dataset_tree.newick'
  contains the simulated tree From Didier and Laurin (2024)
  
 - 'Simulated_dataset_fossils.csv'
   contains the  corresponding fossil ages

 - 'Varanopidae_trees_17.newick'
  contains 100 equiparsimonious trees of 17 Varanopidae
  
 - 'Varanopidae_trees_21.newick'
  contains 100 equiparsimonious trees of 21 Varanopidae

 - 'Varanopidae_fossils.csv'
   contains the corresponding fossil ages
   
   
 - folders 'Simulated_dataset_models' and 'Varanopidae_models' contain the models considered in the manuscript


A complete description of the options of the software is given below.


------------
 getML 
------------

--------------------------
REQUIREMENT

	'getML' requires the gsl libraries.

--------------------------
COMPILING

	Just type
	> make getML
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'getML' returns the maximum likelihood and the corresponding parameters


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME

	getML - returns the maximum likelihood and the corresponding parameters wrt a dataset under a model specification
	
SYNOPSIS

	getML [OPTIONS] <tree(s)> <fossil ages> <model> [output File]

DESCRIPTION

	Return the maximum likelihood and the corresponding parameters of the dataset made of <tree(s)> and <fossil ages> under the model <model> and write two files, one containing the max likelihood and the corresponding and on containing the corresponding exact fossil ages.

	Options are
	-o <file name> 
		read the NLopt parameters in a file cotaining them, e.g., ':SPE [0;1] :EXT [0;1] :FOS [0:1] :TRI 10 :TOL 1.E-5 :ITE 10000'
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound> <sampling probability upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates and the sampling probability are uniformly drawn before numerical optimisation (lower bounds are all 0)
	-r <number>
		set the number of replicas for the numerical optimisation
	-h
		display help

EXAMPLE

./getML -r 50 working_directory/data/Varanopidae_dataset_trees_21.newick working_directory/data/Varanopidae_dataset_fossils.csv working_directory/data/Varanopidae_dataset_models/model_M0.txt results

