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

	getAIC - eturns the maximum likelihood and the corresponding parameters wrt a dataset under a model specification
	
SYNOPSIS

	getAIC [OPTIONS] <tree(s)> <fossil ages> <model specification> [output File]

DESCRIPTION

	Return the maximum likelihood and the corresponding parameters of the dataset made of <tree(s)> and <fossil ages> under the model specification <model specification>.

	Options are
	-w <speciation rate width> <extinction rate width> <fossilization rate width> <sampling probability width>
		set the widths of the sliding windows used for sampling the speciation, extinction and fossilization rates and of the sampling probability during the MCMC
	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound> <sampling probability upper bound>
		set the upper bounds of the interval in which the speciation, extinction and fossilization rates and the sampling probability are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)
	-f <proportion>
		set the proportion of moves in the parameters in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)
	-a <speciation proportion> <extinction proportion> <fossilization proportion>
		set the relation proportion of moves in the speciation, the extinction and the fossilzation rates (thus in the sampling probability)
	-s <number>
		set the number of samples required to estimate the maximum likelihood
	-r <number>
		set the random seed
	-h
		display help

EXAMPLE

./getAIC -s 10000 -f 0.2  -a 0.25 0.25 0.25 -w 0.05 0.05 0.05 0.25 -i 0.5 0.5 0.5 1.  ../data/Simulated_dataset_tree.newick ../data/Simulated_dataset_fossils.csv ../data/Simulated_Dataset_model_spec/Simul_model_spec_0.txt 

