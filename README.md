## GenLouvain Version 2.2
### released July 2019

Please cite this code as
    Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla, and Peter J. Mucha,    
    *"A generalized Louvain method for community detection implemented
    in MATLAB,"* https://github.com/GenLouvain/GenLouvain (2011-2019).




## Contents:

This package consists of the main `genlouvain.m` file which calls a number of
subroutines implemented as mex functions. Source code for the mex files is
included in the "MEX_SRC" directory. Pre-compiled executables for 64bit Mac,
Windows, and Linux systems are included in the private directory. It also
includes `iterated_genlouvain.m` which iteratively applies `genlouvain` on the
output partition of the previous run with optional post-processing. Functions
to compute modularity matrices and to post-process partitions are included in
the "HelperFunctions" directory. The post-processing functions solve optimal
assignment problems using code by Markus Buehren (included in the "Assignment"
directory and available at https://uk.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem/content/assignmentoptimal.m).



## Installation instructions:

Make sure that the "GenLouvain" folder and all its subfolders are on the
MATLAB path to ensure that all dependencies between functions are accessible.

If the mex executables for your system are not in the private directory, you
will need to compile these files on your system by running the `compile_mex.m`
script from the "MEX_SRC" directory (check the mex documentation in your MATLAB).
If you would like to share these compiled files with other users, email them to
Peter Mucha (mucha@unc.edu).



## Changes from previous versions:

#### Support for multiple aspects
Version 2.2 of GenLouvain adds support for multilayer networks with multiple
aspects (see "multiaspect.m" in "HelperFunctions").

#### Additional randomization:
Version 2.1 of GenLouvain also a implements a new 'moverandw' option which chooses
moves at random with a probability proportional to the increase in the quality
function. This is in addition to the 'moverand' option from Version 2.0 which chooses
moves uniformly at random from all possible moves that improve the quality function.

#### Increased speed:
Version 2.1 removes quadratic bottlenecks that could become noticeable for very large
networks (millions of nodes). The mex functions have also been optimized further.

#### Generate modularity matrices:
Version 2.1 includes a folder "HelperFunctions" with functions to
generate different types of monolayer and multilayer modularity matrices.

#### Iterated GenLouvain with postprocessing:
Includes `iterated_genlouvain` which iteratively restarts `genlouvain` with the output
partition of the previous run (with optional post-processing). Post-processing functions
for ordered and unordered multilayer partitions that increase the value of the quality
function without changing partitions on each layer are included in "HelperFunctions".
"HelperFunctions" also includes functions that compute "persistence" for ordered and
unordered multilayer networks.



## Usage:

1.  generate a modularity matrix for your network (see `doc('HelperFunctions')`)

2.  use `genlouvain` or `iterated_genlouvain` to obtain a partition that approximately
    optimizes the corresponding modularity-like quality function

3.  ideally repeat step 2 multiple times to check that the output is consistent between
    randomizations

The genlouvain.m function uses different methods for computing the change in
modularity, depending on whether the modularity matrix is provided as a sparse
matrix or not. Depending on the amount of sparsity in the modularity matrix, it may
be faster to convert it to a full matrix.

More extensive documentation and example use of this code is provided online
(http://netwiki.amath.unc.edu/GenLouvain) and in the individual functions (e.g., see
`doc('genlouvain')` and `doc('iterated_genlouvain')`).

***IMPORTANT NOTE:***
When using the multilayer quality function in Mucha et al. 2010, we recommend
using `iterated_genlouvain` with 'moverandw' and the appropriate post-processing
function (i.e., `postprocess_ordinal_multilayer` for an ordered multilayer
network and `postprocess_categorical_multilayer` for an unordered multilayer network)
for better results.

## Acknowledgments:
 A special thank you to Stephen Reid, whose greedy.m code was the
 original version that has over time developed into the present code.

 Thank you also to Dani Bassett, Jesse Blocher, Mason Porter and Simi
 Wang for inspiring improvements to the code.

## References:

Mucha, P. J., Richardson, T., Macon, K., Porter, M. A. & Onnela, J.-P.
Community structure in time-dependent, multiscale, and multiplex networks.
Science 328, 876-878 (2010).

## License:

The codes included in this directory are provided for broad use under
a minor (last line) modification of the "FreeBSD License" (see License.txt)

-------------------------------------------------------------------------------------
***Notes on OCTAVE compatibility:***

The compile_mex.m script from the MEX_SRC directory creates OCTAVE .mex files
when run from OCTAVE.

If you are trying to use this from the old 3.4.0 .app bundle version of OCTAVE for
Mac, you will need to fix OCTAVE's build configuration first (or you may want to
consider upgrading to a recent 3.8.x version where this seems to work out of the
box):

1. Ensure that the environment variables CXX and DL_LD point to a C++ compiler
	installed on your system (e.g. by running
		`setenv(‘CXX’,’/usr/bin/g++’)`
		`setenv(‘DL_LD’,’/usr/bin/g++’)`
	where ‘/usr/bin/g++’ may need to be replaced with the path to your compiler
	depending on your system configuration).

2. Include the ‘-arch i386’ option in CXXFLAGS and LDFLAGS by running
		`setenv('CXXFLAGS',[getenv('CXXFLAGS'),' -arch i386'])`
		`setenv('LDFLAGS',[getenv('LDFLAGS'),' -arch i386'])`
	to create 32bit binaries.

3. Change line 52 of
	/Applications/Octave.app/Contents/Resources/include/octave-3.4.0/octave/mexproto.h
	from `#include <cstdlib>` to `#include <stdlib.h>` to
	avoid a conflict from including two different versions of the standard
	library.

4. Finally run `compile_mex` to compile the binaries.

-------------------------------------------------------------------------------------
