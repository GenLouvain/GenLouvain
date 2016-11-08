README - GENLOUVAIN VERSION 2.0
released July 2014

Please cite this code as
    Inderjit S. Jutla, Lucas G. S. Jeub, and Peter J. Mucha,    
    "A generalized Louvain method for community detection implemented
    in MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2011-2014).

This package consists of the main genlouvain.m file which calls a number of
subroutines implemented as mex functions.

Source code for the mex files is included in the MEX_SRC directory.
Pre-compiled executables for 64bit Mac, Windows, and Linux systems are included
in the private directory.

If the mex executables for your system are not in the private
directory, you will need to compile these files on your system by running the 
compile_mex.m script from the MEX_SRC directory (check the mex documentation in 
your MATLAB). If you would like to share these compiled files with other users, 
email them to Peter Mucha (mucha@unc.edu).

The genlouvain.m function uses different methods for computing the change in
modularity, depending on whether the modularity matrix is provided as a sparse 
matrix or not. Version 2.0 of genlouvain.m should be much faster than 
version 1.2 of genlouvain.m (no mex) for full matrices and faster than 
genlouvainmex.m for sparse matrices.

Version 2.0 of genlouvain also implements a new 'randmove' option which choses 
moves uniformly at random from all possible moves that improve modularity 
(instead of choosing the move that maximally improves modularity).

More extensive example use of this code is provided online
(http://netwiki.amath.unc.edu/GenLouvain).

Acknowledgments:
 A special thank you to Stephen Reid, whose greedy.m code was the
 original version that has over time developed into the present code, 
 and Marya Bazzi for noticing the problematic behavior of genlouvain for
 ordinal interslice coupling and contributing code that developed into the 
 'randmove' option.
 Thank you also to Dani Bassett, Jesse Blocher, Mason Porter and Simi
 Wang for inspiring improvements to the code.

——————————————————————————————————————————————————————————————————————————————————
Notes on OCTAVE compatibility:

The compile_mex.m script from the MEX_SRC directory creates OCTAVE .mex files 
when run from OCTAVE. 

If you are trying to use this from the old 3.4.0 .app bundle version of OCTAVE for 
Mac, you will need to fix OCTAVE's build configuration first (or you may want to 
consider upgrading to a recent 3.8.x version where this seems to work out of the 
box):

1. Ensure that the environment variables CXX and DL_LD point to a C++ compiler 
	installed on your system (e.g. by running 
		setenv(‘CXX’,’/usr/bin/g++’)
		setenv(‘DL_LD’,’/usr/bin/g++’)
	where ‘/usr/bin/g++’ may need to be replaced with the path to your compiler
	depending on your system configuration).

2. Include the ‘-arch i386’ option in CXXFLAGS and LDFLAGS by running
		setenv('CXXFLAGS',[getenv('CXXFLAGS'),' -arch i386'])
		setenv('LDFLAGS',[getenv('LDFLAGS'),' -arch i386'])
	to create 32bit binaries.

3. Change line 52 of 
	/Applications/Octave.app/Contents/Resources/include/octave-3.4.0/octave/mexproto.h
	from “#include <cstdlib>” to “#include <stdlib.h>” (without quotes) to 
	avoid a conflict from including two different versions of the standard 
	library.

4. Finally run compile_mex to compile the binaries.
——————————————————————————————————————————————————————————————————————————————————


The codes included in this directory are provided for broad use under
a minor (last line) modification of the "FreeBSD License" (see License.txt)
