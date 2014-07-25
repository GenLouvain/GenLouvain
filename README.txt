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

More extensive example use of this code is provided online
(http://netwiki.amath.unc.edu/GenLouvain).

The codes included in this directory are provided for broad use under
a minor (last line) modification of the "FreeBSD License":

Copyright (c) 2012-2014, Inderjit S. Jutla, Lucas G. S. Jeub, and Peter J. Mucha
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of entities associated with the authors.
