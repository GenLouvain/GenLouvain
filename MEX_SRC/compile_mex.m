% Compile mex
% different options for 32bit and 64bit Matlab
ext=mexext;
switch ext
    case {'mexw32','mexglx','mexmac','mexmaci'} %32bit
        arraydims='-compatibleArrayDims';
    case {'mexw64','mexa64','mexmaci64'} %64bit
        arraydims='-largeArrayDims';
    otherwise %potentially new architectures in the future
        warning('unknown mexext %s, assuming 64bit',ext)
        arraydims='-largeArrayDims';
end
mkdir('../private');
setenv('CXXFLAGS',[getenv('CXXFLAGS'),' -std=c++11 -O4']);
if exist('OCTAVE_VERSION','builtin')
    mex -DOCTAVE -Imatlab_matrix metanetwork_reduce.cpp matlab_matrix/full.cpp matlab_matrix/sparse.cpp group_index.cpp
    mex -DOCTAVE -Imatlab_matrix group_handler.cpp matlab_matrix/full.cpp matlab_matrix/sparse.cpp group_index.cpp
    mex -DOCTAVE ../Assignment/assignmentoptimal.c
else
    mex(arraydims,'-Imatlab_matrix','metanetwork_reduce.cpp', 'matlab_matrix/full.cpp', 'matlab_matrix/sparse.cpp', 'group_index.cpp')
    mex(arraydims,'-Imatlab_matrix', 'group_handler.cpp', 'matlab_matrix/full.cpp', 'matlab_matrix/sparse.cpp', 'group_index.cpp')
    mex(arraydims,'../Assignment/assignmentoptimal.c')
end

movefile(['metanetwork_reduce.',ext],['../private/metanetwork_reduce.',ext]);
movefile(['group_handler.',ext],['../private/group_handler.',ext]);
movefile(['assignmentoptimal.',ext],['../Assignment/assignmentoptimal.',ext]);
