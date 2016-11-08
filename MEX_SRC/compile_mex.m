% Compile mex
ext=mexext;
mkdir('../private');
setenv('CXXFLAGS',[getenv('CXXFLAGS'),' -std=c++11']);
if exist('OCTAVE_VERSION','builtin')
    
    
    mex -DOCTAVE -Imatlab_matrix metanetwork_reduce.cpp matlab_matrix/full.cpp matlab_matrix/sparse.cpp group_index.cpp

    mex -DOCTAVE -Imatlab_matrix group_handler.cpp matlab_matrix/full.cpp matlab_matrix/sparse.cpp group_index.cpp
else
    mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -O4" -Imatlab_matrix metanetwork_reduce.cpp matlab_matrix/full.cpp matlab_matrix/sparse.cpp group_index.cpp

    mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -O4" -Imatlab_matrix group_handler.cpp matlab_matrix/full.cpp matlab_matrix/sparse.cpp group_index.cpp
end

movefile(['metanetwork_reduce.',ext],['../private/metanetwork_reduce.',ext]);
movefile(['group_handler.',ext],['../private/group_handler.',ext]);
