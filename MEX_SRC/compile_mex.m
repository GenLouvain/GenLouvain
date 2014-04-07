% Compile mex
ext=mexext;
mkdir('../private')

mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -std=c++0x" -Imatlab_matrix metanetwork_reduce.cpp matlab_matrix/full.cpp matlab_matrix/sparse.cpp group_index.cpp
movefile(['metanetwork_reduce.',ext],['../private/metanetwork_reduce.',ext]);

mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -std=c++0x" -Imatlab_matrix group_handler.cpp matlab_matrix/full.cpp matlab_matrix/sparse.cpp group_index.cpp
movefile(['group_handler.',ext],['../private/group_handler.',ext]);