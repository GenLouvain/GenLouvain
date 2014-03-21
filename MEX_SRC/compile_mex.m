% Compile mex
ext=mexext;
mkdir('../private')

mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -std=c++11" metanetwork_reduce.cpp full.cpp sparse.cpp group_index.cpp
movefile(['metanetwork_reduce.',ext],['../private/metanetwork_reduce.',ext]);

mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -std=c++11" group_handler.cpp full.cpp sparse.cpp group_index.cpp
movefile(['group_handler.',ext],['../private/group_handler.',ext]);