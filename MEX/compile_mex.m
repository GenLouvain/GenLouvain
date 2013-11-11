% Compile mex

mex -largeArrayDims metanetwork_reduce.cpp full.cpp sparse.cpp group_index.cpp

mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -std=c++0x" group_handler.cpp full.cpp sparse.cpp group_index.cpp
