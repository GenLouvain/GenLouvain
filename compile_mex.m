% Compile mex

mex -largeArrayDims metanetwork_reduce.cpp full.cpp sparse.cpp group_index.cpp

mex -largeArrayDims group_handler.cpp full.cpp sparse.cpp group_index.cpp
