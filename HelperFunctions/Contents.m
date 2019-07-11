% HELPERFUNCTIONS
%
% Monolayer modularity:
%
%   modularity                         - returns monolayer Newman-Girvan modularity matrix for network given by adjacency matrix A, matrix version
%   modularity_f                       - returns monolayer Newman-Girvan modularity matrix for undirected network given by adjacency matrix A, function handle version
%   modularitydir_f                    - returns monolayer Leicht-Newman modularity matrix for directed network given by adjacency matrix A, function handle version
%   bipartite                          - returns monolayer Barber modularity matrix for undirected bipartite networks, matrix version
%   bipartite_f                        - returns monolayer Barber modularity matrix for undirected bipartite networks, function handle version
%
%
% Multilayer modularity (categorical coupling):
%
%   multicat                           - returns multilayer Newman-Girvan modularity matrix for unordered undirected layers, matrix version
%   multicat_f                         - returns multilayer Newman-Girvan modularity matrix for unordered undirected layers, function handle version
%   multicatbipartite                  - returns multilayer Barber modularity matrix for unordered undirected bipartite networks, matrix version
%   multicatbipartite_f                - returns multilayer Barber modularity matrix for unordered undirected bipartite networks, function handle version
%
%
% Multilayer modularity (ordinal coupling):
%
%   multiord                           - returns multilayer Newman-Girvan modularity matrix for ordered layers, matrix version
%   multiord_f                         - returns multilayer Newman-Girvan modularity matrix for ordered undirected layers, function handle version
%   multiorddir_f                      - returns multilayer Leicht-Newman modularity matrix for ordered directed layers, function handle version
%   multiordbipartite                  - returns multilayer Barber modularity matrix for ordered undirected bipartite networks, matrix version
%   multiordbipartite_f                - returns multilayer Barber modularity matrix for ordered undirected bipartite networks, function handle version
%
%
% Multilayer modularity with multiple aspects (supports mix of ordinal and categorical coupling)
%
%   multicat                           - returns multilayer Newman-Girvan modularity matrix for a multiaspect multilayer network with a mix of ordered and categorical coupling
%
% Postprocessing functions:
%
%   postprocess_categorical_multilayer - post-process an unordered multilayer partition
%   postprocess_ordinal_multilayer     - post-process an ordered multilayer partition
%
%
% Persistence:
%
%   categorical_persistence            - computes the persistence of an unordered multilayer partition
%   ordinal_persistence                - computes the persistence of an ordered multilayer partition
%
%
% Sorting functions for visualization:
%
%   sort_categorical                   - reorders nodes and layers to emphasize persistent structure in an unordered multilayer partition
%   sort_ordinal                       - reorders nodes to emphasize persistent structure in an ordered multilayer partition
