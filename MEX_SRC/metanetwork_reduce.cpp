//
//  metanetwork_reduce.cpp
//  metanetwork_reduce
//
//  Created by Lucas Jeub on 22/11/2012.
//
// usage:
//
//  [output]=metanetwork_reduce('function_handle',input)
//
//  implemented functions are 'assign', 'move', 'return'
//
//      assign: takes a group vector as input and uses it to initialise the "group_index"
//
//
//      reduce: takes a column of the modularity matrix as
//              input
//
//              returns reduced column, where the i's entry is the sum of the original modularity
//              matrix over all nodes in group i
//
//
//      nodes: takes a group and returns the matlab index of all nodes in this group
//
//
// Version: 2.2.0
// Date: Thu 11 Jul 2019 12:25:43 CEST


#include "mex.h"

#include "matlab_matrix.h"
#include "group_index.h"
#include <unordered_map>
#include <cstring>
#include <string>

#ifndef OCTAVE
    #include "matrix.h"
#endif

using namespace std;
static group_index group;
static vector<double> mod_reduced=vector<double>();
static bool return_sparse;

enum func {ASSIGN, REDUCE, NODES, RETURN};
static const unordered_map<string, func> function_switch({ {"assign", ASSIGN}, {"reduce", REDUCE}, {"nodes", NODES}, {"return", RETURN} });

//metanetwork_reduce(handle, varargin)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if (nrhs>0) {
        //get handle to function to perform
        mwSize strleng = mxGetM(prhs[0])*mxGetN(prhs[0])+1;
        char * handle;
        handle=(char *) mxCalloc(strleng, sizeof(char));
        
        if (mxGetString(prhs[0],handle,strleng)) {
            mexErrMsgIdAndTxt("metanetwork_reduce:handle:string", "handle needs to be a string");
        }
        


        //switch on handle
        if (function_switch.count(handle)>0) {
            switch (function_switch.at(handle)) {
                    
                case ASSIGN: {
                    //assign new group structure for aggregation
                    if (nrhs!=2) {
                        mexErrMsgIdAndTxt("metanetwork_reduce:assign", "assign needs 1 input argument");
                    }
                    group=prhs[1];
                    //zero out mod_reduced for next iteration
                    mod_reduced=vector<double>(group.n_groups,0);
                    return_sparse=false;
                    break;
                }
                    
                case REDUCE: {
                    //add modularity contributions to current column
                    if (nrhs!=2||nlhs!=0) {
                        mexErrMsgIdAndTxt("metanetwork_reduce:reduce", "reduce needs 1 input and no output argument");
                    }
                    if (mxIsDouble(prhs[1])) {
                        if (mxIsSparse(prhs[1])) { //sparse modularity input
                            return_sparse=true;
                            sparse mod_s(prhs[1]);
                            if (mod_s.m==group.n_nodes) {
                                for (mwIndex i=0; i<mod_s.nzero(); i++) {
                                    mod_reduced[group.nodes[mod_s.row[i]]]+=mod_s.val[i];
                                }
                            }
                            else {
                                mexErrMsgIdAndTxt("metanetwork_reduce:reduce:mod", "input modularity matrix has wrong size");
                            }
                        }
                        else {//full modularity input
                            full mod_d(prhs[1]);
                            if (mod_d.m==group.n_nodes) {
                                for (mwIndex i=0; i<group.n_groups; i++) {
                                    for (list<mwIndex>::iterator it=group.groups[i].begin(); it!=group.groups[i].end(); ++it) {
                                        for (full::rowiterator rit=mod_d.rowit(*it); rit!=mod_d.rowit(*it+1); ++rit) {
                                            mod_reduced[i]+=*rit;
                                        }
                                    }
                                }
                            }
                            else {
                                mexErrMsgIdAndTxt("metanetwork_reduce:reduce:mod", "input modularity matrix has wrong size");
                            }
                        }
                    }
                    break;
                }
                    
                case NODES: {
                    //return matlab indeces of nodes in group i
                    if (nrhs!=2||nlhs<1) {
                        mexErrMsgIdAndTxt("metanetwork_reduce:nodes", "nodes needs 1 input and 1 output argument");
                    }
                    full nodes=group.index(*mxGetPr(prhs[1])-1);
                    nodes.export_matlab(plhs[0]);
                    break;
                }
                    
                case RETURN: {
                    //return modularity contributions and reset
                    if (nrhs!=1|nlhs!=1) {
                        mexErrMsgIdAndTxt("metanetwork_reduce:return", "return needs 1 output argument and no input arguments");
                    }
                    if (return_sparse) {
                        sparse mod_out=mod_reduced;
                        mod_out.export_matlab(plhs[0]);
                    }
                    else {
                        full mod_out=mod_reduced;
                        mod_out.export_matlab(plhs[0]);
                    }
                    //zero out mod_reduced for next iteration
                    for (vector<double>::iterator it=mod_reduced.begin(); it!=mod_reduced.end(); ++it) {
                        *it=0;
                    }
                    return_sparse=false;
                    break;
                }
                    
                default: {
                    mexErrMsgIdAndTxt("metanetwork_reduce:switch","switch implementation error");
                    break;
                }
            }
        } else {
            mexErrMsgIdAndTxt("group_handler:handle", "invalid handle");
        }
    } else {
        mexErrMsgIdAndTxt("metanetwork_reduce:handle", "need a handle to function");
    }
}
