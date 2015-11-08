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
//  Last modified by Lucas Jeub on 25/07/2014


#include "mex.h"


#include "matlab_matrix.h"
#include "group_index.h"
#include <unordered_map>
#include <cstring>

#ifndef OCTAVE
    #include "matrix.h"
#endif

using namespace std;
static group_index group;
static vector<double>* mod_reduced =new vector<double>();
static bool return_sparse;

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
        if (!strcmp(handle, "assign")) {
            if (nrhs!=2) {
                mexErrMsgIdAndTxt("metanetwork_reduce:assign", "assign needs 1 input argument");
            }
            group=prhs[1];
            delete mod_reduced;
            mod_reduced= new vector<double>(group.n_groups,0);
            //zero out mod_reduced for next iteration

            return_sparse=false;
        }
        else if (!strcmp(handle, "reduce")){
            if (nrhs!=2||nlhs!=0) {
                mexErrMsgIdAndTxt("metanetwork_reduce:reduce", "reduce needs 1 input and no output argument");
            }
            if (mxIsDouble(prhs[1])) {
                
                if (mxIsSparse(prhs[1])) {
                    return_sparse=true;
                    sparse mod_s(prhs[1]);
                    
                    if (mod_s.m==group.n_nodes) {

                        for (mwIndex i=0; i<mod_s.nzero(); i++) {
                            mod_reduced->operator[](group.nodes[mod_s.row[i]])+=mod_s.val[i];
                        }

                    }
                    else {
                        mexErrMsgIdAndTxt("metanetwork_reduce:reduce:mod", "input modularity matrix has wrong size");
                    }
                }
                else {
                    full mod_d(prhs[1]);
                    if (mod_d.m==group.n_nodes) {

                        for (mwIndex i=0; i<group.n_groups; i++) {
                            for (list<mwIndex>::iterator it=group.groups[i].begin(); it!=group.groups[i].end(); ++it) {
                                for (full::rowiterator rit=mod_d.rowit(*it); rit!=mod_d.rowit(*it+1); ++rit) {
                                    mod_reduced->operator[](i)+=*rit;
                                }
                            }
                        }

                    }
                    else {
                        mexErrMsgIdAndTxt("metanetwork_reduce:reduce:mod", "input modularity matrix has wrong size");
                    }
                }
            }
        }
        else if (!strcmp(handle, "nodes")) {
            if (nrhs!=2||nlhs<1) {
                mexErrMsgIdAndTxt("metanetwork_reduce:nodes", "nodes needs 1 input and 1 output argument");
            }
            full nodes=group.index(*mxGetPr(prhs[1])-1);
            nodes.export_matlab(plhs[0]);
        }
        else if (!strcmp(handle, "return")) {
            if (nrhs!=1|nlhs!=1) {
                mexErrMsgIdAndTxt("metanetwork_reduce:return", "return needs 1 output argument and no input arguments");
            }
            if (return_sparse) {
                mwSize nmax=0;
                for (vector<double>::iterator it=mod_reduced->begin(); it!=mod_reduced->end(); ++it) {
                    if (*it!=0) {
                        ++nmax;
                    }
                }
                sparse mod_out(group.n_groups,1,nmax);
                mwIndex c=0;
                mwIndex i=0;
                mod_out.col[0]=0;
                for (vector<double>::iterator it=mod_reduced->begin(); it!=mod_reduced->end(); ++it) {
                    if (*it!=0) {
                        mod_out.row[c]=i;
                        mod_out.val[c]=*it;
                        ++c;
                    }
                    ++i;
                }
                mod_out.col[1]=c;
                mod_out.export_matlab(plhs[0]);
            }
            else {
                full mod_out(group.n_groups,1);
                vector<double>::iterator it=mod_reduced->begin();
                for (mwIndex i=0; i<group.n_groups; ++i) {
                    mod_out[i]=*it;
                    ++it;
                }
                mod_out.export_matlab(plhs[0]);
            }
            //zero out mod_reduced for next iteration
            for (vector<double>::iterator it=mod_reduced->begin(); it!=mod_reduced->end(); ++it) {
                *it=0;
            }
            return_sparse=false;
         }
    }
    else {
        mexErrMsgIdAndTxt("metanetwork_reduce:handle", "need a handle to function");
    }
}