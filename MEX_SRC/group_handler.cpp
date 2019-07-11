//
//  group_handler.cpp
//  group_handler
//
//  Created by Lucas Jeub on 21/11/2012.
//
// usage:
//
//  [output]=group_handler('function_handle',input)
//
//  implemented functions are 'assign', 'move', 'moverand', 'return'
//
//      assign: takes a group vector as input and uses it to initialise the "group_index"
//
//
//      move:   takes a node index and the corresponding column of the modularity matrix as
//              input
//
//              moves node to group with maximum improvement in modularity (stays in same group
//              if no improvement possible)
//
//              returns improvement if given an output argument
//
//
//      moverand:   takes a node index and the corresponding column of the modularity matrix as
//              input
//
//              moves node to random group with improvement in modularity (stays in same group
//              if no improvement possible). Chooses from possible moves uniformly at random.
//
//              returns improvement if given an output argument
//
//      moverandw:   takes a node index and the corresponding column of the modularity matrix as
//              input
//
//              moves node to random group with improvement in modularity (stays in same group
//              if no improvement possible). Chooses from possible moves with probability
//              proportional to increase in modularity.
//
//              returns improvement if given an output argument
//
//
//      return: outputs the community assignment for all nodes as a tidy group vector, that is
//              e.g. S = [1 2 1 3] rather than S = [3 1 3 2]
//
//
// Version: 2.2.0
// Date: Thu 11 Jul 2019 12:25:43 CEST


#include "group_handler.h"

using namespace std;

static group_index group;
//switch on handle
enum func {ASSIGN, MOVE, MOVERAND, MOVERANDW, RETURN};
static const unordered_map<string, func> function_switch({ {"assign", ASSIGN}, {"move", MOVE}, {"moverand", MOVERAND}, {"moverandw", MOVERANDW}, {"return", RETURN} });

//group_handler(handle, varargin)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    if (nrhs>0) {
        //get handle to function to perform
        mwSize strleng = mxGetM(prhs[0])*mxGetN(prhs[0])+1;
        char * handle;
        handle=(char *) mxCalloc(strleng, sizeof(char));
        
        if (mxGetString(prhs[0],handle,strleng)) {
            mexErrMsgIdAndTxt("group_handler:handle", "handle needs to be a string");
        }
                
        if (function_switch.count(handle)) {
            switch (function_switch.at(handle)) {
                case ASSIGN: {
                    if (nrhs!=2) {
                        mexErrMsgIdAndTxt("group_handler:assign", "assign needs 1 input argument");
                    }
                    group=prhs[1];
                    break;
                }
                    
                case MOVE: {
                    if (nrhs !=3) {
                        mexErrMsgIdAndTxt("group_handler:move", "move needs 2 input arguments");
                    }
                    double dstep;
                    mwIndex node=((mwIndex) * mxGetPr(prhs[1]))-1;
                    if (mxIsSparse(prhs[2])) {
                        sparse mod_s(prhs[2]);
                        dstep = move(group, node, mod_s);
                    } else {
                        full mod_d(prhs[2]);
                        dstep = move(group, node, mod_d);
                    }
                    //output improvement in modularity
                    if (nlhs>0) {
                        plhs[0]=mxCreateDoubleScalar(dstep);
                    }
                    break;
                }
                    
                case MOVERAND: {
                    if (nrhs!=3) {
                        mexErrMsgIdAndTxt("group_handler:moverand", "move needs 2 input arguments");
                    }
                    double dstep;
                    mwIndex node=((mwIndex) * mxGetPr(prhs[1]))-1;
                    if (mxIsSparse(prhs[2])) {
                        sparse mod_s(prhs[2]);
                        dstep = moverand(group, node, mod_s);
                    } else {
                        full mod_d(prhs[2]);
                        dstep = moverand(group, node, mod_d);
                    }
                    
                    //output improvement in modularity
                    if (nlhs>0) {
                        plhs[0]=mxCreateDoubleScalar(dstep);
                    }
                    break;
                }
                    
                case MOVERANDW: {
                    if (nrhs!=3) {
                        mexErrMsgIdAndTxt("group_handler:moverandw", "move needs 2 input arguments");
                    }
                    double dstep;
                    mwIndex node=((mwIndex) * mxGetPr(prhs[1]))-1;
                    if (mxIsSparse(prhs[2])) {
                        sparse mod_s(prhs[2]);
                        dstep = moverandw(group, node, mod_s);
                    } else {
                        full mod_d(prhs[2]);
                        dstep = moverandw(group, node, mod_d);
                    }
                    
                    //output improvement in modularity
                    if (nlhs>0) {
                        plhs[0]=mxCreateDoubleScalar(dstep);
                    }
                    break;
                }
                    
                case RETURN: {
                    if (nlhs>0) {
                        group.export_matlab(plhs[0]);
                    }
                    else {
                        mexErrMsgIdAndTxt("group_handler:return", "need ouput argument to return");
                    }
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
        mexErrMsgIdAndTxt("group_handler:input", "need handle to function");
    }
}


//find possible moves
set_type possible_moves(group_index & g, mwIndex node, const sparse & mod){
    set_type unique_groups(g.n_groups);
    unique_groups.insert(g.nodes[node]);
    //add nodes with potential positive contribution to unique_groups
    for(mwIndex i=0; i<mod.nzero(); ++i){
        if(mod.val[i]>0){
            unique_groups.insert(g.nodes[mod.row[i]]);
        }
    }
    return unique_groups;
}


set_type possible_moves(group_index & g, mwIndex node, const full & mod){
    set_type unique_groups(g.n_groups);
    unique_groups.insert(g.nodes[node]);
    //add nodes with potential positive contribution to unique_groups
    for(mwIndex i=0; i<g.n_nodes; ++i){
        if(mod.get(i)>0){
            unique_groups.insert(g.nodes[i]);
        }
    }
    return unique_groups;
}


//calculates changes in modularity for full modularity matrix
map_type mod_change(group_index &g, const full & mod, set_type & unique_groups, mwIndex current_node){
    mwIndex current_group=g.nodes[current_node];
    map_type mod_c;
    double mod_current= mod[current_node];
    for (set_type::iterator it1=unique_groups.begin(); it1!=unique_groups.end(); ++it1) {
        double mod_c_group=0;
        for(list<mwIndex>::iterator it2=g.groups[*it1].begin(); it2!=g.groups[*it1].end(); ++it2){
            mod_c_group+=mod.get(*it2);
        }
        mod_c[*it1]=mod_c_group;
    }
    mod_c[current_group]-=mod_current;
    mod_current=mod_c[current_group];
    for (set_type::iterator it=unique_groups.begin(); it!=unique_groups.end(); ++it) {
        mod_c[*it]-=mod_current;
    }
    
    return mod_c;
}


//calculates changes in modularity for sparse modularity matrix
map_type mod_change(group_index & g, const sparse & mod, set_type & unique_groups, mwIndex current_node){
    mwIndex current_group=g.nodes[current_node];
    map_type mod_c;
    double mod_current=mod.get(current_node, 0);
    for(set_type::iterator it=unique_groups.begin(); it!=unique_groups.end();++it){
        mod_c[*it]=0;
    }
    //calculate changes in modularity
    for (mwIndex i=0; i<mod.nzero(); ++i) {
        if (unique_groups.count(g.nodes[mod.row[i]])) {
            mod_c[g.nodes[mod.row[i]]]+=mod.val[i];
        }
    }
    mod_c[current_group]-=mod_current;
    mod_current=mod_c[current_group];
    for (set_type::iterator it=unique_groups.begin(); it!=unique_groups.end(); ++it) {
        mod_c[*it]-=mod_current;
    }
    return mod_c;
}

//find moves that improve modularity
move_list positive_moves(set_type & unique_groups, map_type & mod_c){
    move_list moves;
    for(set_type::iterator it=unique_groups.begin();it!=unique_groups.end();++it){
        if(mod_c[*it]>NUM_TOL){
            moves.first.push_back(*it);
            moves.second.push_back(mod_c[*it]);
        }
    }
    return moves;
}


