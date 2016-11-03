//
//  group_handler.h
//  group_handler
//
//  Created by Lucas Jeub on 21/11/2012.
//
//  Last modified by Lucas Jeub on 25/07/2014

#ifndef __group_handler__group_handler__
#define __group_handler__group_handler__

#include "mex.h"

#ifndef OCTAVE
    #include "matrix.h"
#endif

#include "matlab_matrix.h"
#include "group_index.h"
#include <cstring>
#include <unordered_map>
#include <set>
#include <vector>
#include <random>
#include <ctime>
#include <utility>


typedef std::unordered_map<mwIndex, double> map_type;
//typedef std::map<mwIndex, double> map_type;
//typedef std::vector<double> map_type;

//map for unique possible moves
struct unique_group_map {
    unique_group_map();
    unique_group_map(mwSize n);
    std::vector<bool> ismember;
    std::vector<mwIndex> members;
    
    bool count(mwIndex i);
    void insert(mwIndex i);
    
    typedef std::vector<mwIndex>::iterator iterator;
    iterator begin();
    iterator end();
};

typedef unique_group_map set_type;


typedef std::pair<std::vector<mwIndex>,std::vector<double>> move_list;


//move node to group with most improvement in modularity
double move(group_index & g, mwIndex node, const mxArray * mod);

//move node to random group that increases modularity
double moverand(group_index & g, mwIndex node, const mxArray * mod);

//move node to random group with probability proportional to increase in modularity
double moverandw(group_index & g, mwIndex node, const mxArray * mod);


set_type possible_moves(group_index & g, mwIndex node, sparse & mod);

set_type possible_moves(group_index & g, mwIndex node, full & mode);

map_type mod_change(group_index &g, sparse &mod,set_type & unique_groups,mwIndex current_node);

map_type mod_change(group_index &g, full & mod, set_type & unique_groups, mwIndex current_node);

move_list positive_moves(set_type & unique_groups, map_type & mod_c);

unique_group_map::unique_group_map(mwSize n) : ismember(std::vector<bool>(n,false)) {}

bool unique_group_map::count(mwIndex i) {return ismember[i];}

void unique_group_map::insert(mwIndex i) {
    if (!ismember[i]) {
        ismember[i] = true;
        members.push_back(i);
    }
}

unique_group_map::iterator unique_group_map::begin() { return members.begin(); }

unique_group_map::iterator unique_group_map::end() { return members.end(); }




#endif /* defined(__group_handler__group_handler__) */
