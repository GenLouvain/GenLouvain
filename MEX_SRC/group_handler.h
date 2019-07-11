//
//  group_handler.h
//  group_handler
//
//  Created by Lucas Jeub on 21/11/2012.
//
// Version: 2.2.0
// Date: Thu 11 Jul 2019 12:25:43 CEST

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

#define NUM_TOL 1e-10


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


//move node to most optimal group
template<class M> double move(group_index & g, mwIndex node, const M & mod);

//move node to random group that increases modularity
template<class M> double moverand(group_index & g, mwIndex node, const M & mod);

//move node to random group with probability proportional to increase in modularity
template<class M> double moverandw(group_index & g, mwIndex node, const M & mod);

set_type possible_moves(group_index & g, mwIndex node, const sparse & mod);

set_type possible_moves(group_index & g, mwIndex node, const full & mod);

map_type mod_change(group_index &g, const sparse &mod,set_type & unique_groups,mwIndex current_node);

map_type mod_change(group_index &g, const full & mod, set_type & unique_groups, mwIndex current_node);

move_list positive_moves(set_type & unique_groups, map_type & mod_c);


//implement unique_group_map (quick membership check and insertion of elements, quick iteration over members, unordered)
unique_group_map::unique_group_map() : ismember(std::vector<bool>()) {}
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

//implement template functions
template<class M> double move(group_index & g, mwIndex node, const M & mod){
    set_type unique_groups=possible_moves(g, node, mod);
    map_type mod_c=mod_change(g, mod, unique_groups, node);
    
    //find best move
    double mod_max=0;
    double d_step=0;
    mwIndex group_move=g.nodes[node]; //stay in current group if no improvement
    for(set_type::iterator it=unique_groups.begin();it!=unique_groups.end();++it){
        if(mod_c[*it]>mod_max){
            mod_max=mod_c[*it];
            group_move=*it;
        }
    }
    
    //move current node to most optimal group
    if(mod_max>NUM_TOL){
        g.move(node,group_move);
        d_step+=mod_max;
    }
    return d_step;
}


//set up random engine
std::default_random_engine generator((unsigned int)time(0));

//move node to random group increasing modularity
template<class M> double moverand(group_index & g, mwIndex node, const M & mod){
    set_type unique_groups=possible_moves(g, node, mod);
    map_type mod_c=mod_change(g, mod, unique_groups, node);
    
    //find modularity increasing moves
    move_list mod_pos=positive_moves(unique_groups, mod_c);
    
    // move node to a random group that increases modularity
    double d_step=0;
    if (!mod_pos.first.empty()) {
        std::uniform_int_distribution<mwIndex> randindex(0,mod_pos.first.size()-1);
        mwIndex randmove=randindex(generator);
        g.move(node,mod_pos.first[randmove]);
        d_step=mod_pos.second[randmove];
    }
    return d_step;
}


//move to random group with probability proportional to increase in modularity
template<class M> double moverandw(group_index & g, mwIndex node, const M & mod){
    set_type unique_groups=possible_moves(g, node, mod);
    map_type mod_c=mod_change(g, mod, unique_groups, node);
    
    //find modularity increasing moves
    move_list mod_pos=positive_moves(unique_groups, mod_c);
    
    //move node to a random group that increases modularity with probability proportional to the increase
    double d_step=0;
    if (!mod_pos.first.empty()) {
        std::discrete_distribution<mwIndex> randindex(mod_pos.second.begin(),mod_pos.second.end());
        mwIndex randmove=randindex(generator);
        g.move(node,mod_pos.first[randmove]);
        d_step=mod_pos.second[randmove];
    }
    return d_step;
}


#endif /* defined(__group_handler__group_handler__) */
