//
//  matlab_matrix.h
//  matlab_matrix
//
//  Created by Lucas Jeub on 24/10/2012
//
//  Implements thin wrapper classes for full and sparse matlab matrices
//
//
//  Last modified by Lucas Jeub on 25/07/2014






#ifndef MATLAB_MATRIX_H
#define MATLAB_MATRIX_H

#define inf std::numeric_limits<double>::infinity();


#include <limits>
#include <iterator>

#include "mex.h"

#ifndef OCTAVE
    #include "matrix.h"
#endif

struct full;

struct sparse{
	sparse();
	sparse(mwSize m, mwSize n, mwSize nmax);
	sparse(const sparse &matrix);
	sparse(const mxArray *matrix);

	~sparse();
	
	sparse & operator = (const sparse & matrix);
    
    sparse & operator = (const full & matrix);
	
	sparse & operator = (const mxArray *matrix);
	
	
	/*operations*/
	/*pointwise division*/

	sparse operator / (const sparse & B);
	sparse operator / (const full & B);

	mwSize nzero() const { return col[n];}
    
    double get(mwIndex i, mwIndex j);
	
	void export_matlab(mxArray * & out);
	
	mwSize m;
	mwSize n;
	mwSize nmax;
	mwIndex *row;
	mwIndex *col;
	double *val;

	private:
	
	bool export_flag;
};


struct full{
	full();
	full(mwSize m, mwSize n);
	full(const full &matrix);
	full(const mxArray * matrix);
	
	~full();
    
	void export_matlab(mxArray * & out);
	
	full & operator = (const full & matrix);
    
    full & operator = (const sparse & matrix);
	
	full & operator = (const mxArray * matrix);
	
	double & get(mwIndex i, mwIndex j);
    double get(mwIndex i,mwIndex j) const;
	double & get(mwIndex i);
    double & operator [] (mwIndex i);
    double get(mwIndex i) const;
    double operator [] (mwIndex i) const;
    
	
	full operator / (const sparse &B);
	full operator / (const full &B);

	mwSize m;
	mwSize n;
	
	double *val;
	
private:
	
	bool export_flag;
    
public:
    class rowiterator : public std::iterator<std::random_access_iterator_tag, double, mwSignedIndex> {
        double * p;
        mwSize m;
        mwSize n;
        mwSignedIndex rowpos;
    public:
        rowiterator() : p(nullptr), m(0) {}
        rowiterator(double * init, mwSize _m, mwSize _n, mwSignedIndex _rowpos) : p(init), m(_m), n(_n) ,rowpos(_rowpos) {}
        rowiterator & operator=(const rowiterator & it) {p=it.p; m=it.m; n=it.n; rowpos=it.rowpos; return *this;}
        bool operator == (const rowiterator & it) {return rowpos==it.rowpos;}
        bool operator != (const rowiterator & it) {return rowpos!=it.rowpos;}
        bool operator < (const rowiterator & it) {return rowpos<it.rowpos;}
            bool operator <= (const rowiterator & it) {return rowpos<=it.rowpos;}
        bool operator > (const rowiterator & it) {return rowpos>it.rowpos;}
        bool operator >= (const rowiterator & it) {return rowpos>=it.rowpos;}
        double & operator*();
        rowiterator & operator ++ (){rowpos++;return *this;}
        rowiterator operator ++ (int) {rowiterator tmp(*this); operator++(); return tmp;}
        rowiterator & operator -- () {rowpos--;return *this;}
        rowiterator operator -- (int) {rowiterator tmp(*this); operator--(); return tmp;}
        rowiterator & operator + (mwSignedIndex i) {rowpos+=i; return *this;}
        mwSignedIndex operator + (const rowiterator & it) {return it.rowpos+rowpos;}
        rowiterator & operator - (mwSignedIndex i) {rowpos-=i; return *this;}
        mwSignedIndex operator - (const rowiterator & it) {return it.rowpos-rowpos;}
        double & operator [] (mwSignedIndex i);
        double operator [] (mwSignedIndex i) const;
    };
    
    rowiterator rowit(mwIndex i);
    rowiterator rowit(mwIndex i,mwIndex j);
    
    typedef double * coliterator;
    
    coliterator colit(mwIndex i) {return val+i;}
    coliterator colit(mwIndex i,mwIndex j) {return val+(i+j*m);}
};


#endif
