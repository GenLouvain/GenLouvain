//
//  multiord_bipartite.cpp
//  multiord_bipartite
//
//  Created by Lucas Jeub on 30/11/2012.
//
//

#include "matlab_matrix.h"
#include "mex.h"
#include "matrix.h"
#include <cmath>
#include <vector>

using namespace std;

//b=multiord_bipartite(i,A,k,d,mm,omega,gamma)
void mexFunction(int nlhs, mxArray * plhs[], int nrhs,const mxArray * prhs[]){
    //error checking
    if (nrhs!=7) {
        mexErrMsgIdAndTxt("multiord_bipartite:input", "7 input arguments required");
    }
        
    //get input arguments
    mwIndex i=*mxGetPr(prhs[0])-1; //C-index
    full k(prhs[2]);
    full d(prhs[3]);
    full mm(prhs[4]);
    double omega=*mxGetPr(prhs[5]);
    double gamma=*mxGetPr(prhs[6]);

    //sizes
    mwSize T=k.n;
    mwSize m=k.m;
    mwSize n=d.m;
    mwSize N=m+n;
    
    mwIndex s=i/N;
    mwIndex ii=i-s*N;
    
    double twom=mm.get(s);
    
    //get slice of matrix
    sparse A(mxGetCell(prhs[1],s));
    
    if (ii<m) {
        vector<double> modv(n,0);
        sparse modularity(N*T, 1, n+2);
        
        double kk=k.get(ii,s);
        for (mwIndex it = 0; it<n; it++) {
            modv[it]=-gamma*kk*d.get(it, s)/twom;
        }
        
        for (mwIndex it =0; it<A.n; it++) {
            mwIndex it2=A.col[it];
            while((A.row[it2]<=ii)&&it2<A.col[it+1]){
                if (A.row[it2]==ii) {
                    modv[it]+=A.val[it2];
                }
                it2++;
            }
        }
        
        mwIndex count=0;
        //interslice coupling
        if (N<=i) {
            modularity.row[0]=i-N;
            modularity.val[0]=omega;
            count++;
        }
        //intraslice coupling
        for (mwIndex it=0; it<n; it++) {
            if (modv[it]!=0) {
                modularity.row[count]=m+s*N+it;
                modularity.val[count]=modv[it];
                count++;
            }
        }
        //interslice coupling
        if ((N+i)<(N*T)) {
            modularity.row[count]=N+i;
            modularity.val[count]=omega;
            count++;
        }
        
        modularity.col[0]=0;
        modularity.col[1]=count;
        
        if (nlhs>0) {
                  modularity.export_matlab(plhs[0]);
        }
    }
    else {
        vector<double> modv(m,0);
        sparse modularity(N*T,1,m+2);
        
        double dd=d.get(ii-m,s);
        for (mwIndex it=0; it<m; it++) {
            modv[it]=-gamma*k.get(it,s)*dd/twom;
        }
        
        for (mwIndex it=A.col[ii-m]; it<A.col[ii-m+1]; it++) {
            modv[A.row[it]]+=A.val[it];
        }
        
        mwIndex count=0;
        //interslice coupling
        if (N<=i) {
            modularity.row[0]=i-N;
            modularity.val[0]=omega;
            count++;
        }
        
        //intraslice coupling
        for (mwIndex it=0; it<m; it++) {
            if (modv[it]!=0) {
                modularity.row[count]=it+s*N;
                modularity.val[count]=modv[it];
                count++;
            }
        }
        
        //interslice coupling
        if ((N+i)<(N*T)) {
            modularity.row[count]=N+i;
            modularity.val[count]=omega;
            count++;
        }
        
        modularity.col[0]=0;
        modularity.col[1]=count;
        
        if (nlhs>0) {
            modularity.export_matlab(plhs[0]);
        }
    }
}
