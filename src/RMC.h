//Header for Reversible Markov Chain estimation.

#include<vector>
#include<algorithm>

#include<R.h>
#include<Rmath.h>
#include<cppad/cppad.hpp> 
#include<cppad/lu_solve.hpp>

#include "ap.h"
#include "lbfgs.h"

using CppAD::AD;           // use AD as abbreviation for CppAD::AD
using std::vector;         // use vector as abbreviation for std::vector
using CppAD::Value;

//////////////////////////////////////////
//Classes
//////////////////////////////////////////

class RMCdata {
	public:
		RMCdata();
		~RMCdata();
		void setVals( int, int, int, int, double *, int *, int *);
//	private:
		int nCats, nCovars, nTot, nChains;
		vector< vector<double> > X;
//		vector< vector<int> > t_phiID;
		vector<int> nObs, states;//, piID;
};

class PARdata {
	public:
		PARdata();
		~PARdata();
		void setVals( int *, int *, int, int);
		vector< vector<int> > phiIndex;  //note that this is transposed to what gets passed in
		vector<int> piIndex;
		int nPhi, nPi, nPars;
};

///////////////////////////////////////////
//Headers
///////////////////////////////////////////

extern "C" void topTailChains( int *, int *, int *, int *);
extern "C" void RMC( double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, double *, double *, double *, int *, int *, double *);
vector< AD<double> > ploglike( const vector< AD<double> > &, const RMCdata &, const PARdata &, double );
AD<double> chainploglike( const vector< AD<double> > &, const RMCdata &, const PARdata &, int, int);
vector< AD<double> > calcPiVec( const vector< AD<double> > &, const RMCdata &, const PARdata &, int);
AD<double> calcSinglePhi( const vector< AD<double> > &, const RMCdata &, const PARdata &, int, int);
AD<double> first_logl( const vector< AD<double> > &, const RMCdata &, const PARdata &, int);

void funcgrad(ap::real_1d_array, double &, ap::real_1d_array &);


