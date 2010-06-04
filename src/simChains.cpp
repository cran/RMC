#include"RMC.h"

extern "C" 	void simChains( int *, int *, double *, int *, int *, int *, int *, double *, int *, int *, double *);
int multidraw( int, vector<double>);
extern "C" void testMult( int *, int *, int *, double *);


extern "C" {
	void simChains( int *Rstates, int *Rnobs, double *RX, int *nCats, int *nCovars, int *nTot, int *nChains, double *Rpars, int *RphiID, int *RpiID, double *Rmarg)
{

	GetRNGstate();

	RMCdata dat;
	dat.setVals( *nCats, *nCovars, *nTot, *nChains, RX, Rnobs, Rstates);
	PARdata parINF;
	parINF.setVals( RphiID, RpiID, *nCovars, *nCats);

	vector< AD<double> > params( parINF.nPars, 0);
	for( size_t i=0; i<parINF.nPars; i++)
		params.at( i) = Rpars[i];

	int start = 0;
	int stop = 0;
	vector< AD<double> > piVec( dat.nCats, 0);
	AD<double> phi=0.;
	vector<double> TPcol( dat.nCats, 0.);
	
	for( size_t i=0; i<dat.nChains; i++){
		start = stop;
		stop += dat.nObs.at(i);
		piVec = calcPiVec( params, dat, parINF, 0);
		double sum=0.;
		for( size_t j=0; j<dat.nCats; j++){
			phi = calcSinglePhi( params, dat, parINF, 0, j+1);	
			TPcol.at(j) = Value( piVec.at(j)) / Value( phi);
			sum += TPcol.at( j);
		}
		for( size_t j=0; j<dat.nCats; j++){
			TPcol.at( j) = TPcol.at(j)/sum;
			Rmarg[j] = TPcol.at(j);
		}
		dat.states.at( start) = multidraw( dat.nCats, TPcol);
		for( size_t j=start+1; j<stop; j++){
			piVec = calcPiVec( params, dat, parINF, j);
			phi = calcSinglePhi( params, dat, parINF, j, dat.states.at( j-1) );
			for( size_t k=0; k<dat.nCats; k++)
				TPcol.at(k) = Value( phi) * Value( piVec.at( k));
			TPcol.at( dat.states.at( j-1)-1) += 1-Value( phi);
			dat.states.at( j) = multidraw( dat.nCats, TPcol);		
		}
	}
	for( size_t i=0; i<dat.nTot; i++)
		Rstates[i] = dat.states.at( i);
	
	PutRNGstate();
}
}

int multidraw( int nCats1, vector<double> probs) 
{
	bool flag = true;
	double temp = 0;
	size_t place = 0;
	vector<double> cumprobs;
	cumprobs.push_back( probs.at( 0));
	for( size_t ii=1; ii<probs.size(); ii++)
	cumprobs.push_back( cumprobs.at( ii-1) + probs.at( ii));
	GetRNGstate();
	temp = runif( 0., 1.);
	PutRNGstate();
	while( flag == true){
		if( temp <= cumprobs.at( place))
			return( place+1);
		if( place == nCats1 -1)
			return( nCats1);
		place++;
	}	
}

	
extern "C" {
	void testMult( int *Draws, int *n, int *m, double *Rprobs){
		vector<double> probs;
		for( size_t i=0; i<*m; i++)
			probs.push_back( Rprobs[i]);
		for( size_t i=0; i<*n; i++)
			Draws[i] = multidraw( *m, probs);
	}
}	

