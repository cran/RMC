//Code for calculating the predicted values from the conditional Markov model.

#include"RMC.h"

extern "C" void getCondPreds( double *, double *, int *, int *, int *, int *, double *, int *, int *, int *, int *);

extern "C" {
	void getCondPreds( double *R_probs, double *pars, int *nCats, int *nCovars, int *nTot, int *nChains, double *R_X, int *R_phiID, int *R_nObs, int *R_states, int *R_piID)
{
	RMCdata dat;
	dat.setVals( *nCats, *nCovars, *nTot, *nChains, R_X, R_nObs, R_states);
	PARdata parINF;
	parINF.setVals( R_phiID, R_piID, *nCovars, *nCats);
	
	vector< AD<double> > params( parINF.nPars, 0);
	for( int i=0; i<parINF.nPars; i++)
		params.at( i) = pars[i];
	
	vector< vector<double> > probs;
	vector<double> row( dat.nCats,0) , prevDist( dat.nCats, 0);
	vector< AD<double> > piVec( dat.nCats, 0);
	AD<double> phi=0.;
	
	int start = 0;	//first element of each chain
	int stop = dat.nObs.at( 0);	//last element of each chain
	for( int i=0; i<dat.nChains; i++){	//iterator over chain
	   row.assign( dat.nCats, -9);
		probs.push_back( row);
		for( int j=start+1; j<stop; j++){	//element within chain, not the first...
			row.assign( dat.nCats, -9);
			if( ( dat.states.at( j) != -9) && ( dat.states.at( j-1) != -9)){
				piVec = calcPiVec( params, dat, parINF, j);
				phi = calcSinglePhi( params, dat, parINF, j-1, dat.states.at( j-1));
				for( int k=0; k<dat.nCats; k++){
				  row.at( k) = Value( phi * piVec.at( k));
				  if( k == dat.states.at( j-1))
					 row.at( k) += 1-Value( phi);
				}
			}
			probs.push_back( row);
		}
		if( i < dat.nChains - 1){
		  start += dat.nObs.at( i);
		  stop += dat.nObs.at( i+1);
		}
	}
	for( size_t j=0; j<probs.at(1).size(); j++)//dat.nCats; j++)
		for( size_t i=0; i<probs.size(); i++)//dat.nTot; i++)
			R_probs[j*dat.nTot+i] = probs.at( i).at( j);
}
}

