//Code for filling in the missing covariate values for the second and susequent stages of the RMC analysis

#include"RMC.h"

extern "C" void MVfillCpp( double *, double *, int *, int *, int *, int *, double *, int *, int *, int *, int *);

extern "C" {
	void MVfillCpp( double *R_probs, double *pars, int *nCats, int *nCovars, int *nTot, int *nChains, double *R_X, int *R_phiID, int *R_nObs, int *R_states, int *R_piID)
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
	vector< AD<double> > piVec( dat.nCats, 0), phiVec( dat.nCats, 0);
	AD<double> phi=0.;
	double Tpel=0;
	
	int start = 0;
	int stop = dat.nObs.at( 0);
	for( int i=0; i<dat.nChains; i++){
		for( int j=start; j<stop; j++){
			row.assign( dat.nCats, 0.);		
			if( dat.states.at( j) != -9)			//current observation IS NOT -9		
				row.at( dat.states.at( j) -1) = 1;
			else{									//current observation IS -9
				piVec = calcPiVec( params, dat, parINF, j);			//used to be j-1
				if( dat.states.at( j-1) == -9){		//Internal -9
					for( int k=0; k<dat.nCats; k++)
						phiVec.at( k) = calcSinglePhi( params, dat, parINF, j-1, k+1);		//used to be j-1
					for( int k=0; k<dat.nCats; k++){
//						row.at( k) = 0;
						Tpel = 0;
						for( int l=0; l<dat.nCats; l++){
							Tpel = Value( phiVec.at( l) * piVec.at( k));
							if( k == l)
								Tpel += 1-Value( phiVec.at( l));
							row.at( k) += Tpel * prevDist.at( l);
						}
					}
					for( int k=0; k<dat.nCats; k++)
						prevDist.at( k) = row.at( k);
				}
				else{								//First -9 in block
					phi = calcSinglePhi( params, dat, parINF, j, dat.states.at( j-1));		//used to be j-1
					for( int k=0; k<dat.nCats; k++){
						row.at( k) = Value( piVec.at( k) * phi);
						if( dat.states.at( j-1)-1 == k)
							row.at( k) += Value( 1-phi);
						prevDist.at( k) = row.at( k);
					}
				}
			}
			probs.push_back( row);
		}
		if( i < dat.nChains - 1){
			start += dat.nObs.at( i);
			stop += dat.nObs.at( i+1);		
		}
	}

	for( int j=0; j<dat.nCats; j++)
		for( int i=0; i<dat.nTot; i++)
			R_probs[j*dat.nTot+i] = probs.at( i).at( j);
}
}
	

