//Predicting stationary dist at a single point for first factor of the model

#include"RMC.h"

extern "C" void predPtRMC( double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, double *, int *);
vector< AD<double> > statDist( const vector< AD<double> > &, const RMCdata &, const PARdata &, int);


vector< AD<double> > statDist( const vector< AD<double> > &params5, const RMCdata &dat5, const PARdata &par_info5, int index5)
{
	vector< AD<double> > phiVec5( dat5.nCats, 0);
	vector< AD<double> > piVec5( dat5.nCats, 0);
	vector< AD<double> > dist( dat5.nCats, 0);
	AD<double> sum = 0.;
	
	piVec5 = calcPiVec( params5, dat5, par_info5, index5);
	for( size_t j=0; j<dat5.nCats; j++)
		phiVec5.at( j) = calcSinglePhi( params5, dat5, par_info5, index5, j+1);
	
	for( size_t i=0; i<dat5.nCats; i++){
		dist.at( i) = piVec5.at( i) / phiVec5.at( i);
		sum += dist.at( i);
	}
	for( size_t i=0; i<dat5.nCats; i++)
		dist.at( i) = dist.at( i) / sum;

	return( dist);	
}



extern "C" {
	void predPtRMC( double *Rpars, double *Rpts, int *Rncats, int *Rncovars, int *Rntot, int *RphiID, int *RpiID, int * RdumbStates, double *ppts, double *area, double *Rjac, int *RdoCov)
{
	//dummy variables for initialising (re-use previous code)
	int intDumb[1] = {87};

	RMCdata predDat;
	predDat.setVals(*Rncats, *Rncovars, *Rntot, 1, Rpts, intDumb, RdumbStates);
	PARdata predPar_info;
	predPar_info.setVals( RphiID, RpiID, *Rncovars, *Rncats);

	//Parameter vector - independent
	vector< AD<double> > params( predPar_info.nPars);
	vector<double> predVal;
	for( size_t i=0; i<predPar_info.nPars; i++){
		params.at( i) = Rpars[i];
		predVal.push_back( Rpars[i]);
	}
	
	if( *RdoCov==1)
		CppAD::Independent( params);
					
	//prediction matrix - dependent (row wise)
	vector< vector< AD<double> > > dist;

	for( size_t i=0; i<predDat.nTot; i++)
		dist.push_back( statDist( params, predDat, predPar_info, i));
		
	AD<double> sum=0;
	vector< AD<double> > areaPred;
	
	for( size_t i=0; i<predDat.nCats; i++){
		sum = 0;
		for( size_t j=0; j<predDat.nTot; j++)
			sum += dist.at(j).at(i);
		areaPred.push_back( sum / predDat.nTot);
	}

	vector<double> jac(params.size() * predDat.nCats);	
	if( *RdoCov==1){
		CppAD::ADFun<double> fpred(params, areaPred);
		jac = fpred.Jacobian(predVal);
	}

	//"output" array
	for( size_t i=0; i<predDat.nTot; i++)
		for( size_t j=0; j<predDat.nCats; j++)
			ppts[i+j*predDat.nTot] = Value( dist.at(i).at(j));
	for( size_t i=0; i<predDat.nCats; i++)
		area[i] = Value( areaPred.at(i));
	if( *RdoCov==1){
		for( size_t i=0; i<predDat.nCats; i++)
			for( size_t j=0; j<params.size(); j++)
				Rjac[j*predDat.nCats+i] = jac.at( i*params.size()+j);
	}


}
}

