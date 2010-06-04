 //Script for estimating parameters for the Reversible Markov Chain
//Employs the lBFGS Quasi Newton routine from ALGLIB

#include"RMC.h"

//////////////////////////////////////////
//Naughty naughty global variables
//////////////////////////////////////////

CppAD::ADFun<double> F;

//////////////////////////////////////////
//Function definintions
//////////////////////////////////////////

extern "C" { 
	void topTailChains( int *c_states9, int *c_nobs, int *c_nnobs, int *idvec)
	{
		int ii, tester, place=0;
		
		for( ii = 0; ii < *c_nnobs; ii++){
			if( ii!=0)
				place += c_nobs[ii-1];
			tester = place;
			while( c_states9[tester] == -9){
				idvec[tester] = 1;
				tester += 1;
			}
			tester = place + c_nobs[ii] -1;
			while( c_states9[tester] == -9){
				idvec[tester] = 1;
				tester -= 1;
			}
		}
	}
}

extern "C" {
	void RMC( double *pars, double *R_plogl, double *R_scores, double *R_hessian, int *R_calcHes, int *R_conv, int *R_iters, int *nCats, int *nCovars, int *nTot, int *nChains, double *R_X, int *R_phiID, int *R_nObs, int *R_states, int *R_piID, int *m, double *epsg, double *epsf, double *epsx, int *max_iter, int *quiet, double *penalty)
	{
		RMCdata dat;
		dat.setVals( *nCats, *nCovars, *nTot, *nChains, R_X, R_nObs, R_states);
		PARdata par_info;
		par_info.setVals( R_phiID, R_piID, *nCovars, *nCats);
//		nPars = par_info.nPars;								//nPars is global
		
		//Independent variables
		vector< AD<double> > params( par_info.nPars);
		for( size_t i=0; i<par_info.nPars; i++)
			params.at( i) = pars[i];
		CppAD::Independent( params);
		
		//Dependent parameters
		vector< AD<double> > pll( 1, 0.);
	
		//calculating ploglikelihood
		pll = ploglike( params, dat, par_info, *penalty);
			
		//Definition of the function
		F.Dependent( pll);									//F is a global ADFun 	
		
		ap::real_1d_array x;
		x.setcontent( 1, par_info.nPars, pars);
		
		if( *quiet==1)
			lbfgsminimize( par_info.nPars, *m, x, *epsg, *epsf, *epsx, *max_iter, *R_conv, *quiet);
		else
			lbfgsminimize( par_info.nPars, *m, x, *epsg, *epsf, *epsx, *max_iter, *R_conv);
		
		vector<double> esti;
		for( size_t i=0; i<par_info.nPars; i++){
			esti.push_back( x(i+1));
			pars[i] = x( i+1);
		}		

		*R_plogl = F.Forward( 0, esti).at( 0);
		vector<double> tempScore( par_info.nPars, -99);
		tempScore = F.RevOne( esti, 0);
		vector<double> temphes( par_info.nPars*par_info.nPars, 0);
		if( *R_calcHes==1)
			temphes = F.Hessian( esti, 0);

		for( size_t i=0; i<par_info.nPars; i++){
			R_scores[i] = tempScore.at( i);
			for( size_t j=0; j<par_info.nPars; j++)
				R_hessian[j*par_info.nPars+i] = temphes.at( j*par_info.nPars+i);
		}
	
	}
}

void funcgrad(ap::real_1d_array x, double &f, ap::real_1d_array &g)
{
	vector<double> newX;
	vector<double> ftemp( 1);
	int dim = x.gethighbound() - x.getlowbound() + 1;
	vector<double> gtemp( dim);
	
	
	for( size_t i=1; i<=dim; i++)
		newX.push_back( x(i));
		
	ftemp = F.Forward( 0, newX);
	gtemp = F.RevOne( newX, 0);
	
	f = 0 - ftemp.at( 0);
	for( size_t i=0; i<dim; i++)
		g(i+1) = 0 - gtemp.at( i);
}


vector< AD<double> > ploglike( const vector< AD<double> > &params1, const RMCdata &dat1, const PARdata &par_info1, double penalty1)
{
	vector< AD<double> > templl( dat1.nChains, 0.), pll1( 1, 0.);
	
//	pll1.push_back( 0.);
	
	//Calculating the likelihood for each chain
	int start = 0;
	AD<double> penTerm = 0.0;
	int stop = dat1.nObs.at( 0);
	for( size_t i=0; i<dat1.nChains; i++){
		templl.at( i) = 0;
		templl.at( i) = chainploglike( params1, dat1, par_info1, start, stop);

		if ( i < dat1.nChains-1){
			start += dat1.nObs.at( i);
			stop += dat1.nObs.at( i+1);
		}
	}
	for( size_t i=0; i<dat1.nChains; i++)
		pll1.at( 0) += templl.at( i);
	if( penalty1>0) {
		for( size_t i=0; i<par_info1.nPars; i++)
			penTerm += params1.at( i)*params1.at( i) / ( 2 * penalty1);
		pll1.at( 0) -= penTerm;
//		Rprintf( "%f", (double)penTerm);
	}
	return( pll1);
}

AD<double> first_logl( const vector< AD<double> > &params5, const RMCdata &dat5, const PARdata &par_info5, int index5)
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
	dist.at( dat5.states.at( index5)-1) = dist.at( dat5.states.at( index5)-1) / sum;
	return( log( dist.at( dat5.states.at( index5)-1)));
}

AD<double> chainploglike( const vector< AD<double> > &params2, const RMCdata &dat2, const PARdata &par_info2, int start2, int stop2)
{
	AD<double> pllsum=0., phi=0., temp=0., Tpel=0.;
	vector< AD<double> > piVec( dat2.nCats, 0), phiVec( dat2.nCats, 0), prevDist( dat2.nCats, 0), newDist( dat2.nCats, 0);
	
	pllsum += first_logl( params2, dat2, par_info2, start2);
	for( size_t ii=start2+1; ii<stop2; ii++){
		temp = 0;
		piVec = calcPiVec( params2, dat2, par_info2, ii);		//used to be ii-1
		if( dat2.states.at( ii) != -9){							//current observation IS NOT -9
			if( dat2.states.at( ii-1) != -9){					//Not first observation after -9
				phi = calcSinglePhi( params2, dat2, par_info2, ii, dat2.states.at( ii-1));		//used to be ii-1
				temp = phi * piVec.at( dat2.states.at( ii)-1);
				if( dat2.states.at( ii) == dat2.states.at( ii-1))
					temp += (1-phi);
				pllsum += log( temp);
			}
			else{																			//First observation after -9
				temp = 0;
				for( size_t jj=0; jj<dat2.nCats; jj++){
					phi = calcSinglePhi( params2, dat2, par_info2, ii, jj+1);		//used to be ii-1
					if( jj == dat2.states.at( ii)-1)
						temp += ( ( 1-phi) + phi*piVec.at( dat2.states.at( ii)-1)) * prevDist.at( jj);
					else
						temp += phi*piVec.at( dat2.states.at( ii)-1)*prevDist.at( jj);
				}
				pllsum += log( temp);
			}
		}
		else{																				//current observation IS -9
			if( dat2.states.at( ii-1) == -9){									//Internal -9
				for( size_t jj=0; jj<dat2.nCats; jj++)
					phiVec.at( jj) = calcSinglePhi( params2, dat2, par_info2, ii, jj+1);		//used to be ii-1
				for( size_t jj=0; jj<dat2.nCats; jj++){
					newDist.at( jj) = 0;
					Tpel = 0;
					for( size_t kk=0; kk<dat2.nCats; kk++){
						Tpel = phiVec.at( kk) * piVec.at( jj);
						if( jj == kk)
							Tpel += 1-phiVec.at( kk);
						newDist.at( jj) += Tpel * prevDist.at( kk);
					}
				}
				for( size_t jj=0; jj<dat2.nCats; jj++)
					prevDist.at( jj) = newDist.at( jj);
			}
			else{																	//First -9 in a block
				phi = calcSinglePhi( params2, dat2, par_info2, ii, dat2.states.at( ii-1));		//used to be ii-1
				for( size_t jj=0; jj<dat2.nCats; jj++){
					prevDist.at( jj) = piVec.at( jj) * phi;
					if( dat2.states.at( ii-1)-1 == jj)
						prevDist.at( jj) += 1-phi;
				}
			}
		}
	}
	return( pllsum);

}
			
vector< AD<double> > calcPiVec( const vector< AD<double> > &params3, const RMCdata &dat3, const PARdata &par_info3, int row)
{
	vector< AD<double> > resVec( dat3.nCats, 0);
	vector< AD<double> > lp( dat3.nCats, 0), explp( dat3.nCats, 0);
	AD<double> sumexplp=0;
		
	for( size_t ii=0; ii<dat3.nCats-1; ii++){
		for( size_t jj=0; jj<par_info3.piIndex.size(); jj++)
			lp.at( ii) += dat3.X.at( row).at( par_info3.piIndex.at( jj)) * params3.at( par_info3.nPhi + ii*par_info3.nPi + jj);
		explp.at( ii) = exp( lp.at( ii));
		sumexplp += explp.at( ii);
	}
	
	resVec.at( 0) = 1;
	for( size_t ii=0; ii<dat3.nCats-1; ii++){
		resVec.at( ii+1) = explp.at( ii) / ( 1+sumexplp);
		resVec.at( 0) -= resVec.at( ii+1);
	}
	
	return( resVec);
}	

AD<double> calcSinglePhi( const vector< AD<double> > &params3, const RMCdata &dat3, const PARdata &par_info3, int row, int fromState)
{
	AD<double> lp=0, phi=0;
	int temp=0;
	
	for( size_t i=0; i<(fromState-1); i++)
		temp += par_info3.phiIndex.at( i).size();
	for( size_t jj=0; jj<par_info3.phiIndex.at( fromState-1).size(); jj++)
		lp += dat3.X.at( row).at( par_info3.phiIndex.at( fromState-1).at( jj)) * params3.at( temp + jj);
	phi = 1 - 1/(1+exp( lp));
	
	return phi;	
}
		
		

//RMCdata functions

RMCdata::RMCdata(){}
RMCdata::~RMCdata(){}

void RMCdata::setVals( int R_nCats1, int R_nCovars1, int R_nTot1, int R_nChains1, double *R_X1, int *R_nObs1, int *R_states1)
{
	nCats = R_nCats1;
	nCovars = R_nCovars1;
	nTot = R_nTot1;
	nChains = R_nChains1;
			
	for( size_t i=0; i<nTot; i++){
		X.push_back( vector<double>());
		for( size_t j=0; j<nCovars; j++)
			X.at(i).push_back( R_X1[j*nTot+i]);
		states.push_back( R_states1[i]);
	}
	
	for( size_t i=0; i<nChains; i++)
		nObs.push_back( R_nObs1[i]);
}


//PARdata functions

PARdata::PARdata(){}
PARdata::~PARdata(){}

void PARdata::setVals( int *R_phiID1, int *R_piID1, int R_nCovars1, int R_nCats1)
{
	nPi = 0;
	nPhi = 0;
	
	for( size_t i=0; i<R_nCovars1; i++){
		if( R_piID1[i] == 1){
			nPi++;
			piIndex.push_back( i);
		}
	}
	for( int i=0; i<R_nCats1; i++){
		phiIndex.push_back( vector<int>());
		for( int j=0; j<R_nCovars1; j++){
			if( R_phiID1[i*R_nCovars1+j] == 1){
				nPhi++;
				phiIndex.at( i).push_back( j);
			}
		}
	}
	
	nPars = nPhi + nPi*(R_nCats1-1);
}




