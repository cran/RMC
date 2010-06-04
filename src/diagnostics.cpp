#include"RMC.h"

//###################################
//####		Function Headers		####
//###################################

extern "C" void diagnostics( int *, int *, double *, int *, int *, int *, int *, double *, int *, int *, int *, int *, double *, int *, double *, int *, int *);
void calcJumpDist( const RMCdata &, const PARdata &, const vector< AD<double> > &, double *, int *);
int findFirst( int, const RMCdata &, int, int, vector< vector<int> > &);
int findLast( int, const RMCdata &, int, int, vector< vector<int> > &);
double calcTel( const vector< AD<double> > &, const RMCdata &, const PARdata &, const int, int first);
void summURV( const vector< vector<double> > &, const int &, int *, double *);
void summJumps( const vector< vector<int> > &, const int &, int *);
double urvFirstPatch( const vector< AD<double> > &, const RMCdata &, const PARdata &, int, int);
double urvLastPatch( const vector< AD<double> > &, const RMCdata &, const PARdata &, int, int);

//###################################
//####	Function Definitions		####
//###################################

extern "C" {
	void diagnostics( int *Rstates, int *Rnobs, double *RX, int *nCats, int *nCovars, int *nTot, int *nChains, double *Rpars, int *RphiID, int *RpiID, int *includeEnds, int *jumps, double *patch, int *lens, double *RjumpDists, int *RdistCuts, int *doDist)
{
	RMCdata dat;
	dat.setVals( *nCats, *nCovars, *nTot, *nChains, RX, Rnobs, Rstates);
	PARdata parINF;
	parINF.setVals( RphiID, RpiID, *nCovars, *nCats);
	
	vector< AD<double> > params( parINF.nPars, 0);
	for( size_t i=0; i<parINF.nPars; i++)
		params.at( i) = Rpars[i];
	
	vector< vector<double> > urv( dat.nCats, vector<double>());
	vector< vector<int> > cross( dat.nCats, vector<int> ( dat.nCats, 0));	
	
	size_t inPatch = 1, notNA = 1;
	int start, stop = 0;
	int patchStart, patchStop;
	int j = 0;
	
	double tempURV = 0.;
	double curpdf=0., prevpdf=0.;
	
	double dTP, prevTP;
	double tmp;
	
	for( size_t i=0; i<dat.nChains; i++){
		start = stop;
		stop += dat.nObs.at( i);
		patchStart = findFirst( *includeEnds, dat, start, stop, cross);			//find first patch
		patchStop = findLast( *includeEnds, dat, patchStart, stop, cross);		//find last patch
		if( patchStart < stop){
			tmp = urvFirstPatch( params, dat, parINF, start, patchStart);
			if( tmp != -9)
				urv.at( dat.states.at( start)-1).push_back( tmp);
			tmp = urvLastPatch( params, dat, parINF, patchStop, stop);
			if( tmp != -9)
				urv.at( dat.states.at( stop-1)-1).push_back( tmp);
		}
		while( ( patchStart < patchStop)){// & ( patchStart < stop)){
			inPatch = 1;
			notNA = 1;
			dTP = calcTel( params, dat, parINF, patchStart, 1);
			prevTP = dTP;
			curpdf = (1-dTP);
			tempURV = (1-dTP)/2;
			j = 1;
			while( inPatch==1){
				if( ( dat.states.at( patchStart+j) != -9) && ( dat.states.at( patchStart+j-1) != -9)){
					if( dat.states.at( patchStart+j) == dat.states.at( patchStart+j-1)){
						inPatch = 1;
						prevTP = dTP;
						dTP = calcTel( params, dat, parINF, patchStart+j, 0);
						prevpdf = curpdf;
						curpdf = prevpdf * prevTP * ( 1-dTP) / ( 1-prevTP);
						tempURV = tempURV + 0.5 * prevpdf + 0.5 * curpdf;
						cross.at( dat.states.at( patchStart+j-1)-1).at( dat.states.at( patchStart+j)-1)++;
						j++;
					}
					else{
						inPatch = 0;
						urv.at( dat.states.at( patchStart)-1).push_back( tempURV);
					}		
					notNA = 1;
				}
				else{
					notNA = 0;
					patchStart = findFirst( 0, dat, patchStart+j, patchStop, cross);
					inPatch = 0;
				}
			}
			if( notNA == 1){
				patchStart += j;
				cross.at( dat.states.at( patchStart-1)-1).at( dat.states.at( patchStart)-1)++;
			}		
		}
	}	
	summURV( urv, dat.nCats, lens, patch);
	summJumps( cross, dat.nCats, jumps);
	if( *doDist==1)
		calcJumpDist( dat, parINF, params, RjumpDists, RdistCuts);	
}
}



void summJumps( const vector< vector<int> > &cross1, const int &n_cats, int *jumps1){
	int place = 0;
	for( size_t i=0; i<n_cats; i++)
		for( size_t j=0; j<n_cats; j++){
			jumps1[place] = cross1.at(i).at(j);
			place++;
		}
}

void summURV( const vector< vector<double> > & urv1, const int &n_cats, int *lens1, double *patch1){
	int place = 0;
	
	for( size_t i=0; i<n_cats; i++){
		lens1[i] = urv1.at(i).size();
		for( size_t j=0; j<urv1.at(i).size(); j++){
			patch1[place] = urv1.at(i).at(j);
			place++;
			}
	}
}

double calcTel( const vector< AD<double> > &params1, const RMCdata &dat1, const PARdata &parINF1, const int locat, int first){
	vector< AD<double> > piVec( dat1.nCats, 0);
	AD<double> phi=0.;
	double el=0.;
	
	if( first!=1){
		piVec = calcPiVec( params1, dat1, parINF1, locat);		
		phi = calcSinglePhi( params1, dat1, parINF1, locat, dat1.states.at( locat-1));	
		el = 1 - Value( phi ) + Value( phi) * Value( piVec.at( dat1.states.at( locat-1) - 1 ) );
	}
	else{
		piVec = calcPiVec( params1, dat1, parINF1, locat);
		phi = calcSinglePhi( params1, dat1, parINF1, locat, dat1.states.at( locat-1+1));
		el = 1 - Value( phi ) + Value( phi) * Value( piVec.at( dat1.states.at( locat-1+1) - 1 ) );
	}

	
	return( el);
}

double urvFirstPatch( const vector< AD<double> > &params1, const RMCdata &dat1, const PARdata &parINF1, int firstStart1, int firstStop1){

	double dTP, prevTP, pdf;
	double culProb=0;
	double urv1;

	dTP = calcTel( params1, dat1, parINF1, firstStop1-1, 1);
	prevTP = dTP;
	pdf = 1- dTP;
	culProb += pdf;
	for( size_t i=( firstStop1-1); i>firstStart1; i--){
		if( dat1.states.at( i-1) == -9)
			return( -9);
		dTP = calcTel( params1, dat1, parINF1, i-1, 1);
		pdf = pdf * prevTP * ( 1-dTP) / ( 1-prevTP);
		culProb += pdf;
		prevTP = dTP;
	}
	if( dat1.states.at( firstStart1) == -9)
		return( -9);
	urv1 = culProb + ( 1-culProb) / 2;
	
	return( urv1);
}

double urvLastPatch( const vector< AD<double> > &params1, const RMCdata &dat1, const PARdata &parINF1, int lastStart1, int lastStop1){

	double dTP, prevTP, pdf;
	double culProb=0;
	double urv1;

	if( dat1.states.at( lastStart1+1)==-9)
		return( -9);
	dTP = calcTel( params1, dat1, parINF1, lastStart1+1, 1);
	prevTP = dTP;
	pdf = 1 - dTP;
	culProb += pdf;
	for( size_t i=(lastStart1+1); i<(lastStop1-1); i++){
		if( dat1.states.at( i) == -9)
			return( -9);
		dTP = calcTel( params1, dat1, parINF1, i, 1);
		pdf = pdf * prevTP * (1-dTP) / (1-prevTP);
		culProb += pdf;
		prevTP = dTP;
	}
	if( dat1.states.at( lastStop1-1) == -9)
		return( -9);
	urv1 = culProb + ( 1-culProb) / 2;

	return( urv1);
}


int findFirst( int RincludeEnds, const RMCdata &dat1, int start1, int stop1, vector< vector<int> > &cross1){
	size_t flag=0;
	int newst = start1;
	
	if( RincludeEnds!=1)
		while( flag==0){
			if( newst < stop1 - 1){
				if( ( dat1.states.at( newst) != -9) && ( dat1.states.at( newst+1) != -9)){
					cross1.at( dat1.states.at( newst)-1).at( dat1.states.at( newst+1)-1)++;
					if( dat1.states.at( newst) != dat1.states.at(newst+1))
						flag=1;
				}
				newst++;
			}
			else{
				newst = stop1;
				flag=1;
			}
		}
	return( newst);
}
				
int findLast( int RincludeEnds, const RMCdata &dat1, int start1, int stop1, vector< vector<int> > &cross1){
	size_t flag=0;
	int newend = stop1-1;
	
	if( RincludeEnds!=1){
		while( flag==0){
			if( newend >= start1+1){		//Changed +1
				if( ( dat1.states.at( newend) != -9) && ( dat1.states.at( newend-1) != -9)){
					cross1.at( dat1.states.at( newend-1)-1).at( dat1.states.at( newend)-1)++;
					if( dat1.states.at( newend) != dat1.states.at( newend-1)){
						flag = 1;
						cross1.at( dat1.states.at( newend-1)-1).at( dat1.states.at( newend)-1)--;	//fixing up double counting (occurs later)
					}
				}
				newend--;
			}
			else{
				newend = start1;
				flag = 1;
			}
		}
	}
	return( newend);
}
					
void calcJumpDist( const RMCdata & dat1, const PARdata & parINF1, const vector< AD<double> > & params1, double *RjumpDists, int *RdistCuts)
{
	vector< vector< vector<double> > > jumpDist( dat1.nCats, vector< vector<double> >(dat1.nCats, vector<double>(1, -9.)));
	vector<double> pD;
		
	int start = 0;
	int stop = dat1.nObs.at( 0);
	vector< AD<double> > piVec( dat1.nCats, 0);
	AD<double> phi=0.;
	double dTP=0., sumPi=0.;
	int curState=0, siz=0;
	int place=1;
	
	for( size_t i=0; i<dat1.nChains; i++){
		for( size_t j=start+1; j<stop; j++){
			if( ( dat1.states.at( j-1) != -9) && ( dat1.states.at( j) != -9) && ( dat1.states.at( j-1) != dat1.states.at( j))){
				curState = dat1.states.at( j-1);
				piVec = calcPiVec( params1, dat1, parINF1, j);
				sumPi = 0;
				for( size_t k=0; k<dat1.nCats; k++)
					if( k!=( curState-1))
						sumPi += Value( piVec.at( k));
				for( size_t k=0; k<dat1.nCats; k++){
					if( k!=( curState-1)){
						dTP = Value( piVec.at( k)) / sumPi;
						if( jumpDist.at( curState-1).at( k).size()==1){
							jumpDist.at( curState-1).at( k).push_back( 1-dTP);
							jumpDist.at( curState-1).at( k).push_back( dTP);
						}
						else{
							siz = jumpDist.at( curState-1).at(k).size();
							pD = jumpDist.at( curState-1).at(k);
							jumpDist.at(curState-1).at(k).at(1) = pD.at(1)*(1-dTP);
							jumpDist.at( curState-1).at(k).push_back( pD.at( siz-1) * dTP);
							for( size_t l=2; l<siz; l++)
								jumpDist.at(curState-1).at(k).at(l) = pD.at( l-1)*dTP + pD.at(l)*(1-dTP);
						}
					}
				}
			}
		}
		if( i+1 < dat1.nChains){
			start = stop + 1;
			stop += dat1.nObs.at( i+1);
		}
	}	
	place = 0;
	for( size_t i=0; i<dat1.nCats; i++)
		for( size_t j=0; j<dat1.nCats; j++){
			RdistCuts[i*dat1.nCats+j] = (int)jumpDist.at( i).at( j).size();
			for( size_t k=0; k<RdistCuts[i*dat1.nCats+j]; k++){
				RjumpDists[place] = jumpDist.at( i).at( j).at( k);
				place++;
			}
		}
}



























