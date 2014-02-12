
// EM step of GPA

#include <Rcpp.h>

RcppExport SEXP cppEM( 
	SEXP gwasPval, SEXP annMat, SEXP useAnn, SEXP pleiotropyH0, SEXP empiricalNull, 
	SEXP maxIter, SEXP stopping, SEXP epsStopLL,
	SEXP binaryMat, SEXP betaAlpha, SEXP betaAlphaNull, SEXP pis, SEXP q1, 
	SEXP lbPi, SEXP lbBetaAlpha, SEXP lbQ, SEXP vDigit, SEXP verbose ) { 

	using namespace Rcpp;
	    
	//////////////////////////////////////////////////////
	//
    // initialize variables
    //
    //////////////////////////////////////////////////////
	
	int iter = as<int>( maxIter );
	int empiricalNullVal = as<int>( empiricalNull );
	int useAnnVal = as<int>( useAnn );
	int verboseVal = as<int>( verbose );
	int pleiotropyH0Val = as<int>( pleiotropyH0 );
	
	int stoppingRule = as<int>( stopping );
	double epsStopLLVal = as<double>( epsStopLL );
	//double epsStopPEVal = as<double>( epsStopPE );
	double vDigitVal = as<double>( vDigit );
	
	double lbPiVal = as<double>( lbPi );
	double lbBetaAlphaVal = as<double>( lbBetaAlpha );
	double lbQVal = as<double>( lbQ );
	
	Rcpp::NumericMatrix binaryMatMat( binaryMat );
	Rcpp::NumericMatrix gwasPvalMat( gwasPval );
	Rcpp::NumericMatrix annMatMat( annMat );   
	
	Rcpp::NumericVector betaAlphaVec( Rcpp::clone(betaAlpha) );
	Rcpp::NumericVector betaAlphaNullVec( Rcpp::clone(betaAlphaNull) );
	Rcpp::NumericVector pisVec( Rcpp::clone(pis) );	
	Rcpp::NumericMatrix q1Mat( Rcpp::clone(q1) );
	
	int nBin = gwasPvalMat.nrow();
	int nGWAS = gwasPvalMat.ncol();
	int nComp = binaryMatMat.nrow();
	int nAnn = 1;
	if ( useAnnVal == 1 ) {
		nAnn = annMatMat.ncol();
	}
	for ( int j = 0; j < nGWAS; j++ ) {
		betaAlphaNullVec[j] = 1;
	}
	
	Rcpp::NumericMatrix logPvalMat( nBin, nGWAS );
	for ( int j = 0; j < nBin; j++ ) {
    	for ( int k = 0; k < nGWAS; k++ ) {  		
	    	logPvalMat( j, k ) = log( gwasPvalMat( j, k ) );
    	}
	}
	    
	//////////////////////////////////////////////////////
	//
    // EM algorithm
    //
    //////////////////////////////////////////////////////
	
	Rcpp::NumericMatrix zMat( nBin, nComp );
	
	Rcpp::NumericMatrix piMat( iter, nComp );
	Rcpp::NumericMatrix betaAlphaMat( iter, nGWAS );
	Rcpp::NumericMatrix betaAlphaNullMat( iter, nGWAS );
	Rcpp::NumericVector loglik( iter );	
	Rcpp::NumericMatrix q1Array( iter, nAnn * nComp );
	
	for ( int g = 0; g < nComp; g++ ) {  
		piMat( 0, g ) = pisVec[g];
	}
	for ( int k = 0; k < nGWAS; k++ ) {
		betaAlphaMat( 0, k ) = betaAlphaVec[k];
		betaAlphaNullMat( 0, k ) = betaAlphaNullVec[k];
	}
	loglik[0] = -1.0 / 0.0;
	if ( useAnnVal == 1 ) {
		for ( int d = 0; d < nAnn; d++ ) {
			for ( int g = 0; g < nComp; g++ ) {  
				q1Array( 0, nComp * d + g ) = q1Mat( d, g );
			}
		}
	}
	//double newDiff = 0;
	//double maxDiff = 0;
	//double relativeDiffNumerator = 0;
	//double relativeDiffDenominator = 0;
	double relativeDiff = 0;
	double AitkenValNext = 0;
	double AitkenValCurrent = 0;	
	double AitkenNext = 0;
	double AitkenCurrent = 0;	
    double AitkenDiff = 0;
	
	Rcpp::NumericMatrix zMargMat( nBin, nGWAS );
	Rcpp::NumericVector sumZvec( nComp );
	Rcpp::NumericVector sumZmargVec( nGWAS );
	
	Rcpp::NumericMatrix betaDist( nBin, nGWAS );
	Rcpp::NumericMatrix betaNullDist( nBin, nGWAS );
	
	Rcpp::NumericVector piMargVec( nGWAS );
	
	double denom = 0;	
	double loglikRow = 0;
	double loglikElement = 0;
	
	double sumPi = 0;
	
	int ii;
    
    if ( verboseVal >= 3 ) {
        Rcout << "------------ iteration:" << 1 << "------------" << std::endl;
        
		Rcout << "GWAS combination: ";
		for ( int g = 0; g < nComp; g++ ) {  
			for ( int k = 0; k < nGWAS; k++ ) {
				Rcout << binaryMatMat( g, k );
			}
			Rcout << " ";
		}
		Rcout << std::endl;
		
		Rcout << "pi: ";
		for ( int g = 0; g < nComp; g++ ) {  
			Rcout << round( vDigitVal * pisVec[g] ) / vDigitVal << " ";
		}
		Rcout << std::endl;
		
		Rcout << "alpha: ";
		for ( int k = 0; k < nGWAS; k++ ) {
			Rcout << round( vDigitVal * betaAlphaVec[k] ) / vDigitVal << " ";
		}
		Rcout << std::endl;
		
		if ( empiricalNullVal == 1 ) { 
			Rcout << "alpha0: ";
			for ( int k = 0; k < nGWAS; k++ ) {
				Rcout << round( vDigitVal * betaAlphaNullVec[k] ) / vDigitVal << " ";
			}
			Rcout << std::endl;
		}
		
		if ( useAnnVal == 1 ) {
			Rcout << "q: " << std::endl;
			for ( int d = 0; d < nAnn; d++ ) {
				for ( int g = 0; g < nComp; g++ ) {  
					Rcout << round( vDigitVal * q1Mat( d, g ) ) / vDigitVal << " ";
				}
				Rcout << std::endl;
			}
		}		
    }
	
	// pre-calculate Beta distribution
	
	for ( int k = 0; k < nGWAS; k++ ) {
		for ( int j = 0; j < nBin; j++ ) {
			betaDist( j, k ) = betaAlphaVec[k] * pow( gwasPvalMat( j, k ), ( betaAlphaVec[k] - 1 ) );
		}
	}
	
	if ( empiricalNullVal == 1 ) { 
		for ( int k = 0; k < nGWAS; k++ ) {
			for ( int j = 0; j < nBin; j++ ) {
				betaNullDist( j, k ) = betaAlphaNullVec[k] * pow( gwasPvalMat( j, k ), ( betaAlphaNullVec[k] - 1 ) );
			}
		}
	}
    
	for ( ii = 1; ii < iter; ii++ ) {
        if ( verboseVal >= 3 ) {
            Rcout << "------------ iteration:" << (ii+1) << "------------" << std::endl;
        }
	    
    	//////////////////////////////////////////////////////
    	//
	    // E step
	    //
	    //////////////////////////////////////////////////////
		
		for ( int j = 0; j < nBin; j++ ) {
			for ( int g = 0; g < nComp; g++ ) {  	    
				    
		        // mixing proportion
		        
		        zMat( j, g ) = pisVec[g];
		        
		        // emission for GWAS
		    
			    for ( int k = 0; k < nGWAS; k++ ) {
			        if ( binaryMatMat( g, k ) == 1 ) {
			        	// signal SNP
			        	
			        	zMat( j, g ) *= betaDist( j, k );
		        	} else if ( empiricalNullVal == 1 ) {
			        	// null SNP, if null is empirically estimated
			        	
			        	zMat( j, g ) *= betaNullDist( j, k );
		        	}
		        }
			        	
		    	// if annotation data exists, incorporate q0 & q1
		    	
		    	if ( useAnnVal == 1 ) {
		        	for ( int d = 0; d < nAnn; d++ ) {
			        	if ( annMatMat( j, d ) == 1 ) {
			        		zMat( j, g ) *= q1Mat( d, g );
		        		} else {
			        		zMat( j, g ) *= ( 1 - q1Mat( d, g ) );
		        		}
		        	}
	        	}
	    	}
	    	
	    	// normalize Z
	    	
	    	denom = 0;
	    	for ( int g = 0; g < nComp; g++ ) { 
		    	denom += zMat( j, g );
	    	}
	    	for ( int g = 0; g < nComp; g++ ) { 
	    		zMat( j, g ) = zMat( j, g ) / denom;
			}
	    }
	    
    	// marginal Z
    	
		for ( int j = 0; j < nBin; j++ ) {
    		for ( int k = 0; k < nGWAS; k++ ) {  				
	    		zMargMat( j, k ) = 0;
	    		for ( int g = 0; g < nComp; g++ ) {  
	    			if ( binaryMatMat( g, k ) == 1 ) {
		    			zMargMat( j, k ) += zMat( j, g );
	    			}
    			}
			}
    	}
	    
    	//////////////////////////////////////////////////////
    	//
	    // M step
	    //
	    //////////////////////////////////////////////////////	
	    
	    // pre-computation
	    
	    for ( int g = 0; g < nComp; g++ ) {  
        	sumZvec[g] = 0;
		    for ( int j = 0; j < nBin; j++ ) {
			    sumZvec[g] += zMat( j, g );
		    }
	    }
	    
	    for ( int k = 0; k < nGWAS; k++ ) {  
        	sumZmargVec[k] = 0;        	
		    for ( int j = 0; j < nBin; j++ ) {
			    sumZmargVec[k] += zMargMat( j, k );
		    }
	    }
            
        // M step: update pi1
        
        if ( pleiotropyH0Val == 0 ) {
	        // EM algorithm for original GPA model
	        
	        for ( int g = 0; g < nComp; g++ ) {  
	        	pisVec[g] = sumZvec[g] / nBin;
	    	}
    	} else {        
	        // EM algorithm for GPA model under H0 of pleitropy test	
	        
	        for ( int k = 0; k < nGWAS; k++ ) {
		        piMargVec[k] = sumZmargVec[k] / nBin;
	        }
	        
			for ( int g = 0; g < nComp; g++ ) {  
				pisVec[g] = 1;
				for ( int k = 0; k < nGWAS; k++ ) {
					if ( binaryMatMat( g, k ) == 0 ) {
						pisVec[g] *= ( 1 - piMargVec[k] );
					} else {
						pisVec[g] *= piMargVec[k];
					}
				}
			}
		}
    
        // M step: update Beta parameters (signal)
		
        for ( int k = 0; k < nGWAS; k++ ) {
	        denom = 0;
	        for ( int j = 0; j < nBin; j++ ) {	        	
		        denom += zMargMat( j, k ) * ( -logPvalMat( j, k ) );
        	}
        	betaAlphaVec[k] = sumZmargVec[k] / denom;
        	
        	// check the range of alpha: 0 < alpha <= 1
        	
        	if ( betaAlphaVec[k] < lbBetaAlphaVal ) {
	        	betaAlphaVec[k] = lbBetaAlphaVal;
        	}
        	if ( betaAlphaVec[k] > 1 - lbBetaAlphaVal ) {
	        	betaAlphaVec[k] = 1 - lbBetaAlphaVal;
        	}
        }
        
        // M step: update Beta parameters (null)
		
        if ( empiricalNullVal == 1 ) {
	        for ( int k = 0; k < nGWAS; k++ ) {
		        denom = 0;
		        for ( int j = 0; j < nBin; j++ ) {	        	
			        denom += ( 1 - zMargMat( j, k ) ) * ( -logPvalMat( j, k ) );
	        	}
	        	betaAlphaNullVec[k] = ( nBin - sumZmargVec[k] ) / denom;
        	
	        	// check the range of alpha: alpha <= alpha_null <= 1
	        	
	        	if ( betaAlphaNullVec[k] < betaAlphaVec[k] ) {
		        	betaAlphaNullVec[k] = betaAlphaVec[k];
	        	}
	        	if ( betaAlphaNullVec[k] > 1 - lbBetaAlphaVal ) {
		        	betaAlphaNullVec[k] = 1 - lbBetaAlphaVal;
	        	}
	        }
        }
		
		// M step: if annotation data exists, update q1
        
		if ( useAnnVal == 1 ) {	        
	        for ( int g = 0; g < nComp; g++ ) {		        
		        for ( int d = 0; d < nAnn; d++ ) {
			        q1Mat( d, g ) = 0;
			        for ( int j = 0; j < nBin; j++ ) {
			        	q1Mat( d, g ) += zMat( j, g ) * annMatMat( j, d );
		        	}
		        	q1Mat( d, g ) /= sumZvec[g];
		        }
	        }
        }
	
		// update Beta distribution
		
		for ( int k = 0; k < nGWAS; k++ ) {
			for ( int j = 0; j < nBin; j++ ) {
				betaDist( j, k ) = betaAlphaVec[k] * pow( gwasPvalMat( j, k ), ( betaAlphaVec[k] - 1 ) );
			}
		}
		
		if ( empiricalNullVal == 1 ) { 
			for ( int k = 0; k < nGWAS; k++ ) {
				for ( int j = 0; j < nBin; j++ ) {
					betaNullDist( j, k ) = betaAlphaNullVec[k] * pow( gwasPvalMat( j, k ), ( betaAlphaNullVec[k] - 1 ) );
				}
			}
		}
	    
    	//////////////////////////////////////////////////////
    	//
	    // safe guard for M step
	    //
	    //////////////////////////////////////////////////////	
        
        // check the lower bound for pi: pi > 0
	
        sumPi = 0;
        for ( int g = 0; g < nComp; g++ ) {	
			if ( pisVec[g] < lbPiVal ) {
				pisVec[g] = lbPiVal;
			}
			sumPi += pisVec[g];
		}
		for ( int g = 0; g < nComp; g++ ) {	
			pisVec[g] = pisVec[g] / sumPi;
		}
		
		// check the lower bound for q1: 0 < q1 < 1
		
		if ( useAnnVal == 1 ) {
	        for ( int g = 0; g < nComp; g++ ) {		        
		        for ( int d = 0; d < nAnn; d++ ) {		
					if ( q1Mat( d, g ) < lbQVal ) {	
						q1Mat( d, g ) = lbQVal;
					}
					
					if ( q1Mat( d, g ) > 1 - lbQVal ) {
						q1Mat( d, g ) = 1 - lbQVal;
					}
				}
			}
		}
	    
    	//////////////////////////////////////////////////////
    	//
	    // track parameter estimates & log likelihood
	    //
	    //////////////////////////////////////////////////////	
        
        // track parameter estimates
	
		for ( int g = 0; g < nComp; g++ ) {  
			piMat( ii, g ) = pisVec[g];
		}
		for ( int k = 0; k < nGWAS; k++ ) {
			betaAlphaMat( ii, k ) = betaAlphaVec[k];
			betaAlphaNullMat( ii, k ) = betaAlphaNullVec[k];
		}
		if ( useAnnVal == 1 ) {
			for ( int d = 0; d < nAnn; d++ ) {
				for ( int g = 0; g < nComp; g++ ) {  
					q1Array( ii, nComp * d + g ) = q1Mat( d, g );
				}
			}
		}
        
        if ( verboseVal >= 3 ) {
			Rcout << "GWAS combination: ";
			for ( int g = 0; g < nComp; g++ ) {  
				for ( int k = 0; k < nGWAS; k++ ) {
					Rcout << binaryMatMat( g, k );
				}
				Rcout << " ";
			}
			Rcout << std::endl;
			
			Rcout << "pi: ";
			for ( int g = 0; g < nComp; g++ ) {  
				Rcout << round( vDigitVal * pisVec[g] ) / vDigitVal << " ";
			}
			Rcout << std::endl;
			
			Rcout << "alpha: ";
			for ( int k = 0; k < nGWAS; k++ ) {
				Rcout << round( vDigitVal * betaAlphaVec[k] ) / vDigitVal << " ";
			}
			Rcout << std::endl;
			
			if ( empiricalNullVal == 1 ) { 
				Rcout << "alpha0: ";
				for ( int k = 0; k < nGWAS; k++ ) {
					Rcout << round( vDigitVal * betaAlphaNullVec[k] ) / vDigitVal << " ";
				}
				Rcout << std::endl;
			}
			
			if ( useAnnVal == 1 ) {
				Rcout << "q: " << std::endl;
				for ( int d = 0; d < nAnn; d++ ) {
					for ( int g = 0; g < nComp; g++ ) {  
						Rcout << round( vDigitVal * q1Mat( d, g ) ) / vDigitVal << " ";
					}
					Rcout << std::endl;
				}
			}				
        }
        
        // track log likelihood
		
		loglik[ii] = 0;
	    
		for ( int j = 0; j < nBin; j++ ) {
	    	loglikRow = 0;
	    	
			for ( int g = 0; g < nComp; g++ ) {  	    
				    
		        // mixing proportion
		        
		        loglikElement = pisVec[g];
		        
		        // emission for GWAS
		    
			    for ( int k = 0; k < nGWAS; k++ ) {
			        if ( binaryMatMat( g, k ) == 1 ) {
			        	// signal SNP
			        	
			        	loglikElement *= betaDist( j, k );
		        	} else if ( empiricalNullVal == 1 ) {
			        	// null SNP, if null is empirically estimated
			        	
			        	loglikElement *= betaNullDist( j, k );
		        	}
		        }
			        	
		    	// if annotation data exists, incorporate q0 & q1
		    	
		    	if ( useAnnVal == 1 ) {
		        	for ( int d = 0; d < nAnn; d++ ) {
			        	if ( annMatMat( j, d ) == 1 ) {
			        		loglikElement *= q1Mat( d, g );
		        		} else {
			        		loglikElement *= ( 1 - q1Mat( d, g ) );
		        		}
		        	}
	        	}
	        	
	        	loglikRow += loglikElement;
	    	}
	    	
	    	// update log likelihood
	    	
		    loglik[ii] += log( loglikRow );
	    }	
    	
        if ( verboseVal >= 3 ) {
	        Rcout << "Increment in log likelihood: ";
	        Rcout << ( loglik[ii] - loglik[(ii-1)] ) << std::endl;
        }
	    
    	//////////////////////////////////////////////////////
    	//
	    // stopping rule
	    //
	    //////////////////////////////////////////////////////	
        
        // stop iterations if converges
        // - Note #1. Not use safe guard against decreasing log likelihood.
        // - Note #2. Stop based only on log likelihood (not based on parameter estimates)
        
       	if ( stoppingRule == 1 && ii > 10 ) {
	       	// based on absolute difference
	       	
	       	if ( ( loglik[ii] - loglik[(ii-1)] ) < epsStopLLVal ) {
				Rcout << "Info: EM algorithm stops in " << ii << "-th iteration because there is no improvements in log likelihood." << std::endl;
					
				break;
			}
       	} else if ( stoppingRule == 2 && ii > 10 ) {
	       	// based on relative difference
	       	
	       	//relativeDiff = ( loglik[ii] - loglik[(ii-1)] ) / loglik[(ii-1)];
			relativeDiff = ( loglik[ii] - loglik[(ii-1)] ) / loglik[ii];
	       	if ( relativeDiff < 0 ) {
		       	relativeDiff = - relativeDiff;
	       	}
	       	//relativeDiffNumerator = ( loglik[ii] - loglik[(ii-1)] );
	       	//relativeDiffDenominator = loglik[(ii-1)];
	       	//if ( relativeDiffDenominator < 0 ) {
		    //   	relativeDiffDenominator = - relativeDiffDenominator;
	       	//}
	       	//relativeDiff = relativeDiffNumerator / relativeDiffDenominator;
	       	
	       	if ( relativeDiff < epsStopLLVal ) {
				Rcout << "Info: EM algorithm stops in " << ii << "-th iteration because there is no improvements in log likelihood." << std::endl;
					
				break;
			}
    	} else if ( stoppingRule == 3 && ii > 10 ) {
			// based on Aitken stopping rule
			
			AitkenValNext = ( loglik[ii] - loglik[(ii-1)] ) / ( loglik[(ii-1)] - loglik[(ii-2)] );
			AitkenValCurrent = ( loglik[(ii-1)] - loglik[(ii-2)] ) / ( loglik[(ii-2)] - loglik[(ii-3)] );
			
			AitkenNext = loglik[(ii-1)] + ( loglik[ii] - loglik[(ii-1)] ) / ( 1 - AitkenValNext );
			AitkenCurrent = loglik[(ii-2)] + ( loglik[(ii-1)] - loglik[(ii-2)] ) / ( 1 - AitkenValCurrent );
			
		    AitkenDiff = AitkenNext - AitkenCurrent;
		    if ( AitkenDiff < 0 ) {
			    AitkenDiff = - AitkenDiff;
		    }
		    if ( AitkenDiff < epsStopLLVal ) {
		        
				Rcout << "Info: EM algorithm stops in " << ii << "-th iteration because there is no improvements in log likelihood." << std::endl;
					
				break;
			}
		} else if ( ii > 10 ) {
			Rcout << "Info: Inappropriate stopping rule specification. Will be based on relative difference in log likelihood." << std::endl;
			stoppingRule = 2;
		}
			
    }
    
	// report if EM algorithm hits the max # iterations
	
	if ( ii == iter ) {
		Rcout << "Info: EM algorithm stops because it reaches the specified maximum number of iterations (maxIter = " << iter << ")." << std::endl;
	}
    
    // prepare the list to return
	
	Rcpp::List emResult;
	
	emResult["Z"] = zMat;
	emResult["Zmarg"] = zMargMat;	
	
	emResult["pis"] = pisVec;
	emResult["betaAlpha"] = betaAlphaVec;
	emResult["betaAlphaNull"] = betaAlphaNullVec;
	emResult["q1"] = q1Mat;
	
	emResult["piMat"] = piMat;
	emResult["betaAlphaMat"] = betaAlphaMat;
	emResult["betaAlphaNullMat"] = betaAlphaNullMat;
	emResult["loglik"] = loglik;
	emResult["q1Array"] = q1Array;
	
	emResult["maxIter"] = ii;
	
	return( emResult );	
}
