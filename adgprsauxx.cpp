/**
  \functions to create adgprs input files (implementation) for reservoir optimization research 
  \ Author : Mehrdad Gharib Shirangi, April 2013, 
  \ Last modified: 2015
*/
#include "adgprsauxx.hpp"
#include <cmath>
#include <fstream>
#include <algorithm>

/* Functions to find max and min element in a vector, starting from index n */
double max_vect(vector<double> vVector, int n)
{
  double dMax = -1e20;
  for(vector<double>::size_type aa = n; aa != vVector.size(); ++aa)
    if (dMax < vVector[aa]) {
      dMax = vVector[aa];
    }
  return dMax;

}
double min_vect(vector<double> vVector, int n)
{
  double dMin = 1e20;
  for(vector<double>::size_type aa = n; aa != vVector.size(); ++aa)
    if (dMin > vVector[aa]) {
      dMin = vVector[aa];
    }
  return dMin;

}

/* Function to return the minimum distance between a set of wells (in 2D reservoir) 
 * On September 01 2013, the function is extended to 3D models by Mehrdad Gharib Shirangi
 * */
double minDistBetweenWells( vector<vector<int> >& AllWells , vector<double> wellBinVars, int Nz)
{
  if (Nz > 1 ) {
    const int nWells = AllWells.size();
    vector<double> dist_vec(nWells*(nWells-1)/2, 0);
    int nCoords = 3;
    double value = 0;
    for ( int i = 0; i < nWells-1; i++ ) {
	if (wellBinVars[i] < - 0.1){  // an injector well 
           for ( int j = i+1; j < nWells; j++ ) {
            	value = 0;
		if (wellBinVars[j] < - 0.1) { // an injector well
            	   for ( int k = 0; k < 2; k++ ) {
                	value += ( AllWells[i][k] - AllWells[j][k] ) * ( AllWells[i][k] - AllWells[j][k] );
            	   }
		} else if (wellBinVars[j] >  0.1) { // a producer well
		   int j1 =  AllWells[j][1];
		   int j2 =  AllWells[j][2];
		   if  (AllWells[i][1] > max(j1, j2) || AllWells[i][1] < min(j1, j2)  ) {
		        double d1 = 0; 
			d1 = ( AllWells[i][1] - j1) * ( AllWells[i][1] - j1 );
		        double d2 = 0; 
			d2 = ( AllWells[i][1] - j2) * ( AllWells[i][1] - j2 );
			value = min(d1,d2);	
		
			value += ( AllWells[i][0] - AllWells[j][0] ) * ( AllWells[i][0] - AllWells[j][0] ); 
		   } else {
		  	value = (AllWells[i][0] - AllWells[j][0] ) * (AllWells[i][0] - AllWells[j][0] );
		   }
                } else {
		     value = 100000000;
		}
            	double ell = i*( nWells-(double)(i+1)/2 )+j-i-1;
            	dist_vec[(int)ell] = sqrt(value);	
//          	cout << "distance between well # " << i+1 << " (inj) and well # " << j+1 << " equals "<< dist_vec[(int)ell] << endl;
           }
      } else if (wellBinVars[i] >  0.1){  // a producer well
           for ( int j = i+1; j < nWells; j++ ) {
                value = 0;
                if (wellBinVars[j] >  0.1) { // a producer  well
                    value = ( AllWells[i][0] - AllWells[j][0] ) * ( AllWells[i][0] - AllWells[j][0] );
                    value += ( AllWells[i][3] - AllWells[j][3] ) * ( AllWells[i][3] - AllWells[j][3] );
                } else if (wellBinVars[j] <  -0.1) { // an injector well
                   int i1 =  AllWells[i][1];
                   int i2 =  AllWells[i][2];
                   if  (AllWells[j][1] > max(i1, i2) || AllWells[j][1] < min(i1, i2)  ) {
                        double d1 = 0;
                        d1 = ( AllWells[j][1] - i1) * ( AllWells[j][1] - i1 );
                        double d2 = 0;
                        d2 = ( AllWells[j][1] - i2) * ( AllWells[j][1] - i2 );
                        value = min(d1,d2);

                        value += ( AllWells[j][0] - AllWells[i][0] ) * ( AllWells[j][0] - AllWells[i][0] );
                   } else {
                        value = (AllWells[i][0] - AllWells[j][0] ) * (AllWells[i][0] - AllWells[j][0] );
                   }
                } else {
                     value = 100000000;
                }

                double ell = i*( nWells-(double)(i+1)/2 )+j-i-1;
                dist_vec[(int)ell] = sqrt(value);
//          	cout << "distance between well # " << i + 1 << " (prod) and well # " << j +1 << " equals "<< dist_vec[(int)ell] << endl;
            }
        } else {
            for ( int j = i+1; j < nWells; j++ ) {
               double ell = i*( nWells-(double)(i+1)/2 )+j-i-1;
               dist_vec[(int)ell] = 10000;
	     }
	}	
    }
    return min_vect(dist_vec,0);
 } else { // 2D
    const int nWells = AllWells.size();
    vector<double> dist_vec(nWells*(nWells-1)/2, 0);
    int nCoords = 3;
    double value = 0;
    for ( int i = 0; i < nWells-1; i++ ) {
        for ( int j = i+1; j < nWells; j++ ) {
            value = 0;
            for ( int k = 0; k < nCoords; k++ ) {
                value += ( AllWells[i][k] - AllWells[j][k] ) * ( AllWells[i][k] - AllWells[j][k] );
            }
            double ell = i*( nWells-(double)(i+1)/2 )+j-i-1;
            dist_vec[(int)ell] = sqrt(value);
        }
    }
    return min_vect(dist_vec,0);
 }
}

/* Function to write the AD GPRS well file */ 
void writeADGPRSwellfile (string wellFileNameExtension, vector<double> wellBinVars, double zeroEps, vector<vector<int> >& wellDefs)
{
    vector<int> prodIndxs;
    vector<int> injIndxs;

    int nWells = 0;
    for (int i = 0; i < wellBinVars.size(); i++) {
	if (wellBinVars[i] > zeroEps) {		// Producers
	    prodIndxs.push_back(i);
	    nWells++;
	} else if (wellBinVars[i] < -zeroEps) {	// Injectors
	    injIndxs.push_back(i);
            nWells++;
	}
    }
    int nprod = prodIndxs.size();
    int ninj  = injIndxs.size();
	string fileName;
	fileName = "simForOptFiles/COMPDAT" + wellFileNameExtension;
	ofstream comp(fileName.c_str());
	comp << "COMPDAT" << endl;
	fileName = "simForOptFiles/WELSPECS" + wellFileNameExtension;
	ofstream specs(fileName.c_str());
	specs << "WELSPECS" << endl;
	fileName = "simForOptFiles/WELLSTRE" + wellFileNameExtension; 
	ofstream stre(fileName.c_str());
	stre << "WELLSTRE" << endl;
/*	stre << "-- composition of injected fluid (gas - oil - water)" << endl; */
        stre << "-- composition of injected fluid (oil - water)" << endl;
    /* Write all producer well information */
    	for (int ii = 0; ii < nprod; ii++) {
      		if (wellBinVars[prodIndxs[ii]] > zeroEps) {
			int ind = prodIndxs[ii];
			int prodConns = max(wellDefs[ind ][1],wellDefs[ind][2]) - min(wellDefs[ind][1],wellDefs[ind][2]) + 1;

		 	for (int jj = 0; jj < prodConns; jj++) { // Horizontal well, increment in y	
			   int k = wellDefs[ prodIndxs[ii] ][3] ;			 
			   if (wellDefs[ind ][1] <= wellDefs[ind ][2]) {
    			      comp << "PRD"<< ii+1 <<"   " << wellDefs[ prodIndxs[ii] ][0]; 
  			      comp << "    " << wellDefs[ prodIndxs[ii] ][1] + jj << "  " << k << "  " << k << " SHUT  *  *  4*  Y/" << endl;
 
		  	   } else {
    			      comp << "PRD"<< ii+1 <<"   " << wellDefs[ prodIndxs[ii] ][0]; 
  			      comp << "    " << wellDefs[ prodIndxs[ii] ][2] + jj << "  " << k << "  " << k << " SHUT  *  *  4*  Y/" << endl; 

			   }
 			}
   			specs << "PRD"<< ii+1 <<"  PRD  " << wellDefs[ prodIndxs[ii] ][0]; 
			specs << "    " << wellDefs[ prodIndxs[ii] ][1]<< " 4* HALT NO /" << endl; 
			//specs << "    " << wellDefs[ prodIndxs[ii] ][1]<< " 4* SHUT NO /" << endl; 
		}
	} 
    	for (int ii = 0; ii < ninj; ii++) {
      		if (wellBinVars[injIndxs[ii]] < -zeroEps) {	
			int ind = injIndxs[ii];
			if (wellDefs[ind ][2] <= wellDefs[ind ][3]) {
    			   comp << "INJ"<< ii+1 <<"   " << wellDefs[ injIndxs[ii] ][0]; 
			   comp << "    " << wellDefs[ ind ][1]<< "  " << wellDefs[ind ][2] << "  "<< wellDefs[ind ][3] << "  SHUT  *  *  4*  Z/" << endl; 		
			} else {
    			   comp << "INJ"<< ii+1 <<"   " << wellDefs[ injIndxs[ii] ][0]; 
			   comp << "    " << wellDefs[ind ][1]<< "  " << wellDefs[ind ][3] << "  "<< wellDefs[ind ][2] << "  SHUT  *  *  4*  Z/" << endl; 		
			}

    			specs << "INJ"<< ii+1 <<"  INJ  " << wellDefs[ injIndxs[ii] ][0]; 
			specs << "    " << wellDefs[ injIndxs[ii] ][1]<< " 4*  * NO /" << endl; 

/*			stre << "INJ"<< ii+1 <<"   0.0 0.0 1.0 /" << endl; */
                        stre << "INJ"<< ii+1 <<"   1.0 0.0 /" << endl;
		}
	}
	comp << "/" << endl;
	comp.close();
	specs << "/" << endl;
        specs.close();
	stre << "/" << endl;
        stre.close();

    return;
}

/*
 * writeCOMPDAT_SPECS_STRE( ext, wellDefs, prodIndxs, injIndxs, Nz,  wellBinVars, zeroEps );
 * Written on 10/23/2014
 */
void writeCOMPDAT_SPECS_STRE( string ext, vector<vector<int> >  wellDefs, vector<int> prodIndxs,  vector<int> injIndxs, int Nz, vector<double> wellBinVars, double zeroEps ) {
    int nprod = prodIndxs.size();
    int ninj  = injIndxs.size();
        string fileName;
        fileName = "simForOptFiles/COMPDAT" + ext;
        ofstream comp(fileName.c_str());
        comp << "COMPDAT" << endl;
        fileName = "simForOptFiles/WELSPECS" + ext;
        ofstream specs(fileName.c_str());
        specs << "WELSPECS" << endl;
        fileName = "simForOptFiles/WELLSTRE" + ext;
        ofstream stre(fileName.c_str());
        stre << "WELLSTRE" << endl;
        stre << "-- composition of injected fluid (water - oil)" << endl;
        /* stre << "-- composition of injected fluid (gas - oil - water)" << endl;
    =========================================================================================== */
    /* Write all producer well information */
	if (Nz == 1) { // 2D
            for (int ii = 0; ii < nprod; ii++) {
                if (wellBinVars[prodIndxs[ii]] > zeroEps) {
                        comp << "PRD"<< ii+1 <<"   " << wellDefs[ prodIndxs[ii] ][0];
                        comp << "    " << wellDefs[ prodIndxs[ii] ][1]<< "      1  1  SHUT  *  *  4*  Z/" << endl;

                        specs << "PRD"<< ii+1 <<"  PRD  " << wellDefs[ prodIndxs[ii] ][0];
                        specs << "    " << wellDefs[ prodIndxs[ii] ][1]<< " 4* HALT NO /" << endl;
                        //specs << "    " << wellDefs[ prodIndxs[ii] ][1]<< " 4* SHUT NO /" << endl;
                }
            }
            for (int ii = 0; ii < ninj; ii++) {
                if (wellBinVars[injIndxs[ii]] < -zeroEps) {
                        comp << "INJ"<< ii+1 <<"   " << wellDefs[ injIndxs[ii] ][0];
                        comp << "    " << wellDefs[ injIndxs[ii] ][1]<< "       1  1  SHUT  *  *  4*  Z/" << endl;

                        specs << "INJ"<< ii+1 <<"  INJ  " << wellDefs[ injIndxs[ii] ][0];
                        specs << "    " << wellDefs[ injIndxs[ii] ][1]<< " *  4*  NO /" << endl;

                        stre << "INJ"<< ii+1 <<"  1.0  0.0 /" << endl;
                        //stre << "INJ"<< ii+1 <<"   0.0 0.0 1.0 /" << endl;
                }
            }
/*     =========================================================================================== */
	} else { // 3D
            for (int ii = 0; ii < nprod; ii++) {
                if (wellBinVars[prodIndxs[ii]] > zeroEps) {
                        int ind = prodIndxs[ii];
                        int prodConns = max(wellDefs[ind ][1],wellDefs[ind][2]) - min(wellDefs[ind][1],wellDefs[ind][2]) + 1;

                        for (int jj = 0; jj < prodConns; jj++) { // Horizontal well, increment in y
                           int k = wellDefs[ prodIndxs[ii] ][3] ;
                           if (wellDefs[ind ][1] <= wellDefs[ind ][2]) {
                              comp << "PRD"<< ii+1 <<"   " << wellDefs[ prodIndxs[ii] ][0];
                              comp << "    " << wellDefs[ prodIndxs[ii] ][1] + jj << "  " << k << "  " << k << " SHUT  *  *  4*  Y/" << endl;

                           } else {
                              comp << "PRD"<< ii+1 <<"   " << wellDefs[ prodIndxs[ii] ][0];
                              comp << "    " << wellDefs[ prodIndxs[ii] ][2] + jj << "  " << k << "  " << k << " SHUT  *  *  4*  Y/" << endl;

                           }
                        }
                        specs << "PRD"<< ii+1 <<"  PRD  " << wellDefs[ prodIndxs[ii] ][0];
                        specs << "    " << wellDefs[ prodIndxs[ii] ][1]<< "  4* HALT NO / " << endl;
                        //specs << "    " << wellDefs[ prodIndxs[ii] ][1]<< "  4* SHUT NO / " << endl;
                }
            }
            for (int ii = 0; ii < ninj; ii++) {
                if (wellBinVars[injIndxs[ii]] < -zeroEps) {
                        int ind = injIndxs[ii];
                        if (wellDefs[ind ][2] <= wellDefs[ind ][3]) {
                           comp << "INJ"<< ii+1 <<"   " << wellDefs[ injIndxs[ii] ][0];
                           comp << "    " << wellDefs[ ind ][1]<< "  " << wellDefs[ind ][2] << "  "<< wellDefs[ind ][3] << "  SHUT  *  *  4*  Z/" << endl;
                        } else {
                           comp << "INJ"<< ii+1 <<"   " << wellDefs[ injIndxs[ii] ][0];
                           comp << "    " << wellDefs[ind ][1]<< "  " << wellDefs[ind ][3] << "  "<< wellDefs[ind ][2] << "  SHUT  *  *  4*  Z/" << endl;
                        }

                        specs << "INJ"<< ii+1 <<"  INJ  " << wellDefs[ injIndxs[ii] ][0];
                        specs << "    " << wellDefs[ injIndxs[ii] ][1]<< " *  4*  NO /" << endl;

                        //stre << "INJ"<< ii+1 <<"   0.0 0.0 1.0 /" << endl;
                        stre << "INJ"<< ii+1 <<"  1.0  0.0 /" << endl;
                }
            }
	}
/*     =========================================================================================== */
        comp << "/" << endl;
        comp.close();
        specs << "/" << endl;
        specs.close();
        stre << "/" << endl;
        stre.close();

}

/* Function to read in vector of permeabilities from file */
void permFileReader (string& inputFilename, vector<double>& permVec)
{
    ifstream permin(inputFilename.c_str());

    string junk;
    //permin >> junk; // Skip first line

    for (int ii = 0; ii < permVec.size(); ii++) {
	permin >> permVec[ii];
    }

    permin.close();
    return;
}

/* Function to compute ro */
double  ComputeRo(double KX, double KY, double KZ, char direction, double DX,  double DY, double DZ)
{
	double ro = 0;
	if (direction == 'X')
		ro = 0.28*sqrt((sqrt(KY/KZ)*DZ*DZ + sqrt(KZ/KY)*DY*DY))/(pow(KY/KZ,0.25) + pow(KZ/KY,0.25));
	else
	if(direction == 'Y' )
		ro = 0.28*sqrt((sqrt(KZ/KX)*DX*DX + sqrt(KX/KZ)*DZ*DZ))/(pow(KZ/KX,0.25) + pow(KX/KZ,0.25));
	else
	if(direction == 'Z' )
		ro = 0.28*sqrt((sqrt(KY/KX)*DX*DX + sqrt(KX/KY)*DY*DY))/(pow(KY/KX,0.25) + pow(KX/KY,0.25));
		
	return ro;
}

/* Function that returns the number injectors, producers and control steps */
void getWellControlInfo(string optFilename, int& nInj, int& nProd, int& nSteps)
{
    FILE *input;
    char rmstr[1500];
    int nSkip[2] = { 12, 3 };   // Keywords to skip, first 12 before reading wells, then 3 more before reading control steps
    input = fopen("optimSetup.in","r");

    for (int i = 0; i < nSkip[0]; i++) {
        fgets(rmstr, 1500, input);      // Skip keyword
        fgets(rmstr, 1500, input);      // and value
    }
    fgets(rmstr, 1500, input);  // NPROD keyword
    fscanf(input, "%d\n", &nProd);
    fgets(rmstr, 1500, input);  // NINJ  keyword
    fscanf(input, "%d\n", &nInj);

    for (int i = 0; i < nSkip[1]; i++) {
        fgets(rmstr, 1500, input);      // Skip keyword
        fgets(rmstr, 1500, input);      // and value
    }
    fgets(rmstr, 1500, input);  // TSTEPS keyword
    fscanf(input, "%d\n", &nSteps);

    fclose(input);
//    cout << endl << "NPROD = " << nProd << ", NINJ = " << nInj << ", and TSTEPS = " << nSteps << endl;
}
