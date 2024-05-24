/*
 * Code postFERN.cpp to run analyses on FERN outputs 
 * 
 * Compile with:
 *     gcc postFERN.cpp -o postFERN -lgsl -lgslcblas -lm -lstdc++
 * 
 * Resulting compiled code can be executed with
 * 
 *     ./postFERN | tee temp.txt
 * 
 *  
 * ----------------------------------
 * To set up a specific calculation:
 * ----------------------------------
 * 
 * 1. Change values of ISOTOPES and SIZE to match that of the network
 * 2. Double check the number of plot steps (lines in output file)
 * 3. Change the values of t0index and tEQindex to match the indice of the calculation being analysed
 * 4. Change input files in main to match the output files from explicitMatrix.cpp
 *    - Input files are found in main
 *	  - Data file for beta decay detection is also in main
 *
 *
 * ------------------------------------------------------------------------------------------
 * SOME SAMPLE NETWORKS:
 * 
 * Network   ISOTOPES     SIZE     networkFile[]               rateLibraryFile[]
 * 3=alpha          3        8     data/network_3alpha.inp     data/rateLibrary_3alpha.data
 * 4-alpha          4       14     data/network_4alpha.inp     data/rateLibrary_4alpha.data
 * alpha           16       48     data/network_alpha.inp      data/rateLibrary_alpha.data
 * pp               7       28     data/network_pp.inp         data/rateLibrary_pp.data
 * main cno         8       22     data/network_cno.inp        data/rateLibrary_cno.data
 * full cno        16      134     data/network_cnoAll.inp     data/rateLibrary_cnoAll.data
 * 48              48      299     data/network_48.inp         data/rateLibrary_48.data
 * 70(C-O)         70      598     data/network_70.inp         data/rateLibrary_70.data
 * 70(4He)         70      598     data/network_70_alpha.inp   data/rateLibrary_70.data
 * 116            116     1135     data/network_116.inp        data/rateLibrary_116.data
 * nova134        134     1566     data/network_nova134.inp    data/rateLibrary_nova134.data
 * 150 (12C-16O)  150     1604     data/network_150.inp        data/rateLibrary_150.data
 * 150 (solar)    150     1604     data/network_150_solar.inp  data/rateLibrary_150.data
 * 194            194     2232     network_194.inp             data/rateLibrary_194.data
 * 268            268     3175     network_268.inp             data/rateLibrary_268.data
 * 365 (12C-16O)  365     4395     data/network_365.inp        data/rateLibrary_365.data
 * 365 (solar)    365     4395     data/network_365_solar.inp  data/rateLibrary_365.data
 * tidalSN_alpha   16       48     data/network_alpha_he4.inp  data/rateLibrary_alpha.data
 * big bang         8       64     data/network_bigbang.inp    data/rateLibrary_bigbang.data
 * 28              28      104     data/network_28.inp         data/rateLibrary_28.data
 * 30P             47      283     data/network_test30P.inp    data/rateLibrary_test30P.data
 * ------------------------------------------------------------------------------------------
 *
 * 
 * AUTHORS:
 * ---------------
 * Nick Brey
 *
 * ----------------
 *
 *		TO DO LIST:
 *	___________________
 *
 *	1. Fix the errors that occur when nan is in an output file (anything larger than alpha network)
 * 		-- in matlab the nan values were set to = 0. gdb shows the nan in the file is an error in writing out (actual X values are of order 1e-18 or lower)
 *  2. Find a way to improve where t0 and tEQ are found (this would remove the need to implement a manual entry for the times to calculate over)
 *		-- t0 is set to step 0 so that could stay, but tEQ could be automated (beta decays could complicate this)
 *		-- maybe look at the fraction plot to see the max % PE/ASY(or QSS) obtain and use the 5th/6th column as a indicator where to set tEQ
 *
 */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <stdbool.h>
#include <ctime>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <tuple>

using namespace std;
//using std::string;

#define ISOTOPES 16                  // Max isotopes in network (e.g. 16 for alpha network)
#define SIZE 134					 // The number of reactions in the network
#define plotSteps 200                // Number of plot output steps (usually set to 200)

#define t0index 0                    // The plotStep # for where the error starts, (where RMS[i] > 0) 
#define tEQindex 175                 // The plotStep # for where the RMS calc will end (where RMS[i] ~0)
#define calcs 4                      // The number of calculations being used for uncertainty calculations (typically 4 between the Fast, Intermediate, Accurate, and Reference cases) 


char fastFile[] = "gnu_out/dataFiles/Alphafast1.data";
char intFile[] = "gnu_out/dataFiles/AlphaInt1.data";
char accFile[] = "gnu_out/dataFiles/AlphaAcc1.data";
char refFile[] = "gnu_out/dataFiles/AlphaRef1.data";

char rateFile[] = "data/rateLibrary_cnoAll.data";

bool DetectBeta = true; 			// Beta decay detection





// The methods and functions are listed below from 1-3
// Main function is at the bottom where the input file names will be edited
//Beta Decay detecttion is included: be sure to use the correct data file above
/*_______________________________________________________________________________________________________________________________________________________
											1: Report the maximum error for each Isotope in the network
Take in the residual arrays
Separate into ISOTOPE number of arrays
Find the max value in each array
_________________________________________________________________________________________________________________________________________________________*/
double method1(std::vector<std::vector<double>>& resF, const std::vector<std::vector<double>>& resI, const std::vector<std::vector<double>>& resA){


 // Iterate over each column
    for (int j = 0; j < ISOTOPES; j++) {
        // Initialize max value for current column
        double maxValsF = resF[0][j];
        double maxValsI = resI[0][j];
        double maxValsA = resA[0][j];

        // Iterate over each row in current column
        for (int i = 0; i < plotSteps; i++) {
            // Update max value if current element is greater
            if (resF[i][j] > maxValsF) {
                maxValsF = resF[i][j];
            }
            if (resI[i][j] > maxValsI) {
                maxValsI = resI[i][j];
            }
            if (resA[i][j] > maxValsA) {
                maxValsA = resA[i][j];
            }
        }
        
        // Print max value for current column
//        printf("FAST: max value for ISO[%d] = %f\n",j,maxValsF);
//        printf("INT: max value for ISO[%d] = %f\n",j,maxValsI);
//        printf("ACC: max value for ISO[%d] = %f\n",j,maxValsA);
    }
return 0;

}

// FUNCTION FOR FINDING TEQ and T0
double getTEquil(std::vector<std::vector<double>>& Xf){

  vector<vector<double>> dX(plotSteps, vector<double>(ISOTOPES)); 
  double sumdX[plotSteps] = {0}; 

  bool equil = false;
  double EQtol = 1e-4;
  int Teq;




    for(int i=1; i < plotSteps; i++){ // start at i = 2 to allow the first step to have been completed.
      for(int j=0; j < ISOTOPES; j++){
    
      	 dX[i][j] = Xf[i][j] - Xf[i-1][j];
       	 sumdX[i] += dX[i][j];

      } // end j-loop

//          printf("Sum of dx[%d] = %f\n",i,abs(sumdX[i]));

    } // END i-loop

  int ndX = sizeof(sumdX) / sizeof(sumdX[0]);
  double* maxValue1 = std::max_element(sumdX, sumdX + ndX);
  double maxdX = *maxValue1;
  double* TargetdX = find(&sumdX[0], sumdX + ndX, maxdX);
  int maxIndexdX = TargetdX - sumdX;


  for(int k = maxIndexdX; k < plotSteps; k++){

  		if(abs(sumdX[k]) < EQtol){

  			equil == true;
  			Teq = k;
  			break;
  		}

  }

  printf("The index with the max dX = %d and has value %f\n",maxIndexdX, maxdX);
  printf("Teq Index = %d\n",Teq);

    return Teq;
 } // END TEquil function


double getT0(std::vector<std::vector<double>>& Xf){

  vector<vector<double>> dX2(plotSteps, vector<double>(ISOTOPES)); 
  double sumdX2[plotSteps] = {0}; 

  double T0tol = 1e-5;
  int Tzero;

  bool diverge = false;

    for(int i=1; i < plotSteps; i++){ // start at i = 2 to allow the first step to have been completed.
      for(int j=0; j < ISOTOPES; j++){
    
       dX2[i][j] = Xf[i][j] - Xf[i-1][j];
       sumdX2[i] += dX2[i][j];
      } // end j-loop
    } // end i-loop

	for(int k = 1; k < plotSteps; k++){

  		if(sumdX2[k] > T0tol){

  			diverge == true;
  			Tzero = k;
  			break;
  		}

  	}	

//  printf("The index with the max dX = %d and has value %f\n",maxIndexdX, maxdX);
  printf("T0 Index = %d\n",Tzero);

    return Tzero;
 }


/*________________________________________________________________________________________________________________________________________________________
											2: Find the time where the error is the highest and report the RMS at that time

Sum the residuals at each step.
Find the step that has the highest sum
At that step, calculate the RMS error
This gives idea of accuracy across the network as a whole at a time where error is highest.
Can use RMS values to describe accuracy with single value 

 *** NOTE *** This method asssumes a single peak of the error 
 If the error is more spread out (not a single peak) then use method 3 where the RMS is calculated over the course of the whole calculation until equilibrium is reached
__________________________________________________________________________________________________________________________________________________________*/
//Take in residual arrays resF, resI, and resA
double method2(const std::vector<std::vector<double>>& resF, const std::vector<std::vector<double>>& resI, const std::vector<std::vector<double>>& resA){

   // Method 2 arrays
   vector<double> ResSqF(ISOTOPES);
   vector<double> ResSqI(ISOTOPES);
   vector<double> ResSqA(ISOTOPES);

// Sum the rows to find total error ar each plotStep output. Result should be 1x200 matrix (200 or # of plotsteps)
// with each entry containing the total error (sum of isotope residuals) for the itertation

double SumF[plotSteps] = {0};
double SumI[plotSteps] = {0};
double SumA[plotSteps] = {0};
    for(int i=0; i < plotSteps; i++){
        for(int j=0; j< ISOTOPES; j++){
                
          SumF[i] += resF[i][j];
          SumI[i] += resI[i][j];
          SumA[i] += resA[i][j];
        }

//            printf("sumF[%d] = %f\n",i,SumF[i]);
    }

// Now that the errors are summed at each step, find the step with the maximum error:          
// FAST
  int nF = sizeof(SumF) / sizeof(SumF[0]);
  double* maxValue1 = std::max_element(SumF, SumF + nF);
  double maxF = *maxValue1;
  double* TargetF = find(&SumF[0], SumF + nF, maxF);
  int maxIndexF = TargetF - SumF;
//  printf("The max error is %f at step %d\n",maxF,maxIndexF);

// INT
  int nI = sizeof(SumI) / sizeof(SumI[0]);
  double* maxValue2 = std::max_element(SumI, SumI + nI);
  double maxI = *maxValue2;
  double* TargetI = find(&SumI[0], SumI + nI, maxI);
  int maxIndexI = TargetI - SumI;
//  printf("The max error is %f at step %d\n",maxI,maxIndexI);

// ACC
  int nA = sizeof(SumA) / sizeof(SumA[0]);
  double* maxValue3 = std::max_element(SumA, SumA + nA);
  double maxA = *maxValue3;
  double* TargetA = find(&SumA[0], SumA + nA, maxA);
  int maxIndexA = TargetA - SumA;
//  printf("The max error is %f at step %d\n",maxA,maxIndexA);


// _______________ CALCULATE THE RMS OF THE ERROR ______________________//
// In j-loop, use the residual array to square each isotopic residual at the time of max_index
// Then sum over all isotopes (j's)

  double SumResSqF;
  double SumResSqI;
  double SumResSqA;

    for(int j=0; j < ISOTOPES; j++){

      ResSqF[j] = pow(resF[maxIndexF][j],2);
      ResSqI[j] = pow(resI[maxIndexI][j],2);
      ResSqA[j] = pow(resA[maxIndexA][j],2);
 
//      	printf("ResSqF[%d] = %f \n",j,ResSqI[j]);

      SumResSqF += ResSqF[j];
      SumResSqI += ResSqI[j];
      SumResSqA += ResSqA[j];
  } 

//  printf("Sum ResSqF = %f\n",SumResSqF);



// Now square root the Sum of the residuals squared for the RMS value at the point with the most error -- 
      double RMS1Fast = sqrt(SumResSqF);
      double RMS1Int = sqrt(SumResSqI);
      double RMS1Acc = sqrt(SumResSqA);

//        printf("RMS Error FAST = %f\n",RMS1Fast);
//        printf("RMS Error INT = %f\n",RMS1Int);
//        printf("RMS Error ACC = %f\n",RMS1Acc);


// return the RMS values for each calculation to main()      
	return 0;

}

/*_____________________________________________________________________________________________________________________________________________________________________________
												 3: Show the average error per timestep over the course of the calculation.

Take in residual arrays.
Calculate the RMS at each step
Find # of steps between where approximatiosn start (t0) and where equilibrium is achieved (tEQ)
Normalize the RMS over the time region of interest (t0->tEQ) 
Get single value for whole network for all time.
_______________________________________________________________________________________________________________________________________________________________________________*/
double method3(std::vector<std::vector<double>>& resF, const std::vector<std::vector<double>>& resI, const std::vector<std::vector<double>>& resA, const std::vector<double>& Tlog){


	vector<vector<double>> ResSqF(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> ResSqI(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> ResSqA(plotSteps, vector<double>(ISOTOPES));

    vector<double> SumResSqF(plotSteps);
    vector<double> SumResSqI(plotSteps);
    vector<double> SumResSqA(plotSteps);

    vector<double> RMSf(plotSteps);
    vector<double> RMSi(plotSteps);
    vector<double> RMSa(plotSteps);

    vector<double> R0f(plotSteps);
    vector<double> R0i(plotSteps);
    vector<double> R0a(plotSteps);

// Loop to square each residual in hte res[][] array -> then sum up all j-components to create a 1d array of Sum[plotSteps] -> sqrt the sum to get the RMS at each step
    for(int i=0; i < plotSteps; i++){
    	for(int j=0; j< ISOTOPES; j++){

    		ResSqF[i][j] = pow(resF[i][j],2);
      		ResSqI[i][j] = pow(resI[i][j],2);
      		ResSqA[i][j] = pow(resA[i][j],2);

//      	printf("ResSqF[%d][%d] = %f \n",i,j,ResSqF[i][j]);

      		SumResSqF[i] += ResSqF[i][j];
      		SumResSqI[i] += ResSqI[i][j];
        	SumResSqA[i] += ResSqA[i][j];

//      		printf("SumResSqF[%d] = %f \n",i,j,SumResSqF[i]);
    	}

    	RMSf[i] = sqrt(SumResSqF[i]);
    	RMSi[i] = sqrt(SumResSqI[i]);
    	RMSa[i] = sqrt(SumResSqA[i]);

//    	printf("RMSf[%d] = %f \n",i,RMSf[i]);

	} // close out i-loop

//                       CALC AVG ERROR PER TIME STEP
//______________________________________________________________________

 double t0 = Tlog[t0index];
 double teq = Tlog[tEQindex];

// printf("t0 = %f   teq = %f\n",t0, teq);
 
 double sumR0f = 0;
 double sumR0i = 0;
 double sumR0a = 0;

	for(int n=t0index; n < tEQindex+1; n++){

		R0f[n] = RMSf[n]*(Tlog[n]-Tlog[n-1]);
		R0i[n] = RMSi[n]*(Tlog[n]-Tlog[n-1]);
		R0a[n] = RMSa[n]*(Tlog[n]-Tlog[n-1]);

//		printf("R0f[%d] = %f  RMSf = %f  Tlog[n] = %f  Tlog[n-1] = %f\n",n, R0f[n], RMSf[n], Tlog[n], Tlog[n-1]);

		sumR0f += R0f[n];
		sumR0i += R0i[n];
		sumR0a += R0a[n];

//		printf("R0f[%d] = %f  RMSf = %f  and Sum = %f\n",n, R0f[n], RMSf[n], sumR0f);	
	}

// printf("sum of R0f = %f\n",sumR0f);

// NORMALIZE THE ERROR/TIME STEP 
 double NR0F = abs(sumR0f/(teq - t0));
 double NR0I = abs(sumR0i/(teq - t0));
 double NR0A = abs(sumR0a/(teq - t0));

 	printf("Normalized RMS-Fast  = %f sumR0f = %f t0 = %f   tEq= %f\n",NR0F,sumR0f, t0, teq);
	printf("Normalized RMS-Int  = %f sumR0i = %f t0 = %f   tEq= %f\n",NR0I,sumR0i, t0, teq);
	printf("Normalized RMS-Acc  = %f sumR0a = %f t0 = %f   tEq= %f\n",NR0A,sumR0a, t0, teq);	 	

return 0;

}

/*______________________________________________________________________________________________________________________________________________
								Uncertainty: Takes in all mass fractions from the different calculations and obtains uncertainty values



__________________________________________________________________________________________________________________________________________________*/
double uncertainty(std::vector<std::vector<double>>& Xf, const std::vector<std::vector<double>>& Xi, const std::vector<std::vector<double>>& Xa, const std::vector<std::vector<double>>& Xr){

//declare arrays/vectors
	vector<vector<double>> Mean(plotSteps, vector<double>(ISOTOPES));
	vector<vector<double>> DevF(plotSteps, vector<double>(ISOTOPES));
	vector<vector<double>> DevI(plotSteps, vector<double>(ISOTOPES));
	vector<vector<double>> DevA(plotSteps, vector<double>(ISOTOPES));
	vector<vector<double>> DevR(plotSteps, vector<double>(ISOTOPES));
	vector<vector<double>> IsoUncert(plotSteps, vector<double>(ISOTOPES));


// Calculate mean, deviation from the mean and uncertainty
for(int i=0; i < plotSteps; i++){
	for(int j=0; j < ISOTOPES; j++){

		Mean[i][j] = (Xf[i][j] + Xi[i][j] + Xa[i][j] + Xr[i][j])/ calcs;

		// Square of the deviation
		DevF[i][j] = pow(Xf[i][j] - Mean[i][j],2);
		DevI[i][j] = pow(Xi[i][j] - Mean[i][j],2);
		DevA[i][j] = pow(Xa[i][j] - Mean[i][j],2);
		DevR[i][j] = pow(Xr[i][j] - Mean[i][j],2);



		IsoUncert[i][j] = sqrt((DevF[i][j] + DevI[i][j] + DevA[i][j] + DevR[i][j]) / (calcs*(calcs-1)));
	}

}

 // Iterate over each column
    for (int j = 0; j < ISOTOPES; j++) {
        // Initialize max value for current column
        double maxUncert = IsoUncert[0][j];


        // Iterate over each row in current column
        for (int i = 0; i < plotSteps; i++) {
            // Update max value if current element is greater
            if (IsoUncert[i][j] > maxUncert) {
                maxUncert = IsoUncert[i][j];
            }
        }
        
        // Print max value for current column
//        printf("Uncertainty: max value for ISO[%d] = %f\n",j,maxUncert);

    }

return 0;
}

/*___________________________________________________________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________________________________________________________


																		BETA DECAY DETECTION

_____________________________________________________________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________________________________________________________*/
bool detectBetaPlus(const std::string& input_string, const std::string& sequence) {
    size_t found = input_string.find(sequence);
    return (found != std::string::npos);
}

bool detectBetaMinus(const std::string& input_string, const std::string& sequence) {
    size_t found = input_string.find(sequence);
    return (found != std::string::npos);
}

void BetaDecays(const char* filename){

    // Initialize variable for loop n: 0->SIZE
    int n = 0;

    // Define vectors that will be used to store the variables from the .data file while reading
    std::vector<int> RGclass(SIZE);
    std::vector<int> RGmemberIndex(SIZE);
    std::vector<int> reaclibClass(SIZE);
    std::vector<int> NumReactingSpecies(SIZE);
    std::vector<int> NumProducts(SIZE);
    std::vector<int> isEC(SIZE);
    std::vector<int> isReverseR(SIZE);
    std::vector<float> Prefac(SIZE);
    std::vector<float> Q(SIZE);
    std::vector<std::vector<int>> reactantZ(SIZE, std::vector<int>(4));
    std::vector<std::vector<int>> reactantN(SIZE, std::vector<int>(4));
    std::vector<std::vector<int>> productZ(SIZE, std::vector<int>(4));
    std::vector<std::vector<int>> productN(SIZE, std::vector<int>(4));
    std::vector<std::vector<int>> ReactantIndex(SIZE, std::vector<int>(4));
    std::vector<std::vector<int>> ProductIndex(SIZE, std::vector<int>(4));
    std::vector<float> P0(SIZE), P1(SIZE), P2(SIZE), P3(SIZE), P4(SIZE), P5(SIZE), P6(SIZE);

    std::vector<std::string> reactionType(SIZE);
    std::vector<int> BetaDecays;

    // Open file for reading
    std::ifstream fr(filename);
    if (!fr.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    // First while loop to loop through all the reactions: Starts @ n=0 and goes to n = SIZE
    while (n < SIZE) {
        int subindex = 0;

        // Lines 1-2 have a set number of input values to be read in
        //char reactionType[20];
        fr >> reactionType[n];
        //std::cout << reactionType[n] << std::endl;

        double reactionValues[10];
        for (int i = 0; i < 10; ++i)
            fr >> reactionValues[i];
        //std::cout << "Reaction values: ";

        //for (int i = 0; i < 10; ++i)
        //    std::cout << reactionValues[i] << " ";
        //    std::cout << std::endl;


        float reactionParams[7];
        for (int i = 0; i < 7; ++i)
            fr >> reactionParams[i];
            //std::cout << "Reaction parameters: ";

        //for (int i = 0; i < 7; ++i)
        //    std::cout << reactionParams[i] << " ";
        //    std::cout << std::endl;

        int numR = reactionValues[3];
        int numP = reactionValues[4];
        std::vector<int> ZreacNum(numR);
        std::vector<int> NreacNum(numR);
        std::vector<int> ZprodNum(numP);
        std::vector<int> NprodNum(numP);
        std::vector<int> reacVec(numR);
        std::vector<int> prodVec(numP);

        // Reacting Proton Number
        for (int i = 0; i < reactionValues[3]; ++i){
            fr >> ZreacNum[i];
        //    printf("i = %d\n",i);
        //    printf("Z-reac = %d\n",ZreacNum[i]);
        }

        // Reacting Neutron Number
        for (int i = 0; i < reactionValues[3]; ++i){
            fr >> NreacNum[i];

        //    printf("i = %d\n",i);
        //    printf("N-reac = %d\n",NreacNum[i]);
        }

        // Product Proton Number
        for (int i = 0; i < reactionValues[4]; ++i){
            fr >> ZprodNum[i];

        //    printf("i = %d\n",i);
        //    printf("Z-prod = %d\n",ZprodNum[i]);
        }

        // Product Neutron Number
        for (int i = 0; i < reactionValues[4]; ++i){
            fr >> NprodNum[i];

        //    printf("i = %d\n",i);
        //    printf("N-prod = %d\n",NprodNum[i]);
        }

        // 7th line (reactant species-vector index )
        for (int i = 0; i < reactionValues[3]; ++i){
            fr >> reacVec[i];

        //    printf("i = %d\n",i);
        //    printf("reactant Vector = %d\n",reacVec[i]);           
        }
       
        // 8th line (product species-vector index)
        for (int i = 0; i < reactionValues[4]; ++i){
            fr >> prodVec[i];

        //    printf("i = %d\n",i);
        //    printf("product Vector = %d\n",prodVec[i]);           
        }

        n++;

    } // END WHILE N < SIZE

// DETECT BETA DECAYS BY LOOKING FOR "e++" or "e+nubar" which are only in the Beta +/- decays respectively.
    std::string sequencePlus = "e++nu";
    std::string sequenceMinus = "e+nubar";        

    for(int i =0; i < SIZE; i++){
        std::string rxn = reactionType[i];


        if (detectBetaPlus(rxn, sequencePlus)) {
            std::cout << "Beta decay \"" << sequencePlus << "\" detected in reaction: " << reactionType[i] << " --- Index # : " << i << " (PLUS) " << std::endl;
        } 
        if (detectBetaMinus(rxn, sequenceMinus)) {
            std::cout << "Beta decay \"" << sequenceMinus << "\" detected in reaction: " << reactionType[i] << " --- Index # : " << i << " (MINUS) " << std::endl;
        } 
      //else {
      //    std::cout << "The sequence \"" << sequencePlus << " or " << sequenceMinus << "\" is not present in the line " << i << std::endl;
      //}
    }


}


/*________________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________________

																MAIN FUNCTION 

__________________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________________*/
// read in file (use edited versions of txt files for now)
int main() {
    //__________________________________________________________________//
    //                           Storage arrays 
    //__________________________________________________________________//

    // Time arrays
    vector<double> Tlog(plotSteps);
    vector<double> T(plotSteps);

    // Mass Fraction arrays
    vector<vector<double>> Xflog(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> Xilog(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> Xalog(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> Xrlog(plotSteps, vector<double>(ISOTOPES));
    // non-log values 
    vector<vector<double>> Xf(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> Xi(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> Xa(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> Xr(plotSteps, vector<double>(ISOTOPES));
    // Residual arrays
    vector<vector<double>> resF(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> resI(plotSteps, vector<double>(ISOTOPES));
    vector<vector<double>> resA(plotSteps, vector<double>(ISOTOPES));
    //_______________________________________________________________________//

    // READ IN FILES AND STORE/MODIFY ARRAYS
    // Open data files to be read in
    char *fastFilePtr = fastFile;
    char *intFilePtr = intFile;
    char *accFilePtr = accFile;
    char *refFilePtr = refFile;

    ifstream fid1(fastFilePtr); 
    ifstream fid2(intFilePtr);    
    ifstream fid3(accFilePtr);
    ifstream fid4(refFilePtr);

    int n = 0;                          // Loop variable set to 0
    int Columns = 7 + ISOTOPES;			// Number of columns in data file (7 columns listd before the isotope Mass Fractions are listed)
//___________________________________________________________________________________________________________________________________________________
    // While loop reads the data into variable arrays (F,I,A,R) 
    while (n < plotSteps){

        // Open/scan Fast calc
        vector<double> F(Columns);
        for (int k = 0; k < Columns; k++){
            fid1 >> F[k];

//       printf("F[%d][%d] = %f\n",n,k,F[k]);
        }

        // Open/scan Int calc
        vector<double> I(Columns);
        for (int k = 0; k < Columns; k++){
            fid2 >> I[k];

//       printf("I[%d][%d] = %f\n",n,k,I[k]);
        }
        
        // Open/scan Acc file
        vector<double> A(Columns);
        for (int k = 0; k < Columns; k++){
            fid3 >> A[k];

//      printf("A[%d][%d] = %f\n",n,k,A[k]);
        }
        
        // Open/scan Ref file
        vector<double> R(Columns);
        for (int k = 0; k < Columns; k++){
            fid4 >> R[k];

//       printf("F[%d][%d] = %f\n",n,k,R[k]);
        }
        
        Tlog[n] = R[0];
//        printf("Tlog[%d] = %f\n",n,Tlog[n]);

        // J-loop only writes the mass fractions into the X arrays
        // (r-ref, a-acc, m-med (or int)). Note values of data are in log(X)
        for (int j = 0; j < ISOTOPES; j++){
            Xflog[n][j] = F[7 + j];
            Xilog[n][j] = I[7 + j];
            Xalog[n][j] = A[7 + j];
            Xrlog[n][j] = R[7 + j];

//            printf("Xilog[%d][%d] = %f\n",n,j,Xilog[n][j]);
        }

         n++;
    } // END WHILE LOOP


/* ------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------  BEGIN ERROR ANALYSIS --------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------*/

// Mass Fraction array X, should now be a matrix [dataPoints x Isotopes] (with X in log values)
// Example if using triple alpha w/ 200 plot outputs array size [200 x 3]
// Array is currently in Log value, so to get true difference, get rid of logs: X = 10^(log(X)) for each array (Xf,Xi,Xa, Xr)
  for(int i=0; i < plotSteps; i++){
    for(int j=0; j< ISOTOPES; j++){

      Xf[i][j] = pow(10,Xflog[i][j]);
      Xi[i][j] = pow(10,Xilog[i][j]);
      Xa[i][j] = pow(10,Xalog[i][j]);
      Xr[i][j] = pow(10,Xrlog[i][j]);

//      printf("Xf[%d][%d} = %f\n",i,j,Xr[i][j]);

// Calculate residuals of each calculation: |X - Xr|
  	  resF[i][j] = abs(Xf[i][j] - Xr[i][j]);
      resI[i][j] = abs(Xi[i][j] - Xr[i][j]);
      resA[i][j] = abs(Xa[i][j] - Xr[i][j]);

//      printf("Xi[%d][%d} = %f\n",i,j,resI[i][j]);

    } //  end j (ISOTOPE) loop

      //T[i] = pow(10,Tlog[i]); // IF NEEDED (likely not used)

  } // end i (plotSteps) loop
//_____________________________________________________________________________________________________________________________

// Call the different methods 

method1(resF, resI, resA);

method2(resF, resI, resA);


// IF USING getT0 and getTEquil then make sure to pass those ints to method 3, these will be the bounds of the error calc. 
// These values won't have to be set up at the top of the code then

//int T0ind = getT0(Xf);
//int EQindex = getTEquil(Xf);
method3(resF, resI, resA, Tlog);

// Uncertainty will be its own function
uncertainty(Xf, Xi, Xa, Xr);


    fid1.close();
    fid2.close();
    fid3.close();
    fid4.close();

	if(DetectBeta == true){
		//BetaDecays("data/rateLibrary_cnoAll.data");
		char *rateFilePtr = rateFile;
    	BetaDecays(rateFilePtr);
	}

    return 0;
}