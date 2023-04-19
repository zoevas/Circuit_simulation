#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parser.h"
#include "mna.h"


//#define PRINT_RESULTS


int main(int argc, char *argv[]){
	int i, j, pos;
	double val, time;
	int tmp_N;
	struct element_t *elem;


	// parse the input file
	if(argc>1 && parse(argv[1])==PARSE_NOT_OK){
		printf("Parsing failed\n");
		return;
		}
	else if(argc<=1){
		printf("input file not specified\n");
		return;
		}
	printf("Parsing %d elements \t\t[OK]\n", elements_count);	


	// rename the nodes so they can be used for DC Modified Nodal Analysis
	nodeMapping(MNA_DC_ANALYSIS);
	tmp_N = N;
	printf("DC node mapping \t\t[OK]\n");

	// solve for DC analysis
	if(optionSPARSE){
		initLinearSystemSparse();
		calculateAbSparse();
		printf("Calculating A sparse matrix \t[OK]\n");
	
		if(optionITER){
			if(optionSPD){
				printf("Using sparse CG\n");
				calculateCGSparse();
                                }
			else{
				printf("Using sparse Bi-CG\n");
				calculateBiCGSparse();
				}
			dc_result = x;
			}
		else if(optionSPD){
			printf("Using sparse Cholesky\n");
			calculateCholeskySparse();
			solveAforCholSparse();
			dc_result = b;
			}
		else{
			printf("Using sparse LU\n");
			calculateLUSparse();
			solveASparse();
			dc_result = b;
			}
			
		}
	else{
		initLinearSystem();
		calculateAb();
		printf("Calculating A matrix \t\t[OK]\n");

		if(optionITER){
			if(optionSPD){
				printf("Using CG\n");
				calculateCG();
				}
			else{
				printf("Using Bi-CG\n");
				calculateBiCG();
				}			
			}
		else if(optionSPD){
			printf("Using Cholesky\n");
			calculateCholesky();
			solveAforChol();
			}
		else{
			printf("Using LU\n");
			calculateLU();
			solveA();
			}
		dc_result = x;
		}

#ifdef PRINT_RESULTS
	printVector(dc_result);
#endif

	// write the DC results to a file
	writeDCResultsToFile();

	// write the results of DC sweeps (if any) to a file
	writeDCsweepToFile();
	
	// free the recources used for the linear system
	if(optionSPARSE) freeLinearSystemSparse();
	else freeLinearSystem();


	// solve for transient analysis
	if(dotTRANflag){

		// rename the nodes so they can be used for Transient Modified Nodal Analysis
		nodeMapping(MNA_TRANSIENT_ANALYSIS);
		if(tmp_N==N){
			same_N = 1;
			printf("DC and TRAN use the same N\n");
			}
		else same_N = 0;
		printf("Transient node mapping \t\t[OK]\n");

		if(optionSPARSE){
			initLinearSystemTransSparse();
			calculateGCSparse();
			printf("Calculating C and G matrices \t[OK]\n");
			
			if(optionMETHOD == METHOD_BE){
				printf("Using sparse BE\n");
				calculateBE_ASparse(dotTRAN.time_step);
				calculateLUSparse();
				for(time=0.0; time<dotTRAN.fin_time; time+=dotTRAN.time_step){
					calculateBE_bSparse(time, dotTRAN.fin_time);
					solveASparse();
					writeTRANResultsToFile(time);
#ifdef PRINT_RESULTS
					printf("    --- time: %lf ---\n", time);
					printVector(b);
					printf("\n");
#endif
					}
				}
			else if(optionMETHOD == METHOD_TR){
				printf("Using sparse TR\n");
				calculateTR_ASparse(dotTRAN.time_step);
				calculateLUSparse();
				for(time=0.0; time<dotTRAN.fin_time; time+=dotTRAN.time_step){
					calculateTR_bSparse(time, dotTRAN.fin_time);
					solveASparse();
					writeTRANResultsToFile(time);
#ifdef PRINT_RESULTS
					printf("    --- time: %lf ---\n", time);
					printVector(b);
					printf("\n");
#endif					
					}
				}
			freeLinearSystemTransSparse();
			}
		else{
			initLinearSystemTrans();
			calculateGC();
			printf("Calculating C and G matrices \t[OK]\n");
						
			if(optionMETHOD == METHOD_BE){
				printf("Using BE\n");
				calculateBE_A(dotTRAN.time_step);
				calculateLU();
				for(time=0.0; time<dotTRAN.fin_time; time+=dotTRAN.time_step){
					calculateBE_b(time, dotTRAN.fin_time);
					solveA();
					writeTRANResultsToFile(time);
#ifdef PRINT_RESULTS
					printf("    --- time: %lf ---\n", time);
					printVector(x);
					printf("\n");
#endif					
					}
				}
			else if(optionMETHOD == METHOD_TR){
				printf("Using TR\n");
				calculateTR_A(dotTRAN.time_step);
				calculateLU();
				for(time=0.0; time<dotTRAN.fin_time; time+=dotTRAN.time_step){
					calculateTR_b(time, dotTRAN.fin_time);
					solveA();
					writeTRANResultsToFile(time);
#ifdef PRINT_RESULTS
					printf("    --- time: %lf ---\n", time);
					printVector(x);
					printf("\n");
#endif
					}
				}
			freeLinearSystemTrans();
			}
		printf("write TRAN results to file \t[OK]\n");
		}
	
	return;
	}

