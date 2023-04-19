#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parser.h"
#include "mna.h"



/************************************************************************************************************************
 * This function renames all the node IDs so they are integers and sequential. In DC it open-circuits the Capacitors.	*
 * Saves the new node IDs to the specially reserved for this job slots of struct element: mna_terminal_a and		*
 * mna_terminal_b. Returns the nodes' count, what we call N in the MNA.							*
 ************************************************************************************************************************/
void nodeMapping(short dc_or_trans){
	int i, j, k, nodes_ctr=0;
	struct element_t *elem;

	// allocate memory for the array that all the useful node IDs will be stored in
	nodes = (char**)malloc(2*elements_count*sizeof(char*));
	MALLOC_RETURN_CHECK(nodes);	

	// always put node 0 to index 0 so find the first 0 node in the list and make nodes[0] point to it
	for(elem=list; elem!=NULL; elem=elem->next){
		if(strcmp(elem->terminal_a, "0")==0){ 
			nodes[nodes_ctr++] = elem->terminal_a;
			break;
			}
		else if(strcmp(elem->terminal_b, "0")==0){ 
			nodes[nodes_ctr++] = elem->terminal_b;
			break;
			}
		}
	
	// copy the node IDs to the nodes array ignoring the duplicates
	for(elem=list; elem!=NULL; elem=elem->next){
		short found_a=0, found_b=0;
		
		if(dc_or_trans==MNA_DC_ANALYSIS && elem->name[0] == 'C') continue;

		for(i=0; i<nodes_ctr; i++){
			if(strcmp(nodes[i], elem->terminal_a)==0) found_a = 1;
			if(strcmp(nodes[i], elem->terminal_b)==0) found_b = 1;
			}
		if(!found_a) nodes[nodes_ctr++] = elem->terminal_a;
		if(!found_b) nodes[nodes_ctr++] = elem->terminal_b;
		}

#ifdef ORDER_NODES
	// make lexicographical node ordering
	for(i=0; i<nodes_ctr-1; i++){
		for(j=i+1; j<nodes_ctr; j++){
			k=0;
			while(nodes[j][k]==nodes[i][k]){ k++; }
				if(nodes[j][k]<nodes[i][k]){
				char *a = nodes[i];
				nodes[i] = nodes[j];
				nodes[j] = a;
				}
			}
		}
#endif

 	// rename every node in the nodes array to it's corresponding index (eg.: rename node "n1_500_10" to "7" if nodes[7]="n1_500_10")
	for(elem=list; elem!=NULL; elem=elem->next){
		int new_node_a, new_node_b;		

		if(dc_or_trans==MNA_DC_ANALYSIS && elem->name[0] == 'C') continue;
	
		for(i=0; i<nodes_ctr; i++){
			if(strcmp(nodes[i], elem->terminal_a)==0) new_node_a = i;
			if(strcmp(nodes[i], elem->terminal_b)==0) new_node_b = i;
			}
		elem->mna_terminal_a = new_node_a;
		elem->mna_terminal_b = new_node_b;
		}

	nodes = (char**)realloc(nodes, nodes_ctr*sizeof(char*));	// keep only the neccessary memory for nodes array
	N = nodes_ctr;
	}




/************************************************************************************************************************
 * This function initializes the A*x=b linear system by allocating and seting to zero the appropriate memory. Also	*
 * calculates K, and initializes the pivoting vector P to P=[0,1,2,3,...].						*
 ************************************************************************************************************************/
void initLinearSystem(){
	int i, j;

	// K is the count of the voltage sources plus the count of the inductors (which are replaced by zero-value voltage sources)
	K = voltage_sources_count + inductors_count;

	A = mallocMatrix();
	b = mallocVector();
	x = mallocVector();

	// allocate memory for pivoting array P
	P = (int*)malloc((N-1+K)*sizeof(int));
	MALLOC_RETURN_CHECK(P);
	for(i=0; i<N-1+K; i++) P[i] = i;

	M = mallocVector();
	}




/************************************************************************************************************************
 * This function initializes the A*x=b linear system when sparse matrices must be used.					*
 ************************************************************************************************************************/
void initLinearSystemSparse(){
	int i, j;

	// K is the count of the voltage sources plus the count of the inductors (which are replaced by zero-value voltage sources)
	K = voltage_sources_count + inductors_count;

	b = mallocVector();
	x = mallocVector();
	M = mallocVector();
	b_dupl = mallocVector();
	}




/************************************************************************************************************************
 * This function initializes the linear system for Transient analysis.							*
 ************************************************************************************************************************/
void initLinearSystemTrans(){
	int i, j;
	FILE *f;

	// K is the count of the voltage sources plus the count of the inductors
	K = voltage_sources_count + inductors_count;

	G = mallocMatrix();
	C = mallocMatrix();
	A = mallocMatrix();
	hC = mallocMatrix();
	
	//b = mallocVector();
	e = mallocVector();
	x = mallocVector();
	if(dotTRANflag && optionMETHOD == METHOD_TR){
		e_prev = mallocVector();
		}
	for(i=0; i<N-1+K; i++) P[i] = i;

	// just so the last file is deleted and not appended
	f = fopen(DEFAULT_TRAN_RESULTS_FILENAME, "w");
	if(f==NULL){
		perror(NULL);
		exit(1);
		}
	fclose(f);
	}




/************************************************************************************************************************
 * This function initializes the linear system for sparse matrices and Transient analysis.				*
 ************************************************************************************************************************/
void initLinearSystemTransSparse(){
	int i, j;
	FILE *f;

	// K is the count of the voltage sources plus the count of the inductors (which are replaced by zero-value voltage sources)
	K = voltage_sources_count + inductors_count;

	b = mallocVector();
	x = mallocVector();
	M = mallocVector();
	//b_dupl = mallocVector();
	e = mallocVector();
	if(dotTRANflag && optionMETHOD == METHOD_TR){
		e_prev = mallocVector();
		}

	// just so the last file is deleted and not appended
	f = fopen(DEFAULT_TRAN_RESULTS_FILENAME, "w");
	if(f==NULL){
		perror(NULL);
		exit(1);
		}
	fclose(f);
	}




/************************************************************************************************************************
 * This function frees all the previously allocated memory regarding the A*x=b linear system.				*
 ************************************************************************************************************************/
void freeLinearSystem(){
	free(A);
	//free(b);
	//free(x);
	free(nodes);
	//free(P);
	free(M);
	}




/************************************************************************************************************************
 * This function frees all the previously allocated memory regarding the (sparse) A*x=b linear system. 			*
 ************************************************************************************************************************/
void freeLinearSystemSparse(){
	if(optionITER)free(x);
	//else free(b);
	free(nodes);
	free(M);
	cs_sfree(Ss);
	cs_nfree(Ns);
	//free(b_dupl);
	}




/************************************************************************************************************************
 * This function frees all the previously allocated memory regarding the transient linear system. 			*
 ************************************************************************************************************************/
void freeLinearSystemTrans(){
	free(G);
	free(C);
	free(A);
	free(hC);
	free(b);
	free(e);
	if(dotTRANflag && optionMETHOD == METHOD_TR){
		free(e_prev);
		}
	free(x);
	free(P);
	free(dc_result);
	}




/************************************************************************************************************************
 * This function frees all the previously allocated memory regarding the (sparse) transient liner system.		*
 ************************************************************************************************************************/
void freeLinearSystemTransSparse(){
	if(optionITER)free(x);
	//else free(b);
	free(nodes);
	free(M);
	cs_sfree(Ss);
	cs_nfree(Ns);
	free(b_dupl);
	free(e);
	if(dotTRANflag && optionMETHOD == METHOD_TR){
		free(e_prev);
		}
	
	}




/************************************************************************************************************************
 * This function just prints a given N-1+K by N-1+K matrix. Can be used to print matrix A.				*
 ************************************************************************************************************************/
void printMatrix(double **matrix){
	int i, j;
	for(i=0; i<N-1+K; i++){
		for(j=0; j<N-1+K; j++){
			if(matrix[i][j]>=0) printf("\t+%4.3lf", matrix[i][j]);
			else printf("\t%4.3lf", matrix[i][j]);
			}
		printf("\n");
		}
	}




/************************************************************************************************************************
 * This function just prints a given N-1+K vector. Can be used to print vector b or x.					*
 ************************************************************************************************************************/
void printVector(double *vector){
	int i;
	for(i=0; i<N-1+K; i++){
		if(vector[i]>=0) printf("\t+%4.10lf\n", vector[i]);
		else printf("\t%4.10lf\n", vector[i]);
		}
	}




/************************************************************************************************************************
 * This function calculates the matrix A and the vector b af the A*x=b linear system.					*
 ************************************************************************************************************************/
void calculateAb(){
	double g;	// the conductivity g of a resistor
	int pos, neg;	// the possitive and negative terminal respectively
	struct element_t *elem;
	int k_ctr=1;	// the voltage sources' counter. Goes from 0 to K

	// check if the circuit is appropriate for Modified Nodal Analysis
	for(elem=list; elem!=NULL; elem=elem->next)
		if(elem->name[0] != 'R' && elem->name[0] != 'L' && elem->name[0] != 'C' && elem->name[0] != 'V' && elem->name[0] != 'I'){
			printf("Circuit must be consisting only by passive elements in order to perform Modified Nodal Analysis.\n");
			exit(1);
			}

	// search the list and add every element's contribution to matrix A and vector b
	for(elem=list; elem!=NULL; elem=elem->next){

		pos = elem->mna_terminal_a;
		neg = elem->mna_terminal_b;

		if(elem->name[0] == 'R'){
			g = 1 / (elem->value);
			
			// the contribution of the resistors to matrix A
			if(pos!=0 && neg!=0){
				A[pos-1][neg-1] += (-1.0) * g;
				A[neg-1][pos-1] += (-1.0) * g;
				A[pos-1][pos-1] += g;
				A[neg-1][neg-1] += g;
				}
			else if(neg==0) A[pos-1][pos-1] += g;
			else if(pos==0) A[neg-1][neg-1] += g;
			}

		else if(elem->name[0] == 'V'){
			
			// the contribution of the voltage sources to matrix A
			if(neg!=0){
				A[neg-1][N-1+k_ctr-1] += (-1.0);
				A[N-1+k_ctr-1][neg-1] += (-1.0);
				}
			if(pos!=0){
				A[pos-1][N-1+k_ctr-1] += 1.0;
				A[N-1+k_ctr-1][pos-1] += 1.0;
				}

			// the contribution of the voltage sources to vector b
			b[N-1+k_ctr-1] += elem->value;

			k_ctr++;
			}

		else if(elem->name[0] == 'I'){

			// the contribution of the current sources to vector b
			if(pos!=0) b[pos-1] += (-1.0) * elem->value;
			if(neg!=0) b[neg-1] += elem->value;				
			}

		// replace Inductors with Voltage sources of 0V voltage
		else if(elem->name[0] == 'L'){
			
			// the contribution to matrix A
			if(neg!=0){
				A[neg-1][N-1+k_ctr-1] += (-1.0);
				A[N-1+k_ctr-1][neg-1] += (-1.0);
				}
			if(pos!=0){
				A[pos-1][N-1+k_ctr-1] += 1.0;
				A[N-1+k_ctr-1][pos-1] += 1.0;
				}

			k_ctr++;

			// since they are of 0 value, they do not contribute to vector b		
			}
		}//end for
	}




/************************************************************************************************************************
 * This function calculates the G and C matrices for transient analysis.						*
 ************************************************************************************************************************/
void calculateGC(){
	double g;	// the conductivity g of a resistor
	double c;	// the capacity c of a capacitor
	int pos, neg;	// the possitive and negative terminal respectively
	struct element_t *elem;
	int k_ctr=1;	// the voltage sources' counter. Goes from 0 to K

	// check if the circuit is appropriate for Modified Nodal Analysis
	for(elem=list; elem!=NULL; elem=elem->next)
		if(elem->name[0] != 'R' && elem->name[0] != 'L' && elem->name[0] != 'C' && elem->name[0] != 'V' && elem->name[0] != 'I'){
			printf("Circuit must be consisting only by passive elements in order to perform Modified Nodal Analysis.\n");
			exit(1);
			}

	// search the list and add every element's contribution to matrix C ang G and vector e
	for(elem=list; elem!=NULL; elem=elem->next){

		pos = elem->mna_terminal_a;
		neg = elem->mna_terminal_b;

		if(elem->name[0] == 'R'){
			g = 1 / (elem->value);
			
			// the contribution of the resistors to matrix G
			if(pos!=0 && neg!=0){
				G[pos-1][neg-1] += (-1.0) * g;
				G[neg-1][pos-1] += (-1.0) * g;
				G[pos-1][pos-1] += g;
				G[neg-1][neg-1] += g;
				}
			else if(neg==0) G[pos-1][pos-1] += g;
			else if(pos==0) G[neg-1][neg-1] += g;
			}

		else if(elem->name[0] == 'V'){
			
			// the contribution of the voltage sources to matrix G
			if(neg!=0){
				G[neg-1][N-1+k_ctr-1] += (-1.0);
				G[N-1+k_ctr-1][neg-1] += (-1.0);
				}
			if(pos!=0){
				G[pos-1][N-1+k_ctr-1] += 1.0;
				G[N-1+k_ctr-1][pos-1] += 1.0;
				}

			k_ctr++;
			}

		else if(elem->name[0] == 'L'){
			
			// the contribution of the inductors to matrices G and C
			if(neg!=0 && pos!=0){
				C[N-1+k_ctr-1][N-1+k_ctr-1] += (-1.0) * elem->value;
				}

			if(neg!=0){
				G[neg-1][N-1+k_ctr-1] += (-1.0);
				G[N-1+k_ctr-1][neg-1] += (-1.0);				
				}
			if(pos!=0){
				G[pos-1][N-1+k_ctr-1] += 1.0;
				G[N-1+k_ctr-1][pos-1] += 1.0;
				}

			k_ctr++;
			}

                else if(elem->name[0]=='C'){
			c = elem->value;

			// the contribution of the capacitors to matrix C
			if(neg!=0 && pos!=0){
				C[neg-1][pos-1] += (-1.0) * c;
				C[pos-1][neg-1] += (-1.0) * c;
				C[pos-1][pos-1] += c;
				C[neg-1][neg-1] += c;
				}
			}
	     }//end for
	}




/************************************************************************************************************************
 * This function calculates the matrix A of the transient linear system, using Backward Euler method.			*
 ************************************************************************************************************************/
void calculateBE_A(double h){
	int i;

	for(i=0; i<N-1+K; i++) setToZero(A[i]);

	multMatrixConstant(C, 1.0/h, hC);

	addMatrices(G, hC, A);
	
	if(same_N) assignVectorToVector(x, dc_result);
	}




/************************************************************************************************************************
 * This function calculates the vector b of the transient linear system, using Backward Euler method.			*
 ************************************************************************************************************************/
void calculateBE_b(double time, double final_time){
	int pos, neg;	// the positive and negative terminal respectively
	int i;
	int k_ctr=1;
	struct element_t *elem;
	double value;    

	setToZero(e);

	// calculate e(t) which changes at every step of transient analysis
	for(elem=list; elem!=NULL; elem=elem->next){
		pos = elem->mna_terminal_a;
		neg = elem->mna_terminal_b;
		if(elem->transient_source!=NULL){
			switch (elem->transient_source->type){
				case TRANSIENT_TYPE_EXP:
					value = EXP_Comp(elem->transient_source->un.transient_exp.i1, elem->transient_source->un.transient_exp.i2, elem->transient_source->un.transient_exp.td1, elem->transient_source->un.transient_exp.tc1, elem->transient_source->un.transient_exp.td2, elem->transient_source->un.transient_exp.tc2, time, final_time);         
					break;		
				case TRANSIENT_TYPE_SIN:	
					value = SIN_Comp(elem->transient_source->un.transient_sin.i1, elem->transient_source->un.transient_sin.ia, elem->transient_source->un.transient_sin.fr, elem->transient_source->un.transient_sin.td, elem->transient_source->un.transient_sin.df, elem->transient_source->un.transient_sin.ph, time, final_time);
					break;
				case TRANSIENT_TYPE_PULSE:
					value = PULSE_Comp(elem->transient_source->un.transient_pulse.i1, elem->transient_source->un.transient_pulse.i2, elem->transient_source->un.transient_pulse.td, elem->transient_source->un.transient_pulse.tr, elem->transient_source->un.transient_pulse.tf, elem->transient_source->un.transient_pulse.pw, elem->transient_source->un.transient_pulse.per, time, final_time);
					break;
				case TRANSIENT_TYPE_PWL:
					value = PWL_Comp(elem->transient_source->un.transient_pwl.t, elem->transient_source->un.transient_pwl.i, elem->transient_source->un.transient_pwl.n, time, final_time);
					break;
				}
			}
		else value = elem->value;

		if(elem->name[0] == 'V'){
			e[N-1+k_ctr-1] += value;
			k_ctr++;
			}
		else if(elem->name[0] == 'L') k_ctr++;
		else if(elem->name[0] == 'I'){
			if(pos!=0) e[pos-1] += (-1.0) * value;
			if(neg!=0) e[neg-1] += value;
			}
		}

	multMatrixVector(hC, x, b);

	addVectors(e, b, b);
	}




/************************************************************************************************************************
 * This function calculates the matrix A of the transient linear system, using Trapezodial method.			*
 ************************************************************************************************************************/
void calculateTR_A(double h){
	int i;

	for(i=0; i<N-1+K; i++) setToZero(A[i]);

	multMatrixConstant(C, 2.0/h, hC);

	addMatrices(G, hC, A);

	multMatrixConstant(C, (-2.0)/h, hC);

	addMatrices(G, hC, hC);

	multMatrixConstant(hC, (-1.0), hC);


	if(same_N) assignVectorToVector(x, dc_result);

	assignVectorToVector(e, b);
	}




/************************************************************************************************************************
 * This function calculates the vector d of the transient linear system, using Trapezodial method.			*
 ************************************************************************************************************************/
void calculateTR_b(double time, double final_time){
	int pos, neg;	// the positive and negative terminal respectively
	int i;
	int k_ctr=1;
	struct element_t *elem;
	double value;    

	assignVectorToVector(e_prev, e);
	setToZero(e);

	// calculate e(t) which changes at every step of transient analysis
	for(elem=list; elem!=NULL; elem=elem->next){
		pos = elem->mna_terminal_a;
		neg = elem->mna_terminal_b;
		if(elem->transient_source!=NULL){
			switch (elem->transient_source->type){
				case TRANSIENT_TYPE_EXP:
					value = EXP_Comp(elem->transient_source->un.transient_exp.i1, elem->transient_source->un.transient_exp.i2, elem->transient_source->un.transient_exp.td1, elem->transient_source->un.transient_exp.tc1, elem->transient_source->un.transient_exp.td2, elem->transient_source->un.transient_exp.tc2, time, final_time);         
					break;		
				case TRANSIENT_TYPE_SIN:	
					value = SIN_Comp(elem->transient_source->un.transient_sin.i1, elem->transient_source->un.transient_sin.ia, elem->transient_source->un.transient_sin.fr, elem->transient_source->un.transient_sin.td, elem->transient_source->un.transient_sin.df, elem->transient_source->un.transient_sin.ph, time, final_time);
					break;
				case TRANSIENT_TYPE_PULSE:
					value = PULSE_Comp(elem->transient_source->un.transient_pulse.i1, elem->transient_source->un.transient_pulse.i2, elem->transient_source->un.transient_pulse.td, elem->transient_source->un.transient_pulse.tr, elem->transient_source->un.transient_pulse.tf, elem->transient_source->un.transient_pulse.pw, elem->transient_source->un.transient_pulse.per, time, final_time);
					break;
				case TRANSIENT_TYPE_PWL:
					value = PWL_Comp(elem->transient_source->un.transient_pwl.t, elem->transient_source->un.transient_pwl.i, elem->transient_source->un.transient_pwl.n, time, final_time);
					break;
				default:
					printf("error in transient source type\n");
					exit(1);
				}
			}
		else value = elem->value;

		if(elem->name[0] == 'V'){
			e[N-1+k_ctr-1] += value;
			k_ctr++;
			}
		else if(elem->name[0] == 'L') k_ctr++;
		else if(elem->name[0] == 'I'){
			if(pos!=0) e[pos-1] += (-1.0) * value;
			if(neg!=0) e[neg-1] += value;
			}
		}

	multMatrixVector(hC, x, b);
	addVectors(e_prev, b, b);
	addVectors(e, b, b);
	}





/************************************************************************************************************************
 * This function calculates the G and C sparse matrices of the transient linear system.					*
 ************************************************************************************************************************/
void calculateGCSparse(){
	double g;	// the conductivity g of a resistor
	double c;	// the capacitance c of a capacitor
	int pos, neg;	// the possitive and negative terminal respectively
	struct element_t *elem;
	cs  *tran_Cs, *tran_Gs;
        int k_ctr=1;	// the voltage sources' counter. Goes from 0 to K
	int i;


	// check if the circuit is appropriate for Modified Nodal Analysis
	for(elem=list; elem!=NULL; elem=elem->next)
		if(elem->name[0] != 'R' && elem->name[0] != 'L' && elem->name[0] != 'C' && elem->name[0] != 'V' && elem->name[0] != 'I'){
			printf("Circuit must be consisting only by passive elements in order to perform Modified Nodal Analysis.\n");
			exit(1);
			}

	// find the number of non-zeros
	nzC = nzG = 0;
	for(elem=list; elem!=NULL; elem=elem->next){
		if(elem->name[0] == 'C'){
			if(elem->mna_terminal_a!=0 && elem->mna_terminal_b!=0) nzC+=4;
			else nzC+=1;
			}
		else if(elem->name[0] == 'L'){
			nzC+=1;
			if(elem->mna_terminal_a!=0 && elem->mna_terminal_b!=0) nzG+=4;
			else nzG+=1;
			}
		else if(elem->name[0] == 'R'){
			if(elem->mna_terminal_a!=0 && elem->mna_terminal_b!=0) nzG+=4;
			else nzG+=1;
			}
		else if(elem->name[0] == 'V'){
			if(elem->mna_terminal_a!=0 && elem->mna_terminal_b!=0) nzG+=4;
			else nzG+=1;
			}		
		}
	tran_Cs = cs_spalloc(N-1+K, N-1+K, nzC, 1, 1);
	tran_Gs = cs_spalloc(N-1+K, N-1+K, nzG, 1, 1);

	// search the list and add every element's contribution to matrix tran_Cs and tran_Gs
	for(elem=list; elem!=NULL; elem=elem->next){

		pos = elem->mna_terminal_a;
		neg = elem->mna_terminal_b;

		if(elem->name[0] == 'R'){
			g = 1 / (elem->value);
			
			if(pos!=0 && neg!=0){
				cs_entry(tran_Gs, pos-1, neg-1, (-1.0)*g);
                                cs_entry(tran_Gs, neg-1, pos-1, (-1.0)*g);
                                cs_entry(tran_Gs, pos-1, pos-1, g);
                                cs_entry(tran_Gs, neg-1, neg-1, g);
				}
			else if(neg==0) cs_entry(tran_Gs, pos-1, pos-1, g);
			else if(pos==0) cs_entry(tran_Gs, neg-1, neg-1, g);
			}

		else if(elem->name[0] == 'V'){
			
			if(neg!=0){
                                cs_entry(tran_Gs, neg-1, N-1+k_ctr-1, -1.0);
                                cs_entry(tran_Gs, N-1+k_ctr-1, neg-1, -1.0);
				}
			if(pos!=0){
                                cs_entry(tran_Gs, pos-1, N-1+k_ctr-1, 1.0);
                                cs_entry(tran_Gs, N-1+k_ctr-1, pos-1, 1.0);
				}

			k_ctr++;
			}

		else if(elem->name[0] == 'L'){
			
			if(neg!=0){
                                cs_entry(tran_Gs, neg-1, N-1+k_ctr-1, -1.0);
                                cs_entry(tran_Gs, N-1+k_ctr-1, neg-1, -1.0);
				}
			if(pos!=0){
				cs_entry(tran_Gs, pos-1, N-1+k_ctr-1, 1.0);
                                cs_entry(tran_Gs, N-1+k_ctr-1, pos-1, 1.0);
				}
			
			cs_entry(tran_Cs, N-1+k_ctr-1, N-1+k_ctr-1, (-1.0)*elem->value);

			k_ctr++;
			}

		else if(elem->name[0] == 'C'){
			c = elem->value;
			
			if(pos!=0 && neg!=0){
				cs_entry(tran_Cs, pos-1, neg-1, (-1.0)*c);
                                cs_entry(tran_Cs, neg-1, pos-1, (-1.0)*c);
                                cs_entry(tran_Cs, pos-1, pos-1, c);
                                cs_entry(tran_Cs, neg-1, neg-1, c);
				}
			else if(neg==0) cs_entry(tran_Cs, pos-1, pos-1, c);
			else if(pos==0) cs_entry(tran_Cs, neg-1, neg-1, c);
			}
		}//end for
	        
	//Convert the triplet tran_Cs and tran_Gs into a compressed column matrix
	tran_Csc = cs_compress(tran_Cs);
	tran_Gsc = cs_compress(tran_Gs);
	Cs = cs_compress(tran_Cs);
	cs_spfree(tran_Cs);
	cs_spfree(tran_Gs);
	cs_dupl(tran_Csc);
	cs_dupl(tran_Gsc);
	cs_dupl(Cs);
	cs_print(tran_Csc, "pin_tranC.txt", 0);
	cs_print(tran_Gsc, "pin_tranG.txt", 0);
	}




/************************************************************************************************************************
 * This function calculates the sparse matrix A of the transient linear system, using Backward Euler method.		*
 ************************************************************************************************************************/
void calculateBE_ASparse(double h){

	Cs = cs_add(tran_Gsc, tran_Csc, 1.0, 1.0/h);

	tran_Csc = cs_add(tran_Gsc, tran_Csc, 0.0, 1.0/h);
	
	if(same_N) assignVectorToVector(b, dc_result);
	} 




/************************************************************************************************************************
 * This function calculates vector b of the transient linear system, using Backward Euler method.			*
 ************************************************************************************************************************/
void calculateBE_bSparse(double time, double final_time){
	int pos, neg;	// the positive and negative terminal respectively
	int i;
	int k_ctr=1;
	struct element_t *elem;
	double value;    

	setToZero(e);

	// calculate e(t) which changes at every step of transient analysis
	for(elem=list; elem!=NULL; elem=elem->next){
		pos = elem->mna_terminal_a;
		neg = elem->mna_terminal_b;
		if(elem->transient_source!=NULL){
			switch (elem->transient_source->type){
				case TRANSIENT_TYPE_EXP:
					value = EXP_Comp(elem->transient_source->un.transient_exp.i1, elem->transient_source->un.transient_exp.i2, elem->transient_source->un.transient_exp.td1, elem->transient_source->un.transient_exp.tc1, elem->transient_source->un.transient_exp.td2, elem->transient_source->un.transient_exp.tc2, time, final_time);         
					break;		
				case TRANSIENT_TYPE_SIN:	
					value = SIN_Comp(elem->transient_source->un.transient_sin.i1, elem->transient_source->un.transient_sin.ia, elem->transient_source->un.transient_sin.fr, elem->transient_source->un.transient_sin.td, elem->transient_source->un.transient_sin.df, elem->transient_source->un.transient_sin.ph, time, final_time);
					break;
				case TRANSIENT_TYPE_PULSE:
					value = PULSE_Comp(elem->transient_source->un.transient_pulse.i1, elem->transient_source->un.transient_pulse.i2, elem->transient_source->un.transient_pulse.td, elem->transient_source->un.transient_pulse.tr, elem->transient_source->un.transient_pulse.tf, elem->transient_source->un.transient_pulse.pw, elem->transient_source->un.transient_pulse.per, time, final_time);
					break;
				case TRANSIENT_TYPE_PWL:
					value = PWL_Comp(elem->transient_source->un.transient_pwl.t, elem->transient_source->un.transient_pwl.i, elem->transient_source->un.transient_pwl.n, time, final_time);
					break;
				}
			}
		else value = elem->value;

		if(elem->name[0] == 'V'){
			e[N-1+k_ctr-1] += value;
			k_ctr++;
			}
		else if(elem->name[0] == 'L') k_ctr++;
		else if(elem->name[0] == 'I'){
			if(pos!=0) e[pos-1] += (-1.0) * value;
			if(neg!=0) e[neg-1] += value;
			}
		}

	multMatrixVectorSparse(tran_Csc, b, b);

	addVectors(e, b, b);
	}




/************************************************************************************************************************
 * This function calculates the sparse matrix A of the transient linear system, using Trapezodial method.		*
 ************************************************************************************************************************/
void calculateTR_ASparse(double h){
	
	Cs = cs_add(tran_Gsc, tran_Csc, 1.0, 2.0/h);

	tran_Csc = cs_add(tran_Gsc, tran_Csc, (-1.0), 2.0/h);

	if(same_N) assignVectorToVector(b, dc_result);

	assignVectorToVector(e, b_dupl);
	}




/************************************************************************************************************************
 * This function calculates vector b of the transient linear system, using Trapezodial method.				*
 ************************************************************************************************************************/
void calculateTR_bSparse(double time, double final_time){
	int pos, neg;	// the positive and negative terminal respectively
	int i;
	int k_ctr=1;
	struct element_t *elem;
	double value;    

	setToZero(e_prev);
	assignVectorToVector(e_prev, e);
	setToZero(e);

	// calculate e(t) which changes at every step of transient analysis
	for(elem=list; elem!=NULL; elem=elem->next){
		pos = elem->mna_terminal_a;
		neg = elem->mna_terminal_b;
		if(elem->transient_source!=NULL){
			switch (elem->transient_source->type){
				case TRANSIENT_TYPE_EXP:
					value = EXP_Comp(elem->transient_source->un.transient_exp.i1, elem->transient_source->un.transient_exp.i2, elem->transient_source->un.transient_exp.td1, elem->transient_source->un.transient_exp.tc1, elem->transient_source->un.transient_exp.td2, elem->transient_source->un.transient_exp.tc2, time, final_time);         
					break;		
				case TRANSIENT_TYPE_SIN:	
					value = SIN_Comp(elem->transient_source->un.transient_sin.i1, elem->transient_source->un.transient_sin.ia, elem->transient_source->un.transient_sin.fr, elem->transient_source->un.transient_sin.td, elem->transient_source->un.transient_sin.df, elem->transient_source->un.transient_sin.ph, time, final_time);
					break;
				case TRANSIENT_TYPE_PULSE:
					value = PULSE_Comp(elem->transient_source->un.transient_pulse.i1, elem->transient_source->un.transient_pulse.i2, elem->transient_source->un.transient_pulse.td, elem->transient_source->un.transient_pulse.tr, elem->transient_source->un.transient_pulse.tf, elem->transient_source->un.transient_pulse.pw, elem->transient_source->un.transient_pulse.per, time, final_time);
					break;
				case TRANSIENT_TYPE_PWL:
					value = PWL_Comp(elem->transient_source->un.transient_pwl.t, elem->transient_source->un.transient_pwl.i, elem->transient_source->un.transient_pwl.n, time, final_time);
					break;
				default:
					printf("error in transient source type\n");
					exit(1);
				}
			}
		else value = elem->value;

		if(elem->name[0] == 'V'){
			e[N-1+k_ctr-1] += value;
			k_ctr++;
			}
		else if(elem->name[0] == 'L') k_ctr++;
		else if(elem->name[0] == 'I'){
			if(pos!=0) e[pos-1] += (-1.0) * value;
			if(neg!=0) e[neg-1] += value;
			}
		}


	multMatrixVectorSparse(tran_Csc, b, b);

	addVectors(e_prev, b, b);

	addVectors(e, b, b);
	}




/************************************************************************************************************************
 * Computes the [transient_spec] EXP(i1 i2 td1 tc1 td2 tc2). The parameters time and final_time are required for 	*
 * finding which expression to be used. Returns the computed value. If the parameter time is out of time boundaries,	*
 * the function returns a negative value.										*
 ************************************************************************************************************************/
double EXP_Comp(double i1, double i2, double td1, double tc1, double td2, double tc2, double time, double final_time){
	double pow1, pow2;

	if(time>=0.0 && time<=td1) return i1;
	else if(time>td1 && time<=td2){
		pow1 = (-1.0)*(time-td1)/tc1;
		return (i1+(i2-i1)*(1.0-exp(pow1)));
		}
	else if(time>td2 && time<=final_time){
		pow1 = (-1.0)*(time-td1)/tc1;
		pow2 = (-1.0)*(time-td2)/tc2;
		return (i1+(i2-i1)*(exp(pow2)-exp(pow1)));
		}    
	else{
		printf("invalid time in SIN calculation\n");
		exit(1);
		}
	}




/************************************************************************************************************************
 * This function computes the [transient_spec] SIN(i1 ia fr td df ph). The parameters time and final_time are required	*
 * for finding which expression to be used. Returns the computed value. If the parameter time is out of time boundaries	*
 * the function returns a negative value.										*
 ************************************************************************************************************************/
double SIN_Comp(double i1, double ia, double fr, double td, double df, double ph, double time, double final_time){
	double pi = 3.1459;
	double pow;

	if(time>=0.0 && time<=td) return (i1+ia*sin(2*pi*ph/360));
	else if(time>td && time<=final_time){
		pow = (-1.0) * (time-td) * df;  
		return (i1+ia*sin(2*pi*fr*(time-td)+2*pi*ph/360)*exp(pow));
		}
	else{
 		printf("invalid time in EXP calculation\n");
		exit(1);
		}
	}




/************************************************************************************************************************
 * This function computes the [transient_spec] PWL(t1 v1 t2 v2 ... tn vn). The parameters time and final_time are	*
 * required for finding which expression to be used. It returns the computed value. If the parameter time is out of	*
 * time boundaries, the function returns a negative value.								*
 ************************************************************************************************************************/
double PWL_Comp(double *t, double *y, int n, double time, double final_time){
	int i;

	for(i=0; i<n; i++)
		if(t[i]<=time && time<=t[i+1])   
			return linearInterpolation(t[i], t[i+1], y[i], y[i+1], time);
	printf("invalid time in PWL calculation\n");
	exit(1);		
	}




/************************************************************************************************************************
 * This function computes the [transient_spec] PULSE(i1 i2 td tr tf pw per). The parameters time and final_time are	*
 * required for finding which expression to be used. It returns the computed value. If the parameter time is out of	*
 * time boundaries, the function returns a negative value.								*
 ************************************************************************************************************************/
double PULSE_Comp(double i1, double i2, double td, double tr, double tf, double pw, double per, double time, double final_time){
	int k;

	k = (int)(time/per);
	if(time>=k*per && time<(td+k*per)) return i1;
	else if(time>=(td+k*per) && time<(td+tr+k*per)) return linearInterpolation(td+k*per, td+tr+k*per, i1, i2, time);
	else if(time>=(td+tr+k*per) && time<(td+tr+pw+k*per)) return i2;
	else if(time>=(td+tr+pw+k*per) && time<(td+tr+pw+tf+k*per)) return linearInterpolation(td+tr+pw+k*per, td+tr+pw+tf+k*per, i2, i1, time);
	else if(time>=(td+tr+pw+tf+k*per) && time<=(td+per+k*per)) return i1;
	else printf("invalid time in PULSE calculation\n");
	exit(1);	
	}




/************************************************************************************************************************
 * Performs the addition of two matrices A and B, and returns the result to the Matrix res.				*
 ************************************************************************************************************************/
void addMatrices(double  **A, double **B, double **res){
	int i, j;
	
	for(i=0; i<N-1+K; i++)
		for(j=0; j<N-1+K; j++)
			res[i][j] = A[i][j] + B[i][j];
	}




/************************************************************************************************************************
 * Performs the multiplication of a matrix A and a constant, and returns the result to the Matrix res.			*	
 ************************************************************************************************************************/
void multMatrixConstant(double **A, double constant, double **res){
	int i, j;

	for(i=0; i<N-1+K; i++)
		for(j=0; j<N-1+K; j++)
			res[i][j] = A[i][j] * constant;
	}




/************************************************************************************************************************
 * Performs the addition of two vectors x and y, and returns the result to the vector res.				*		
 ************************************************************************************************************************/
void addVectors(double *x , double *y, double *res){
	int i;

	for(i=0; i<N-1+K; i++)
		res[i] = x[i] + y[i];
	}




/************************************************************************************************************************
 *  Assigns the elements of vector y to the vector x.									*		
 ************************************************************************************************************************/
void assignVectorToVector(double *x, double *y){
	int i;
	
	for(i=0; i<N-1+K; i++)
		x[i] = y[i];
	}



/************************************************************************************************************************
 * Given 2 points (x0,y0) and (x1,y0) and an x coordinate of a third point, it calculates the y coordinate of this	*
 * third point using Linear Interpolation. Returns this y coordinate.							*
 ************************************************************************************************************************/
double linearInterpolation(double x0, double x1, double y0, double y1, double x){
	return ((y1-y0)/(x1-x0)) * (x-x0) + y0;
	}




/************************************************************************************************************************
 * This function sets a vector to all zeros.										*
 ************************************************************************************************************************/
void setToZero(double *vector){
	int i;

	for(i=0; i<N-1+K; i++) vector[i] = 0.0;
	}




/************************************************************************************************************************
 * This function allocates memory for a N-1+K vector. Returns a pointer to the allocated memory.			*
 ************************************************************************************************************************/
double *mallocVector(){
	double *vector;
	vector = (double*)malloc((N-1+K)*sizeof(double));
	MALLOC_RETURN_CHECK(vector);
	setToZero(vector);
	return vector;
	}




/************************************************************************************************************************
 * This function allocates memory for a (N-1+K) * (N-1+K) matrix. Returns a pointer to the allocated memory.		*
 ************************************************************************************************************************/
double **mallocMatrix(){
	int i;
	double **matrix;
	matrix = (double**)malloc((N-1+K)*sizeof(double*));
	MALLOC_RETURN_CHECK(matrix);
	for(i=0; i<N-1+K; i++) matrix[i] = mallocVector();
	return matrix;
	}




/************************************************************************************************************************
 * This function calculates the sparse matrix Cs in compressed comumn form and the vector b af the A*x=b linear system.	*
 ************************************************************************************************************************/
void calculateAbSparse(){
	double g;	// the conductivity g of a resistor
	int pos, neg;	// the possitive and negative terminal respectively
	struct element_t *elem;
	cs  *As;
        int k_ctr=1;	// the voltage sources' counter. Goes from 0 to K
	int i;

	// check if the circuit is appropriate for Modified Nodal Analysis
	for(elem=list; elem!=NULL; elem=elem->next)
		if(elem->name[0] != 'R' && elem->name[0] != 'L' && elem->name[0] != 'C' && elem->name[0] != 'V' && elem->name[0] != 'I'){
			printf("Circuit must be consisting only by passive elements in order to perform Modified Nodal Analysis.\n");
			exit(1);
			}

	// find the number of non-zeros
	nz = 0;
	for(elem=list; elem!=NULL; elem=elem->next){
		if((elem->name[0] == 'R')||(elem->name[0] == 'V')){
			if(elem->mna_terminal_a!=0 && elem->mna_terminal_b!=0) nz+=4;
			else nz+=1;
			}
		}
	As = cs_spalloc(N-1+K, N-1+K, nz, 1, 1);

	// search the list and add every element's contribution to matrix As and vector b
	for(elem=list; elem!=NULL; elem=elem->next){

		pos = elem->mna_terminal_a;
		neg = elem->mna_terminal_b;

		if(elem->name[0] == 'R'){
			g = 1 / (elem->value);
			
			// the contribution of the resistors to matrix As
			if(pos!=0 && neg!=0){
				cs_entry(As, pos-1, neg-1, (-1.0)*g);
                                cs_entry(As, neg-1, pos-1, (-1.0)*g);
                                cs_entry(As, pos-1, pos-1, g);
                                cs_entry(As, neg-1, neg-1, g);
				}
			else if(neg==0) cs_entry(As, pos-1, pos-1, g);
			else if(pos==0) cs_entry(As, neg-1, neg-1, g);
			}

		else if(elem->name[0] == 'V'){
			
			// the contribution of the voltage sources to matrix As
			if(neg!=0){
                                cs_entry(As, neg-1, N-1+k_ctr-1, -1.0);
                                cs_entry(As, N-1+k_ctr-1, neg-1, -1.0);
				}
			if(pos!=0){
                                cs_entry(As, pos-1, N-1+k_ctr-1, 1.0);
                                cs_entry(As, N-1+k_ctr-1, pos-1, 1.0);
				}

			// the contribution of the voltage sources to vector b
			b[N-1+k_ctr-1] += elem->value;

			k_ctr++;
			}

		else if(elem->name[0] == 'I'){
			// the contribution of the current sources to vector b
			if(pos!=0) b[pos-1] += (-1.0) * elem->value;
			if(neg!=0) b[neg-1] += elem->value;				
			}

		// replace Inductors with Voltage sources of 0V voltage
		else if(elem->name[0] == 'L'){
			
			// the contribution to matrix As
			if(neg!=0){
                                cs_entry(As, neg-1, N-1+k_ctr-1, -1.0);
                                cs_entry(As, N-1+k_ctr-1, neg-1, -1.0);
				}
			if(pos!=0){
				cs_entry(As, pos-1, N-1+k_ctr-1, 1.0);
                                cs_entry(As, N-1+k_ctr-1, pos-1, 1.0);
				}

			k_ctr++;

			// since they are of 0 value, they do not contribute to vector b		
			}
		}//end for
	        
	//Convert the triplet As into a compressed column matrix c
	Cs = cs_compress(As);
	cs_spfree(As);
	cs_dupl(Cs);
	cs_print(Cs, "pinC.txt", 0);
	for(i=0; i<N-1+K; i++) b_dupl[i] = b[i];
	}




/************************************************************************************************************************
 * This function makes the LU factorization of matrix A. The results are stored back to the input matrix A where the 	*
 * the part over the diagonal is the upper-triangular matrix U and the part below the diagonal is the lower-triangular	*
 * matrix L. The diagonal belongs to U matrix as the diagonal of L is known (all diagonal elements of L are 1).		*
 * Performs partial pivoting and stores the pivoting results to the one-dimentional vector P.				*
 ************************************************************************************************************************/
void calculateLU(){
	int i, j, k, m, tmp_int;
	double x, *tmp_double;

	for(k=0; k<N-1+K; k++){
		x = fabs(A[k][k]);
		m=k;
		// find the largest m for pivoting		
		for(i=k+1; i<N-1+K; i++){
			if(fabs(A[i][k])>x) m=i;
			}
		if(A[m][k]==0){
			printf("Division by zero during LU factorization.\n");
			exit(1);
			}

		// make the pivot in matrix A and also in pivoting vector P			
		tmp_double = A[m];
		A[m] = A[k];
		A[k] = tmp_double;
		tmp_int = P[m];
		P[m] = P[k];
		P[k] = tmp_int;

		// make the LU calculations
		for(i=k+1; i<N-1+K; i++){
			A[i][k] = A[i][k]/A[k][k]; 
			}
		for(i=k+1; i<N-1+K; i++){
			for(j=k+1; j<N-1+K; j++){
				A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
	}




/************************************************************************************************************************
* Same as calculateLU but for sparse matrices.										*
************************************************************************************************************************/
void calculateLUSparse(){
	int i;
	Ss = cs_sqr(2, Cs, 0);
	Ns = cs_lu(Cs, Ss, 1);
	cs_spfree(Cs);
	}




/************************************************************************************************************************
 * Calculates the Cholesky factorization of matrix A. The results are writen back to matrix A. Does not involve vector	*
 * b. Also checks if the matrix is apropriate for the factorization (positive-definite matrix).				*
 ************************************************************************************************************************/
void calculateCholesky(){
	int i, j, k;
	double sum;

	/* //old code -- does not work properly
	for(k=0; k<N-1+K; k++){
		sum=0;
		for(j=0; j<k-1; j++){
			sum += pow(A[k][j], 2.0);
			}
		sum = A[k][k]-sum;
		if(sum<0){
			printf("Matrix is not positive-definite so Cholesky cannot continue.\n");
			exit(1);
			}
		Lkk = A[k][k] = sqrt(sum);
		for(i=k+1; i<N-1+K; i++){
			sum=0;
			for(j=0; j<k-1; j++) sum += A[i][j]*A[k][j];
			A[i][k] = (A[k][i]-sum) / Lkk;
			}
		}*/
 
	for(i=0; i<N-1+K; ++i){
		for(j=0; j<=i; ++j){
			sum = A[i][j];
			A[j][i] = 0;
			for(k=0; k<j; ++k){
				sum -= A[i][k] * A[j][k];
				}
			if(i == j){
				if(sum < 0){
					printf("Matrix is not positive-definite so Cholesky cannot continue.\n");
					exit(1);
					}
				sum = sqrt(sum);
				if(fabs(sum) < EPS){
					printf("Matrix is not positive-definite so Cholesky cannot continue.\n");
					exit(1);
					}
				A[i][j] = sum;
				}
			else A[i][j] = sum / A[j][j];
			}
		}
 
	for(i=0; i<N-1+K; i++){
		for(j=i+1; j<N-1+K; j++){
			A[i][j] = A[j][i];
			}
		}
	}




/************************************************************************************************************************
 * Same as calculateCholesky but for sparse matrices.									*
 ************************************************************************************************************************/
void calculateCholeskySparse(){
	Ss = cs_schol(1, Cs);
	Ns = cs_chol(Cs, Ss);
	cs_spfree(Cs);
	}




/************************************************************************************************************************
 * An iterative method for solving the A*x=b linear system. The parameter double *M represents the preconditioner	*
 * matrix which in our case is the diagonal of matrix A.								*
 ************************************************************************************************************************/
void calculateCG(){
	double *r, *z, *p, *q, *res;
	int i, iter;
	double rho, rho1, alpha, beta, normB; 

	r = mallocVector();
	z = mallocVector();
	p = mallocVector();
	q = mallocVector();
	res = mallocVector();

	for(i=0; i<N-1+K; i++){
		if(A[i][i]!=0) M[i] = 1/A[i][i];
            	else M[i] = 0;
		}
	
	//res=A*x   
	multMatrixVector(A, x, res);
	
	//r=b-A*x
	for(i=0; i<N-1+K; i++) r[i] = b[i]-res[i];

	iter=0;
	normB=normVector(b);
	if(normB==0.0) normB=1.0;

	while((normVector(r)/normB>optionITOL) && iter<N){
		iter=iter+1;
		
		preconditionerSolve(M, r, z);
		rho=multVectorVector(r, z);

		if(iter==1){
			for(i=0;i<N-1+K;i++) p[i] = z[i];
			}
		else{
			beta=rho/rho1;
			for(i=0; i<N-1+K; i++) p[i] = z[i] + beta * p[i];
			}

		rho1=rho;

		//q=A*p
		multMatrixVector(A, p, q); 
		alpha = rho/multVectorVector(p, q);

		for(i=0;i<N-1+K;i++){
			x[i] = x[i] + alpha * p[i];
			r[i] = r[i] - alpha * q[i];
			}
      		}

	free(r);
	free(z);
	free(p);
	free(q);
	free(res);
	}




/************************************************************************************************************************
 * Same as calculateCG but for sparse matrices.										*
 ************************************************************************************************************************/
void calculateCGSparse(){
	double *r, *z, *p, *q, *res;
	int i, iter;
	double rho, rho1, alpha, beta, normB; 

	r = mallocVector();
	z = mallocVector();
	p = mallocVector();
	q = mallocVector();
	res = mallocVector();

	for(i=0; i<N-1+K; i++){
		int p, found;
		found = -1;
		for(p=Cs->p[i]; p<Cs->p[i+1]; p++){
			if(Cs->i[p] == i) found = p;
			}
		if(found>=0) M[i] = 1/Cs->x[found];
		else M[i] = 0.0;
		}

	//res=A*x
	multMatrixVectorSparse(Cs, x, res);
	
	//r=b-A*x
	for(i=0; i<N-1+K; i++) r[i] = b[i]-res[i];

	iter=0;
	normB=normVector(b);
	if(normB==0.0) normB=1.0;

	while((normVector(r)/normB>optionITOL) && iter<N){
		iter=iter+1;
		
		preconditionerSolve(M, r, z);
		rho=multVectorVector(r, z);

		if(iter==1){
			for(i=0;i<N-1+K;i++) p[i] = z[i];
			}
		else{
			beta=rho/rho1;
			for(i=0; i<N-1+K; i++) p[i] = z[i] + beta * p[i];
			}

		rho1=rho;

		//q=A*p
		multMatrixVectorSparse(Cs, p, q); 
		alpha = rho/multVectorVector(p, q);

		for(i=0;i<N-1+K;i++){
			x[i] = x[i] + alpha * p[i];
			r[i] = r[i] - alpha * q[i];
			}
      		}

	free(r);
	free(z);
	free(p);
	free(q);
	free(res);
	}




/***********************************************************************************************************************
 * An iterative method for solving the A*x=b linear system. The parameter double *M represents the preconditioner	*
 * matrix which in our case is the diagonal of matrix A. Can terminate unexpectedly.					*
 ************************************************************************************************************************/
void calculateBiCG(){
	double *r, *r_tilde, *z, *z_tilde, *p, *p_tilde, *q, *q_tilde, *res;
	int i, iter;
	double rho, rho1, alpha, beta, omega, normB;

	r = mallocVector();
	z = mallocVector();
	p = mallocVector();
	q = mallocVector();
	res = mallocVector();
	r_tilde  = mallocVector();
	z_tilde  = mallocVector();
	p_tilde  = mallocVector();
	q_tilde  = mallocVector();

	for(i=0; i<N-1+K; i++){
		if(A[i][i]!=0) M[i] = 1/A[i][i];
            	else M[i] = 0;
		}
	
	//res=A*x   
	multMatrixVector(A, x, res);

	//r=r_tilde=b-A*x
	for(i=0; i<N-1+K; i++) r[i] = r_tilde[i] = b[i] - res[i];

	iter=0;

	normB=normVector(b);
	if(normB==0.0) normB=1.0;
	
	while((normVector(r)/normB>optionITOL) && iter<N){
		iter=iter+1;
		
		preconditionerSolve(M, r, z);
		preconditionerSolve(M, r_tilde, z_tilde);

		rho=multVectorVector(z, r_tilde);
	
		if(fabs(rho) < EPS){
			printf("Bi-CG algorithm failure due to zero rho.\n");
			exit(1);
			}

		if(iter==1){
			for(i=0;i<N-1+K;i++){
				p[i]=z[i];
				p_tilde[i] = z_tilde[i];
				}
			}
		else{
			beta=rho/rho1;
			for(i=0;i<N-1+K;i++){
				p[i] = z[i] + beta * p[i];
				p_tilde[i] = z_tilde[i] + beta * p_tilde[i];
				}
			}
		rho1=rho;

		//q=A*p
		multMatrixVector(A, p, q); 
		multMatrixVectorTrans(A, p_tilde, q_tilde); 
		omega = multVectorVector(p_tilde, q);

		if(fabs(omega) < EPS){
			printf("Bi-CG algorithm failure due to zero omega.\n");
			exit(1);
			}

		alpha=rho/omega;

		for(i=0;i<N-1+K;i++){
			x[i] = x[i] + alpha * p[i];
			r[i] = r[i] - alpha * q[i];
			r_tilde[i] = r_tilde[i] - alpha * q_tilde[i];
			}
		}

	free(r);
	free(r_tilde);
	free(z);
	free(z_tilde);
	free(p);
	free(p_tilde);
	free(q);
	free(q_tilde);
	free(res);
	}




/***********************************************************************************************************************
 * Same as calculateBiCG but for sparse matrices.										*
 ************************************************************************************************************************/
void calculateBiCGSparse(){
	double *r, *r_tilde, *z, *z_tilde, *p, *p_tilde, *q, *q_tilde, *res;
	int i, iter;
	double rho, rho1, alpha, beta, omega, normB;

	r = mallocVector();
	z = mallocVector();
	p = mallocVector();
	q = mallocVector();
	res = mallocVector();
	r_tilde  = mallocVector();
	z_tilde  = mallocVector();
	p_tilde  = mallocVector();
	q_tilde  = mallocVector();

	for(i=0; i<N-1+K; i++){
		int p, found;
		found = -1;
		for(p=Cs->p[i]; p<Cs->p[i+1]; p++){
			if(Cs->i[p] == i) found = p;
			}
		if(found>=0) M[i] = 1/Cs->x[found];
		else M[i] = 0.0;
		}
	
	//res=A*x   
	multMatrixVectorSparse(Cs, x, res);

	//r=r_tilde=b-A*x
	for(i=0; i<N-1+K; i++) r[i] = r_tilde[i] = b[i] - res[i];

	iter=0;

	normB=normVector(b);
	if(normB==0.0) normB=1.0;
	
	while((normVector(r)/normB>optionITOL) && iter<N){
		iter=iter+1;
		
		preconditionerSolve(M, r, z);
		preconditionerSolve(M, r_tilde, z_tilde);

		rho=multVectorVector(z, r_tilde);
	
		if(fabs(rho) < EPS){
			printf("Bi-CG algorithm failure due to zero rho.\n");
			exit(1);
			}

		if(iter==1){
			for(i=0;i<N-1+K;i++){
				p[i]=z[i];
				p_tilde[i] = z_tilde[i];
				}
			}
		else{
			beta=rho/rho1;
			for(i=0;i<N-1+K;i++){
				p[i] = z[i] + beta * p[i];
				p_tilde[i] = z_tilde[i] + beta * p_tilde[i];
				}
			}
		rho1=rho;

		//q=A*p
		multMatrixVectorSparse(Cs, p, q); 
		multMatrixVectorTransSparse(Cs, p_tilde, q_tilde); 
		omega = multVectorVector(p_tilde, q);

		if(fabs(omega) < EPS){
			printf("Bi-CG algorithm failure due to zero omega.\n");
			exit(1);
			}

		alpha=rho/omega;

		for(i=0;i<N-1+K;i++){
			x[i] = x[i] + alpha * p[i];
			r[i] = r[i] - alpha * q[i];
			r_tilde[i] = r_tilde[i] - alpha * q_tilde[i];
			}
		}

	free(r);
	free(r_tilde);
	free(z);
	free(z_tilde);
	free(p);
	free(p_tilde);
	free(q);
	free(q_tilde);
	free(res);
	}




/************************************************************************************************************************
 * Multiplies the vector M with the vector r and the result is stored at vector z.					*
 ************************************************************************************************************************/
void preconditionerSolve(double *M, double *r, double *z){
	int j;
	for(j=0; j<N-1+K; j++) z[j]=M[j]*r[j];
	}




/************************************************************************************************************************
 * Multiplies the vector r with the vector z and returns the result.							*
 ************************************************************************************************************************/
double multVectorVector(double *r, double *z){
	int i;
	double sum = 0.0;

	for(i=0; i<N-1+K; i++) sum += r[i] * z[i];
	return sum;
	}




/************************************************************************************************************************
 * Multiplies the Matrix A with the vector x and the result is stored at vector res.					*
 ************************************************************************************************************************/
void multMatrixVector(double **A, double *x, double *res){
	int i,j;
	double sum;

	for(i=0;i<N-1+K;i++){
		sum=0.0;
		for(j=0;j<N-1+K;j++) sum += A[i][j] * x[j];
		res[i] = sum;
		}
	}



/************************************************************************************************************************
 * Same as multMatrixVector but for sparse matrices.									*
 ************************************************************************************************************************/
void multMatrixVectorSparse(cs *A, double *x, double *y){
	int j, p;
	for(j=0; j<N-1+K; j++) y[j] = 0.0;
	for(j=0; j<N-1+K; j++){
		for(p=A->p[j]; p<A->p[j+1]; p++){
			y[A->i[p]] = y[A->i[p]] + A->x[p] * x[j];
			}
		}
	}




/************************************************************************************************************************
 * Multiplies the Matrix trans(A) with the vector x and the result is stored at vector res.				*
 ************************************************************************************************************************/
void multMatrixVectorTrans(double **A, double *x, double *res){
	int i,j;
	double sum;

	for(i=0;i<N-1+K;i++){
		sum=0.0;
		for(j=0;j<N-1+K;j++) sum += A[j][i] * x[j];
		res[i] = sum;
		}
	}




/************************************************************************************************************************
 * Same as multMatrixVectorTrans but for sparse matrices.								*
 ************************************************************************************************************************/
void multMatrixVectorTransSparse(cs *A, double *x, double *y){
	int j, p;
	for(j=0; j<N-1+K; j++) y[j] = 0.0;
	for(j=0; j<N-1+K; j++){
		for(p=A->p[j]; p<A->p[j+1]; p++){
			y[j] = y[j] + A->x[p] * x[A->i[p]];
			}
		}
	}




/************************************************************************************************************************
 * Returns the second norm of a vector.											*
 ************************************************************************************************************************/
double normVector(double *v){
	int i;
	double res, sum=0.0;

	for(i=0;i<N-1+K;i++) sum += pow(v[i],2.0);

	return sqrt(sum);   
	}




/************************************************************************************************************************
 * Solves the A*x=b linear system by first solving the L*y=b and then the U*x=y. Does not change the content of matrix	*
 * A so it can be used multiple times for different input vectors b. Puts the results in vector x.			*
 ************************************************************************************************************************/
void solveA(){
	int k, j;
	double *y;

	// allocate memory for vector y
	y = mallocVector();

	// solve the L*y=b system
	y[0] = b[P[0]];
	for(k=1; k<N-1+K; k++){
		double sum;
		sum = 0.0;
		for(j=0; j<k; j++){
			sum += A[k][j]*y[j]; 
			}
		y[k] = b[P[k]] - sum;
		}

	// solve the U*x=y system
	for(k=N-1+K-1; k>=0; k--){
		double sum;
		sum = 0.0;
		for(j=N-1+K-1; j>k; j--){
			sum += A[k][j] * x[j]; 
			}
		x[k] = (y[k] - sum)/A[k][k];
		}

	free(y);
	}





/************************************************************************************************************************
 * Solves the A*x=d linear system by first solving the L*y=b and then the Linv*x=y. Does not change the content of 	*
 * matrix A so it can be used multiple times for different input vectors b. Puts the results in vector x.		*
 ************************************************************************************************************************/
void solveAforChol(){
	int k, j;
	double *y;

	// allocate memory for vector y
	y = mallocVector();

	// solve the L*y=b system
	y[0] = b[0]/A[0][0];
	for(k=1; k<N-1+K; k++){
		double sum;
		sum = 0;
		for(j=0; j<k; j++){
			sum += A[k][j]*y[j]; 
			}
		y[k] = (b[k] - sum)/A[k][k];
		}

	// solve the U*x=y system
	for(k=N-1+K-1; k>=0; k--){
		double sum;
		sum = 0;
		for(j=N-1+K-1; j>k; j--){
			sum += A[k][j] * x[j]; 
			}
		x[k] = (y[k] - sum)/A[k][k];
		}

	free(y);
	}




/*************************************************************************************************************************
 * Same as solveA but for sparse matrices.										*
 *************************************************************************************************************************/
void solveASparse(){
	cs_ipvec(Ns->pinv, b, x, N-1+K);
	cs_lsolve(Ns->L, x);
	cs_usolve(Ns->U, x);
	cs_ipvec(Ss->q, x, b, N-1+K);
	}




/*************************************************************************************************************************
 * Same as solveAforChol but for sparse matrices.									*
 ************************************************************************************************************************/
void solveAforCholSparse(){
	cs_ipvec(Ss->pinv, b, x, N-1+K);
	cs_lsolve(Ns->L, x);
	cs_ltsolve(Ns->L, x);
	cs_pvec(Ss->pinv, x, b, N-1+K);
	}




/************************************************************************************************************************
 * For every .DC statement parsed, it makes the DC sweep and writes the results to the DC sweep output file.		*
 ************************************************************************************************************************/
void writeDCsweepToFile(){
	int i, j, k, pos;	
	struct element_t *elem;
	double val;
	double *res;
	FILE *f;

	if(dotDCctr>0){
		f = fopen(DEFAULT_DCSWEEP_RESULTS_FILENAME, "w");
		if(f==NULL){
			perror(NULL);
			exit(1);
			}

		for(i=0; i<dotDCctr; i++){
			pos=0;

			// find the position in vector b that this voltage source affects
			for(elem=list; elem!=NULL; elem=elem->next){
				if(!strcmp(dotDC[i].element_name, elem->name)) break;
				if(elem->name[0] == 'V') pos++;
				}
			if(pos==K){
				printf("invalid element id. DC sweep failed\n");
				exit(1);
				}

			// recalculate the system for every value in the given range
			for(val=dotDC[i].start_value; val<=dotDC[i].end_value; val+=dotDC[i].increment){
				
				b[N-1+pos] = val;

				if(optionSPARSE){
					if(optionITER){
						if(optionSPD) calculateCGSparse();
						else calculateBiCGSparse();
						res = x;
						}
					else{
						for(j=0; j<N-1+K; j++){
							b[j] = b_dupl[j] ;
							x[j] = 0;
							}
						b[N-1+pos] = val;
						
						if(optionSPD) solveAforCholSparse();
						else solveASparse();
						res = b;
						}
					}
				else{
					if(optionITER){
						if(optionSPD) calculateCG();
						else calculateBiCG();
						}
					else{
						if(optionSPD) solveAforChol();
						else solveA();
						}
					res = x;
					}

				// write only the requested (via a .PLOT statement) information	to the file
				fprintf(f, "%s=%lf", dotDC[i].element_name, val);
				for(j=0; j<dotDC[i].plot_node_ctr; j++){

					for(k=0; k<N; k++){
						if(!strcmp(nodes[k], dotDC[i].plot_node_name[j])){								
							fprintf(f, " V(%s)=%lf", dotDC[i].plot_node_name[j], res[k-1]);
							break;
							}
						}
					if(k==N){
						printf("invalid node name in .PLOT statement.");
						exit(1);
						}
					}
				fprintf(f, "\n");
				}
			}
		fclose(f);
		printf("write DC sweep to file \t[OK]\n");
		}		
	}




/************************************************************************************************************************
 * Writes the results to the results output file.									*
 ************************************************************************************************************************/
void writeDCResultsToFile(){
	int i;
	double *res;
	struct element_t *elem;
	FILE *f;

	f = fopen(DEFAULT_DC_RESULTS_FILENAME, "w");
	if(f==NULL){
		perror(NULL);
		exit(1);
		}

	// write the first N-1 lines that contain the calculated voltage of all the nodes
	for(i=0; i<N-1; i++){
		if(optionSPARSE && !optionITER) res = b;
		else res = x;
		fprintf(f, "%s %lf\n", nodes[i+1], res[i]);
		}
	
	// write the final K lines that contain the current through all the voltage sources	
	for(i=N-1, elem=list; elem!=NULL && i<N-1+K; elem=elem->next){
		if(elem->name[0]=='V'){
			if(optionSPARSE) fprintf(f, "%s %lf\n", elem->name, b[i]);
			else fprintf(f, "%s %lf\n", elem->name, res[i]);
			i++;
			}
		}
	fclose(f);
	printf("write DC results to file \t[OK]\n");
	}




/************************************************************************************************************************
 * Writes a line of results of transient analysis to the results output file.						*
 ************************************************************************************************************************/
void writeTRANResultsToFile(double time){
	int i, j;	
	struct element_t *elem;
	FILE *f;

	f = fopen(DEFAULT_TRAN_RESULTS_FILENAME, "a+");
	if(f==NULL){
		perror(NULL);
		exit(1);
		}

	fprintf(f, "t=%lf", time);

	for(i=0; i<dotTRAN.plot_node_ctr; i++){
		for(j=0; j<N; j++){
			if(!strcmp(nodes[j], dotTRAN.plot_node_name[i])){								
				if(optionSPARSE) fprintf(f, " V(%s)=%lf", nodes[j], b[j-1]);
				else fprintf(f, " V(%s)=%lf", nodes[j], x[j-1]);
				break;
				}
			}
		}

	fprintf(f, "\n");

	fclose(f);
	}




