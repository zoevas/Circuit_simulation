#ifndef MNA_H
#define MNA_H

#include "csparse.h"

#define ORDER_NODES

#define MNA_DC_ANALYSIS		0
#define MNA_TRANSIENT_ANALYSIS	1

#define DEFAULT_DC_RESULTS_FILENAME		"DC_RESULTS.TR"
#define DEFAULT_TRAN_RESULTS_FILENAME		"TRAN_RESULTS.TR"
#define DEFAULT_DCSWEEP_RESULTS_FILENAME	"DC_SWEEP_RESULTS.TR"

#define EPS (1e-14)

short same_N;	// turns true if DC and TRANS analysis have calculated the same N


double **A;	// the matrix A of linear system A*x=b
double *x;	// the solution vector
double *b;	// the vector b of linear system A*x=b
double *b_dupl;	// a copy of vector b used in sparse LU and sparse Cholesky
double **hC;	// the matrix C multiplied by constand 1/h or 2/h
 
double *dc_result;
double *e;	// the vector e
double *e_prev;

int *P;		// the pivoting vector P
int N; 		// the count of the nodes used by the Modified Nodal Analysis
int K;		// the count of the voltage sources (including the inductors converted to zero-value voltage sources)
char **nodes;	// an array that keeps pointers to all the node IDs.
double *M;	// the preconditioner matrix
int nz;		// nonzeros count
cs *Cs;		// the compressed sparse matrix A
css *Ss;	// used in sparse methods
csn *Ns;	// used in sparse methods

int nzC;	// nonzeros count of matrix C in transient analysis
int nzG;	// nonzeros count of matrix G in transient analysis
cs *tran_Csc;	// used in sparse methods
cs *tran_Gsc;	// used in sparse methods

double **C;	//array for the capacitors for the transient analysis
double **G;	//array for the transient analysis

void nodeMapping(short dc_or_trans);
void initLinearSystem();
void initLinearSystemSparse();
void initLinearSystemTrans();
void initLinearSystemTransSparse();
void freeLinearSystem();
void freeLinearSystemSparse();
void freeLinearSystemTrans();
void freeLinearSystemTransSparse();
void calculateAb();

void calculateGC();
void calculateBE_A(double h);
void calculateBE_b(double time, double final_time);
void calculateTR_A(double h);
void calculateTR_b(double time, double final_time);

void calculateGCSparse();
void calculateBE_ASparse(double h);
void calculateBE_bSparse(double time, double final_time);
void calculateTR_ASparse(double h);
void calculateTR_bSparse(double time, double final_time);

void calculateAbSparse();
void calculateLU();
void calculateLUSparse();
void calculateCholesky();
void calculateCholeskySparse();
void calculateCG();
void calculateCGSparse();
void calculateBiCG();
void calculateBiCGSparse();

void preconditionerSolve(double *M, double *r, double *z);
double multVectorVector(double *r, double *z);
void multMatrixVector(double **A, double *x, double *res);
void multMatrixVectorSparse(cs *As, double *x, double *y);
void multMatrixVectorTrans(double **A, double *x, double *res);
void multMatrixVectorTransSparse(cs *As, double *x, double *res);
double normVector(double *v);

double EXP_Comp(double i1, double i2, double td1, double tc1, double td2, double tc2, double time, double final_time);
double SIN_Comp(double i1, double ia, double fr, double td, double df, double ph, double time, double final_time);
double PWL_Comp(double *t, double *y, int n, double time, double final_time);
double PULSE_Comp(double i1, double i2, double td, double tr, double tf, double pw, double per, double time, double final_time);

void addMatrices(double  **A, double **B, double **res);
void multMatrixConstant(double **A, double constant, double **res);
void addVectors(double *x , double *y, double *res);
void assignVectorToVector(double *x, double *y);
double linearInterpolation(double x0, double y0, double x1, double y1, double x);

void solveA();
void solveAforChol();
void solveASparse();
void solveAforCholSparse();
void writeDCsweepToFile();
void writeDCResultsToFile();
void writeTRANResultsToFile(double time);

void setToZero(double *vector);
double *mallocVector();
double **mallocMatrix();

void printMatrix(double **matrix);
void printVector(double *vector);

#endif


