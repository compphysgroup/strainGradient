
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#ifndef SUPERLU_H_
#define SUPERLU_H_
#include "slu_ddefs.h"

int luSolver(double *aTemp, int *asubTemp, int *xaTemp, double *rhsTemp, int m, int n, int nnz, double * _sol, int _solverInfo)
{
	/*linear solver variables devlaration*/
	SuperMatrix A, B;
	//NRformat *Astore;
	int      *perm_c; /* column permutation vector */
	int      *perm_r; /* row permutations from partial pivoting */
	SuperMatrix L;      /* factor L */
	SuperMatrix U;      /* factor U */
	SCformat *Lstore; NCformat *Ustore;
	int      nrhs, info;
	int 	_temp;
	mem_usage_t   mem_usage;
	superlu_options_t options;
	SuperLUStat_t stat;
	/*additional solver variables to enable refinement*/
	char           equed[1];
	yes_no_t       equil;
	SuperMatrix    X;
	int            *etree;
	void           *work=0;
	int            lwork;
	int            i;
	double         *rhsx;
	double         *R, *C;
	double         *ferr, *berr;
	double         u, rpg, rcond;
		
	/* Set the default input options:
	options.Fact = DOFACT;
	options.Equil = YES;
	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
	options.Trans = NOTRANS;
	options.IterRefine = NOREFINE;
	options.SymmetricMode = NO;
	options.PivotGrowth = NO;
	options.ConditionNumber = NO;
	options.PrintStat = YES;
	 */
	set_default_options(&options);
	nrhs   = 1;
	/*additional statments to enable refinement*/
	/* Defaults */
	lwork = 0;
	equil = YES;	
	u     = 1.0;
	options.ConditionNumber = YES;
	options.Equil = equil;
	options.DiagPivotThresh = u;
	/* Add more functionalities that the defaults. */
	options.PivotGrowth = YES;    /* Compute reciprocal pivot growth */
	options.ConditionNumber = YES;/* Compute reciprocal condition number */
	options.IterRefine = DOUBLE;  /* Perform double-precision refinement */
    
	/*Initializing A matrix(m x n dimensions)*/
	dCreate_CompCol_Matrix(&A, m, n, nnz, aTemp, asubTemp, xaTemp, SLU_NR, SLU_D, SLU_GE);
	//Astore = A.Store;
	//printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
	
	 
	/*Initializing B matrix(m x 1 dimensions)*/
    dCreate_Dense_Matrix(&B, m, nrhs, rhsTemp, m, SLU_DN, SLU_D, SLU_GE);
    
    /*Initializing column and row permutation vectors*/
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    
    /*additional statments to enable refinement*/
	if ( lwork > 0 ) {work = SUPERLU_MALLOC(lwork); if ( !work ) {ABORT("DLINSOLX: cannot allocate work[]");}}
	if ( !(rhsx = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsx[].");
	dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
	if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
	if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) ) ABORT("SUPERLU_MALLOC fails for R[].");
	if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) ) ABORT("SUPERLU_MALLOC fails for C[].");
	if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) ) ABORT("SUPERLU_MALLOC fails for ferr[].");
	if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) ) ABORT("SUPERLU_MALLOC fails for berr[].");
			
    /*Initialize the statistics variables. */
    StatInit(&stat);
    
    /*Solving system of equations AX=B*/
    dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr, &mem_usage, &stat, &info);
    /*Postprocessing: Extracting solver info*/
    if ( info == 0  || info == n+1 ) {
		/* This is how you could access the solution matrix. */
		double *sol = (double*) ((DNformat*) X.Store)->nzval;
		for(_temp=0; _temp<m; ++_temp){_sol[_temp] =sol[_temp];}


		/* Printout solver information.*/
		if(_solverInfo){
			printf("Recip. condition number = %e\n", rcond);
			if ( options.PivotGrowth == YES ) printf("Recip. pivot growth = %e\n", rpg);
			//if ( options.ConditionNumber == YES ) printf("Recip. condition number = %e\n", rcond);
			if ( options.IterRefine != NOREFINE ) {
				printf("Iterative Refinement:\n");
				printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
				for (i = 0; i < nrhs; ++i) printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i]);
			}
			Lstore = (SCformat *) L.Store;
			Ustore = (NCformat *) U.Store;
			//printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
			//printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
			//printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
			//printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);
			//printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
			if ( options.PrintStat ) StatPrint(&stat);
		}
	}
    else {
    	/*if solver fails*/
    	printf("Recip. condition number = %e\n", rcond);
    	printf("dgssvx() error returns INFO= %d\n", info);
    	if ( info > 0 && lwork == -1 ) printf("** Estimated memory: %d bytes\n", info - n);
    	if ( info <= n ) { /* factorization completes */
    		dQuerySpace(&L, &U, &mem_usage);
    		printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
    	}
    	if ( options.PrintStat ) StatPrint(&stat);
    }
    fflush(stdout);

    /* De-allocate storage */
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);
	SUPERLU_FREE (rhsx);
	SUPERLU_FREE (etree);
	SUPERLU_FREE (R);
	SUPERLU_FREE (C);
	SUPERLU_FREE (ferr);
	SUPERLU_FREE (berr);
	Destroy_SuperMatrix_Store(&X);
    
    return info;
}
#endif
