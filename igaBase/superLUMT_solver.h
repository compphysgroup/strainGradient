
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
#ifndef SUPERLUMT_H_
#define SUPERLUMT_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pdsp_defs.h>

int luSolver(double *aTemp, int *asubTemp, int *xaTemp, double *rhsTemp, int m, int n, int nnz, double * _sol, int _solverInfo)
{
    SuperMatrix   A;
    //NRformat *Astore;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    SuperMatrix   L;       /* factor L */
    SCPformat *Lstore;
    SuperMatrix   U;       /* factor U */
    NCPformat *Ustore;
    SuperMatrix   B;
    int  nrhs, info, b;
    int _temp;
    int  nprocs; /* maximum number of processors to use. */
    int  panel_size, relax, maxsup;
    int  permc_spec;
    trans_t  trans;
    superlu_memusage_t   superlu_memusage;
    nrhs              = 1;
    trans             = NOTRANS;
    nprocs             = 1;
    b                 = 1;
    panel_size        = sp_ienv(1);
    relax             = sp_ienv(2);
    maxsup            = sp_ienv(3);
    //printf("panel_size:%d, relax:%d, maxup:%d, 4Val:%d, 5Val:%d, SL:%u, SU:%u, SRL:%u\n", sp_ienv(1),sp_ienv(2),sp_ienv(3), sp_ienv(4), sp_ienv(5), sp_ienv(6), sp_ienv(7), sp_ienv(8));
    //options
    superlumt_options_t options;
    options.PrintStat = NO;
    if (getenv("OMP_NUM_THREADS")!=NULL){nprocs=atof(getenv("OMP_NUM_THREADS")); printf("num of processors:%u, num of DOF:%u\n", nprocs, m);}
    //creating A matrix
	dCreate_CompCol_Matrix(&A, m, n, nnz, aTemp, asubTemp, xaTemp, SLU_NR, SLU_D, SLU_GE);
    //Astore = A.Store;

    //creating B Matrix (RHS)
    dCreate_Dense_Matrix(&B, m, nrhs, rhsTemp, m, SLU_DN, SLU_D, SLU_GE);

    //creating permutation matrices
    if (!(perm_r = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");
    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering 
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */    	
    permc_spec=1; get_perm_c(permc_spec, &A, perm_c);
    pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);
    
    if (info==0) {
		/* Access the solution matrix. */
		double *sol = (double*) ((DNformat*) B.Store)->nzval;
		for(_temp=0; _temp<m; ++_temp){_sol[_temp] =sol[_temp];}
		if(_solverInfo){
			Lstore = (SCPformat *) L.Store;
			Ustore = (NCPformat *) U.Store;
			printf("#NZ in factor L = %d\n", Lstore->nnz);
			printf("#NZ in factor U = %d\n", Ustore->nnz);
			printf("#NZ in L+U = %d\n", Lstore->nnz + Ustore->nnz - L.ncol);
		}
   }
    else {
    	/*if solver fails*/
    	printf("pdgssv() error returns INFO= %d\n", info);
    	if ( info <= n ) { /* factorization completes */
    		superlu_dQuerySpace(nprocs, &L, &U, panel_size, &superlu_memusage);
    		printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n", superlu_memusage.for_lu/1024/1024, superlu_memusage.total_needed/1024/1024, superlu_memusage.expansions);
        }
    	exit(-1);
    }
    //freeing memory
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);
    return info;
}
#endif