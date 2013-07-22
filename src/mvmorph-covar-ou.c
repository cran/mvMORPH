/*-----------Matrice de covariance pour un processus Ornstein-Uhlenbeck-------------------------*/
/*-mvMORPH 1.0.2 - 2014 - Julien Clavel - julien.clavel@hotmail.fr/julien.clavel@univ-lyon1.fr--*/
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

static void mvmorph_covar_OU(int *nt,
				   double *A, 
				   double *ans, 
			       double *alpha,
				   double *sigma) {
  double sij, ti, tj, T, temp, dat;
  int n = *nt;
  int i, j;

/* Compute the vcv-Ou matrix */
for(i=0; i<n; i++){
	for(j=0; j<=i; j++){
/* allows non-ultrametric trees */
      sij=A[j*n+i];
	  ti = A[i+i*n]-A[j+i*n];
	  tj = A[j+j*n]-sij;
	  T=ti+tj;
	 
	temp=(1-exp(-2.0*alpha[0]*sij))*exp(-1.0*alpha[0]*T)/(2.0*alpha[0]);  
    dat=temp*sigma[0];	
	ans[i*n+j]=dat;
	if (j != i){
	ans[j*n+i]=dat;
	}
  }
}

}

SEXP mvmorph_covar_ou(SEXP A, SEXP alpha, SEXP sigma) {
 int nt; 
	PROTECT(coerceVector(A,REALSXP)); 
	nt=INTEGER(GET_DIM(A))[0];
	SEXP ans;
    PROTECT(ans = allocMatrix(REALSXP,nt,nt)); 
	mvmorph_covar_OU(&nt,REAL(A),REAL(ans),REAL(alpha), REAL(sigma));
  UNPROTECT(2);
  return ans;
}


