/*-Matrice de covariance pour un processus Ornstein-Uhlenbeck multivarie--------------*/
/*-mvMORPH 1.0.2 - 2014 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/
#include "ouch.h"

static void mvmorph_covar_mat_nult (int *nchar, 
                   int *nt,
				   double *bt, 
			       double *lambda, 
			       double *S, 
			       double *sigmasq, 
			       double *V) {
  double *U, *G, *F, *exp1, *exp1l, *exp2, *exp2l;
  double sij, ti, tj, tmp;
  int n = *nchar, nn = *nt;
  int i, j, k, l, s, r;
  U = Calloc(n*n,double);
  G = Calloc(n*n,double);
  F = Calloc(n*n,double);
  exp1 = Calloc(n*n,double);
  exp2 = Calloc(n*n,double);
  exp1l = Calloc(n*n,double);
  exp2l = Calloc(n*n,double);


/* Calcul de la fonction de correlation BBt*/
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      U[i+j*n] = 0;
      for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  U[i+j*n] += S[k+i*n]*sigmasq[k+l*n]*S[l+j*n];
	}
      }
    }
  }

/* Calcul de la matrice de covariance */
  for (i = 0; i < nn; i++) {
    for (j = 0; j <= i; j++) {
	/* si l'arbre n'est pas ultrametrique */
      sij = bt[j+i*nn];
      ti = bt[i+i*nn]-sij;
      tj = bt[j+j*nn]-bt[i+j*nn];

/* preparation du calcul des matrices exponentielles */ 
  for (k = 0; k < n; k++) {
	exp1l[k+k*n] = exp(-lambda[k]*ti);
	exp2l[k+k*n] = exp(-lambda[k]*tj);
	for(l = 0; l < n; l++) {		
		if((k+l*n) != (k+k*n)){
        exp1l[k+l*n] = 0;
		exp2l[k+l*n] = 0;
			}
		}
      }
/* Calcul de G - Gardiner 2004 4.4.47 - Bartoszek et al. 2012 - Najfeld & Havel 1995 */    
 for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	// voir aussi Najfeld & Havel 1995
	  if((lambda[k]+lambda[l]) != 0.0){
	  G[k+l*n] = ((1.0-exp(-1.0*(lambda[k]+lambda[l])*sij))/(lambda[k]+lambda[l]))*U[k+l*n];
	  }else{
	  G[k+l*n] = U[k+l*n]*sij ;
	  }
	}
  }
/* Calcul de SGS - Gardiner 2004 - calcul des matrices exponentielles*/
     for (k = 0; k < n; k++) {
		for (l = 0; l < n; l++) {
		       F[k+l*n] = 0;
			exp1[k+l*n] = 0;
		    exp2[k+l*n] = 0;
			  V[i+nn*(k+n*(j+nn*l))] = 0; 
	          V[j+nn*(k+n*(i+nn*l))] = 0; 
	  for (r = 0; r < n; r++) {
	    for (s = 0; s < n; s++) {
	         F[k+l*n] += S[k+r*n]*G[r+s*n]*S[l+s*n];
		  exp1[k+l*n] += S[k+r*n]*exp1l[r+s*n]*S[l+s*n];
		  exp2[k+l*n] += S[k+r*n]*exp2l[r+s*n]*S[l+s*n];
	    }
	  }
	}
  }
/* Integrale + matrices exponentielles */  
  for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  for (r = 0; r < n; r++) {
	    for (s = 0; s < n; s++) {
	    tmp = exp1[k+r*n]*F[r+s*n]*exp2[l+s*n];
	    V[i+nn*(k+n*(j+nn*l))] += tmp;
	      if (j != i) 
		V[j+nn*(l+n*(i+nn*k))] += tmp;
	    }
	  }
	}
  }
/* Fin de la boucle pour la matrice de covariance
on libère la mémoire */	  
}
}	  
  Free(U);
  Free(G);
  Free(F);
  Free(exp1);
  Free(exp2);
  Free(exp1l);
  Free(exp2l);
}

SEXP mvmorph_covar_mat (SEXP nterm, SEXP bt,SEXP lambda, SEXP S, SEXP sigmasq) {
  int nprotect = 0;
  SEXP V;
  int nchar, nt, vdim[2];
  nt = INTEGER(nterm)[0];
  nchar = GET_LENGTH(lambda);
  vdim[0] = nt*nchar; vdim[1] = vdim[0];
  PROTECT(V = makearray(2,vdim)); nprotect++;
  mvmorph_covar_mat_nult(&nchar,&nt,REAL(bt),REAL(lambda),REAL(S),REAL(sigmasq),REAL(V));
  UNPROTECT(nprotect);
  return V;
}
