/** 
Copyright (c) 2013 Micha Hersch, University of Lausanne, Swiss Institute of Bioinformatics

This file is part of RegNet

RegNet is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#ifdef CVODE_INTEGRATION

#include "rnet.h"
#include <stdio.h>
#include "cvode/llnltyps.h"
#include "cvode/cvode.h" 
#include "cvode/nvector.h" 
#include "cvode/dense.h"
#include "cvode/cvdense.h"

extern "C" {
  N_Vector N_VNew(integer n, void *machEnv);
  void *CVodeMalloc(integer N, RhsFn f, real t0, N_Vector y0, int lmm, int iter,
		    int itol, real *reltol, void *abstol, void *f_data,
		    FILE *errfp, boole optIn, long int iopt[], real ropt[],
		    void *machEnv);	    
  void CVDense(void *cvode_mem, CVDenseJacFn djac, void *jac_data);
  int CVode(void *cvode_mem, real tout, N_Vector yout, real *t, int itask);
  void  N_VFree(N_Vector x);
  void CVodeFree(void *cvode_mem);
}

#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


static void f_cvode(integer N, real t, N_Vector y, N_Vector ydot, void *f_data){
  RegNet *net = (RegNet *) f_data;
  net->updateStates(N_VDATA(y));
  net->getUnconstrainedDerivatives(N_VDATA(ydot));
}

static void Jac(integer N, DenseMat J, RhsFn f, void *f_data, real t,
                N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3){
 RegNet *net = (RegNet *) jac_data;
  net->updateStates(N_VDATA(y));
  net->updateJacobian();
  const Matrix& jac = net->getJacobian();
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      DENSE_ELEM(J,i,j) = jac.At(i,j);
    }
  }
}

static void PrintFinalStats(long int iopt[]){
  printf("\nFinal Statistics.. \n\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nje = %ld\n",
	 iopt[NST], iopt[NFE], iopt[NSETUPS], iopt[DENSE_NJE]);
  printf("nni = %-6ld ncfn = %-6ld netf = %ld\n \n",
	 iopt[NNI], iopt[NCFN], iopt[NETF]);
}

int RegNet::findSteadyState(int maxsteps, float tol){
  float reltol = 1e-2;
  float tol2 = tol*tol;
  long int iopt[OPT_SIZE];
  real ropt[OPT_SIZE];
  N_Vector y, y_old, y_diff,abstol;
  void *cvode_mem;
  int ret = 0;
   resizeNewton();
  int n=nodes.size()-input.size(); // problem size
  y= N_VNew(n, NULL); 
  y_old= N_VNew(n, NULL); 
  y_diff= N_VNew(n, NULL); 
  abstol= N_VNew(n, NULL); 
   N_VConst(1,y);
  N_VConst(1e-2,abstol);
  for(int i=0;i<OPT_SIZE;i++){iopt[i]=0;ropt[i]=0;}
  iopt[MXSTEP] = 30000; //max nb of integration steps
  cvode_mem = CVodeMalloc(n, f_cvode,0,y , BDF, 
			  NEWTON, SV, &reltol, 
			  abstol, (void *) this, stderr, TRUE, 
			  iopt, ropt,NULL);

 if (cvode_mem == NULL) { printf("CVodeMalloc failed.\n"); return 0; }
 
 CVDense(cvode_mem, Jac, this);
 
 float tout,t;
 tout=1;
 for (int iout=0; iout < 10; iout++) {
   N_VAddConst(y,0,y_old);
   int flag = CVode(cvode_mem, tout, y, &t, NORMAL);
   //  printf("%f %f\n",tout,N_VIth(y,nmap(getOutputIndex())));
   if (flag != SUCCESS) {
     //       printf("CVode failed, flag=%d.\n", flag); 
      return 0;
   }
   N_VLinearSum(1,y,-1,y_old,y_diff);//y_diff = y-y_old
   if(N_VDotProd(y_diff,y_diff)<tol2){ret=1;break;}
   tout *=1.5;   
 }

 updateStates(N_VDATA(y));
 N_VFree(y);                  /* Free the y and abstol vectors */
 N_VFree(abstol);   
 CVodeFree(cvode_mem);        /* Free the CVODE problem memory */
 //PrintFinalStats(iopt);       /* Print some final statistics   */

return 1;
}



#endif
