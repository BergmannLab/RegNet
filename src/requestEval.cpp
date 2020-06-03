/** 
    
    Copyright (C) 2013 Micha Hersch, University of Lausanne
    email: micha.hersch@unil.ch

    This program is free software: you can redistribute it and/or modify
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


#include "matrix.h"
#include "mex.h"
#include "evalclient.h"

void mexFunction(int nlhs,mxArray* plhs[],int nrhs, const mxArray *prhs[]){

  const mxArray *mat_in = prhs[0];
  double *p = mxGetPr(mat_in);
  Vector par;
  int n = mxGetM(mat_in)*mxGetN(mat_in);
  float v = 0;
  par.Resize(n);
  for(int i=0;i<n;i++){
    par[i]= p[i];
  }
  char *sockname = mxArrayToString(prhs[1]);
  EvalClient ec;
  v= ec.eval(sockname,&par);
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  double *output = mxGetPr(plhs[0]); 
  output[0]= v;
}
