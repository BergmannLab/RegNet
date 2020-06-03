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

#include "matrix.h"
#include "mex.h"
#include "evalclient.h"

void mexFunction(int nlhs,mxArray* plhs[],int nrhs, const mxArray *prhs[]){

  char *sockname = mxArrayToString(prhs[0]);
  EvalClient ec;
  int dim = ec.init(sockname);
  ec.sendOK();
  ec.closeClient();
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  double *output = mxGetPr(plhs[0]); 
  output[0]= dim;
}
