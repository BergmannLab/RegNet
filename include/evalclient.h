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

#ifndef __EVALCLIENT_H__
#define __EVALCLIENT_H__

#include "Vector.h" 

class EvalClient{
 protected:
  int sock[2];
 public:
  int init(char *sockname);
  void sendOK();
  void sendQuery(int n, float& v,Vector *par=NULL);
  void closeClient();
  void closeServer(char *sockname);
  float eval(char *sockname,Vector *par);
  void getParameterName(char *sockname,int i, char *name);
  void sendNameQuery(int i, char *name);
};

#endif
