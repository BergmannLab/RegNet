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

#ifndef __EDGE_H__
#define __EDGE_H__

#include <string.h>
#include <stdlib.h>
#include "Vector.h"
#include <vector>
class Node;

#define WITH_COACTIVATION

enum EdgeT {NIL, TRANSCRIPTION, DEGRADATION, SEQUESTRATION,HILL,AGGREGATE,DISAGGREGATE,WELL,COACTIVATION};

typedef std::vector<float> ParamList;


class Edge{
 public:
  Node *pto; //the nodes it points to
  Node *pfrom;// edge origin
  Node *paux; //auxillary node
  int nb; // the index of incoming edge on target node
  int param_idx; /*index of parameter values in big parameter table*/
  int npar;			/* number of parameters */
  EdgeT type;
  ParamList* all_parameters;      /**< vector containing all the network parameters */


  static float default_param;

  static float logmin;// for uniform param sampling
  static float logmax;// for uniform param sampling
  static float degr_mean_log;
  static float degr_var_log;

  Edge(){pfrom=pto=paux=NULL;nb=-1;npar=1;}
  Edge(Node *n,int i){pto=n;nb=i;npar=1;pfrom=paux=NULL;}
  Edge(Node *to,Node *from,EdgeT t){pto=to;pfrom=from;type =t;paux=NULL;}
  Edge(Node *to,Node *from,int i=-1){pto=to;nb=i;pfrom=from;npar=1;paux=NULL;}
  int setValue(float v){return setValue(0,v);}
  int setValue(int i,float v);
  int setValue();

  int setValues(const float *p);
  float getValue(int i=0)const
  {float f=0;if(i<npar)f=(*all_parameters)[param_idx+i];return f;}
  int getValues(float *p)const;
  EdgeT getType()const{return type;}
  float compute(float y,const Node *to)const;
  Node *origin()const{return pfrom;};
  Node *target()const{return pto;};
  Node *auxil()const{return paux;};
  void setAux(Node *aux){paux=aux;}
  void setNumber(int n){nb=n;}
  void setParameterIndex(int n){param_idx=n;}
  void initParameters(ParamList& params);
  float getHillDerivative()const;
  void initJacobian(float y,Vector& vec,const Node *to)const;
  void updateJacobianRow(float *row,const std::vector<int>& nmap,const Node *to)const;
  int getNbParameters(){return npar;};
  int randomParameters()
  {for(int i=0;i<npar;i++){setValue(i,randParam());}return npar;}
  int normalRandomParameters()
  {for(int i=0;i<npar;i++){setValue(i,normalRandParam());}return npar;}

  void print(std::ostream& out);
  int hasParamIdx(int i)const
  {if(param_idx<=i && param_idx+npar>i){return i-param_idx+1;}return 0;}

  static float randParam();
  static float normalRandParam();
  static EdgeT edgeType(char *s);
  static void arrowType(EdgeT e,char *s);
};

inline float Edge::randParam(){
  return exp(logmin+RND(logmax-logmin));
}
inline float Edge::normalRandParam(){
  return exp(norm_sample(Edge::degr_mean_log, Edge::degr_var_log));
}

inline EdgeT Edge::edgeType(char *s){
#ifndef HILL_ONLY
  if(!strcmp(s,"->")){return TRANSCRIPTION;}
#else
  if(!strcmp(s,"->")){return HILL;}
#endif
  if(!strcmp(s,"-|")){return DEGRADATION;}
  if(!strcmp(s,"|-|")){return SEQUESTRATION;}
  if(!strcmp(s,"->>")){return HILL;}
  if(!strcmp(s,">>")){return AGGREGATE;}
  if(!strcmp(s,"<<")){return DISAGGREGATE;}
  if(!strcmp(s,"=|")){return WELL;}
  if(!strcmp(s,"=>")){return COACTIVATION;}
  return NIL;
}

inline void Edge::arrowType(EdgeT e, char *s){
  switch(e){
  case TRANSCRIPTION:sprintf(s,"->");break;
  case DEGRADATION:sprintf(s,"-|");break;
  case SEQUESTRATION:sprintf(s,"|-|");break;
  case HILL:sprintf(s,"->>");break;
  case AGGREGATE:sprintf(s,">>");break;
  case DISAGGREGATE:sprintf(s,"<<");break;
  case WELL:sprintf(s,"=|");break;
  case COACTIVATION:sprintf(s,"=>");break;
  default:sprintf(s,"-");break;
  }
}




typedef Edge EdgePt;

#endif
