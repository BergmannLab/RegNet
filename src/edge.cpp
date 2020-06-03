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

#include "edge.h"
#include "node.h"
#include <assert.h>

using namespace std;

float Edge::default_param = -1;
#ifdef MENDOZA
float Edge::logmin = -3;
float Edge::logmax = 2;
#else
float Edge::logmin = -4;
float Edge::logmax = 1;
#endif
float Edge::degr_mean_log = -2;//-3
float Edge::degr_var_log = 0.5;//2



float Edge::compute(float y,const Node *to)const{
  float dy=0;
  float x =0;
  int si;
  switch(type){
  case TRANSCRIPTION: //dy = k*x
    dy = getValue()*pfrom->getState();
    break;
  case COACTIVATION: // dy = -k*x*z
    dy = getValue()*pfrom->getState()*paux->getState();
    break;
  case DEGRADATION: 
  case SEQUESTRATION://dy = -k*x*y  
    dy = - getValue()*y*pfrom->getState();
    break;
  case HILL: //dy =k1*(x*x)/(k2+x*x)
    x = pfrom->getState();
    x= x*x;//HILL_COEFF=2
    dy = getValue(0)*x/(getValue(1)+x);
    break;
  case AGGREGATE:
    si = (pto==to)?1:-1;
    dy = si*getValue()*pfrom->getState()*
      paux->getState();
    break;
  case DISAGGREGATE:
    si = (pfrom==to)?-1:1;
    dy = si*getValue()*pfrom->getState();
    break;
  case WELL:
     dy = -getValue()*pfrom->getState()*paux->getState()*y;//dy = -k*x*z*y
       break;
#ifdef WITH_GROWTH 
  case GROWTH:
    
   break;
#endif
  default:;
  }
  return dy;

}


void Edge::initJacobian(float y,Vector& vec,const Node *to)const{
  int si;
  switch(type){
  case TRANSCRIPTION:
    vec[pfrom->getIndex()]+= getValue();break;
  case COACTIVATION:
    vec[pfrom->getIndex()]+= getValue()*paux->getState();
    vec[paux->getIndex()]+= getValue()*pfrom->getState();
    break;
  case DEGRADATION:
  case SEQUESTRATION:
    vec[pfrom->getIndex()]-= getValue()*y;
    vec[pto->getIndex()]-= getValue()*pfrom->getState();
    break;
  case HILL:
    vec[pfrom->getIndex()] += getHillDerivative();
    break;
  case AGGREGATE:
    si = (pto==to)?1:-1;
    vec[pfrom->getIndex()] += si*getValue()*
      paux->getState();
    vec[paux->getIndex()] += si*getValue()*
      pfrom->getState();
    break;
  case DISAGGREGATE:
    si = (pfrom==to)?-1:1;
    vec[pfrom->getIndex()] += si*getValue();
    break;
  case WELL:
    vec[pfrom->getIndex()]-= getValue()*y*paux->getState();
    vec[pto->getIndex()]-= getValue()*pfrom->getState()*paux->getState();
    vec[paux->getIndex()]-= getValue()*y*pfrom->getState(); 
    break;
 default:;
  }
}

void Edge::updateJacobianRow(float *row,const vector<int>& nmap, const Node *to)const{
  if(!pfrom->isConstrained()){
    int si;//sign
    float y = pto->getState();
    switch(type){
    case TRANSCRIPTION:
      row[nmap[pfrom->getIndex()]]+= getValue();break;
    case COACTIVATION:
      row[nmap[pfrom->getIndex()]]+= getValue()*paux->getState();
      if(!paux->isConstrained()){
	row[nmap[paux->getIndex()]] += getValue()*pfrom->getState();
      }
      break;   
    case DEGRADATION:
    case SEQUESTRATION:
      row[nmap[pfrom->getIndex()]]-= getValue()*y;
      row[nmap[pto->getIndex()]]-= getValue()*pfrom->getState();
      break;
    case HILL:
      row[nmap[pfrom->getIndex()]] += getHillDerivative();
      break;
    case AGGREGATE:
      si = (pto==to)?1:-1;
      row[nmap[pfrom->getIndex()]] += si*getValue()*paux->getState();
      row[nmap[paux->getIndex()]] += si*getValue()*pfrom->getState();
      break;
    case DISAGGREGATE:
      si = (pfrom==to)?-1:1;
      row[nmap[pfrom->getIndex()]] += si*getValue();
    break;
    case WELL:
      row[nmap[pfrom->getIndex()]]-= getValue()*y*paux->getState();
      row[nmap[pto->getIndex()]]-= getValue()*pfrom->getState()*paux->getState();
      row[nmap[paux->getIndex()]]-= getValue()*y*pfrom->getState();    
  break;
    default:;
    }
  }
#ifdef WITH_COACTIVATION
  else if(paux!=NULL){
    if(!paux->isConstrained()){
      switch(type){//TODO: add also other types
      case COACTIVATION:
	row[nmap[paux->getIndex()]] += getValue()*pfrom->getState();
	break;
      default:;
      }
    }
  }
#endif
}


float Edge::getHillDerivative()const{
  float a = getValue(0);
  float kh = getValue(1);
  assert(kh>0);
  float x =  pfrom->getState();
  return a*kh*2*x/((kh+x*x)*(kh+x*x));//HILL_COEFF =2
}

int Edge::setValue(int i,float v){
  if(i<npar){
    (*all_parameters)[param_idx+i]=v;
    return 1;
  }
  else{
    return 0;
  }
}



int Edge::setValues(const float *p){
  memcpy(&((*all_parameters)[param_idx]),p,npar*sizeof(float));
  return npar;
}

int Edge::getValues(float *p)const{
  memcpy(p,&((*all_parameters)[param_idx]),npar*sizeof(float));
  return npar;
}

void Edge::initParameters(ParamList& params){
  all_parameters = &params;
  setParameterIndex(params.size());
  switch(type){
  case TRANSCRIPTION: 
  case COACTIVATION:
  case DEGRADATION:
  case SEQUESTRATION:
  case AGGREGATE:
  case DISAGGREGATE:
  case WELL:
    npar = 1;
    params.push_back(Edge::default_param);
    break;
  case HILL:
    npar = 2;
    params.push_back(Edge::default_param);
    params.push_back(Edge::default_param);
    break;
  default:;
  }
}

void Edge::print(ostream& out){
  char arr[10];
  arrowType(type,arr);
  out<<pfrom->getName();
  //  cout<<auxil()<<endl;
  if (auxil()){
    out<<","<<auxil()->getName();
  }
  out<<" "<<arr<<" "<<pto->getName();
}

