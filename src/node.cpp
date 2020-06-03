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


#include "node.h"
#include <stdlib.h>
#include <time.h>


using namespace std;




float Node::logmin = -4;
float Node::logmax = 1;

float Node::trans_mean_log = -2;//-3
float Node::trans_var_log = 0.5;//2

float Node::src_mean_log = -3;
float Node::src_var_log = 0.5;//1

float Node::degr_mean_log = -2;//-3
float Node::degr_var_log = 0.5;//2


int Node::nb_instances = 0;


#ifndef DISCRETE

#ifndef MENDOZA
#ifndef SIGMOID

float Node::f(float y)const{
  float dy;
  EdgePtList::const_iterator eit;
  if(y<0){y=state;}
  dy = (src()-degrad()*y);//source and degradation
  for(eit=edge_in.begin();eit<edge_in.end();eit++){
    dy += (*eit)->compute(y,this);
  }
  return dy;
}
void Node::updateJacobianRow(float *row, vector<int>& nmap)const{
  EdgePtList::const_iterator eit;
  for(eit=edge_in.begin();eit<edge_in.end();eit++){
    Node *parent = (*eit)->origin();
    if(!parent->isConstrained()){
      row[nmap[parent->getIndex()]]=0;
    }
  }
  row[nmap[getIndex()]]= -degrad();//degradation term
  for(eit=edge_in.begin();eit<edge_in.end();eit++){
    (*eit)->updateJacobianRow(row,nmap,this);
  }    
}


#endif //SIGMOID
#endif //MENDOZA


void Node::initJacobianRow(Vector& vec)const{
  vec.Zero();
  if(!constrained){
    if(!no_degrad){vec[getIndex()] = -degrad();}
    EdgePtList::const_iterator eit;
    float y = state;
    for(eit=edge_in.begin();eit<edge_in.end();eit++){
      (*eit)->initJacobian(y,vec,this);
    }
  }
}

#endif //DISCRETE


int Node::initParameters(ParamList& params){
  int nb=0;
  all_parameters = &params;
  src_idx = params.size();
  params.push_back(randParam());
  nb++;
#ifndef OPTIMIZE_NODE_PARAM
  degrad_idx = params.size();
  params.push_back(randParam());
  nb++;
#endif
  return nb;
}


void Node::randomParams(){
  EdgePtList::const_iterator eit;
  setSource(randParam());
  setDegrad(randParam());
  for(eit=edge_in.begin();eit<edge_in.end();eit++){
    (*eit)->randomParameters();
  }
}

void Node::randomParamsNormal(){
  EdgePtList::const_iterator eit;
  if(!no_source)
    {setSource(exp(norm_sample(Node::trans_mean_log, Node::trans_var_log)));}
  if(!no_degrad)
    {setDegrad(exp(norm_sample(Node::degr_mean_log, Node::degr_var_log)));}
  for(eit=edge_in.begin();eit<edge_in.end();eit++){
    (*eit)->normalRandomParameters();
  }
}

int Node::getParameterNames(vector<string>& pnames)const{
  string s;
  uint n = input.size();
  s.assign(name);
  s.append("_s");
  pnames.push_back(s);
  s.assign(name);
  s.append("_d");
  pnames.push_back(s);
 for(uint i=0;i<n;i++){
   s.assign(input[i]->getName());
    switch(in_edges[i]){ 
    case TRANSCRIPTION:
      s.append("->");break;
    case COACTIVATION:
      s.append("=>");break;
    case DEGRADATION:
      s.append("-|");break;
   case SEQUESTRATION:
     s.append("|-|"); break;
    case HILL:
      s.append("->>");break;
    case WELL:
      s.append("=|");break;
    default:;
    }
    s.append(name);
    pnames.push_back(s);
  }
  return n+2;
}

void Node::printEdgeValues(std::ostream& out)const{
  string s;
  uint n = input.size();
  float sum=0.0;
  for(uint i=0;i<n;i++){
    sum += getParameterValues(i,0);
  }
  sum += src()+degrad();
  for(uint i=0;i<n;i++){
    float val= log(getParameterValues(i,0)/sum);
    s.assign(input[i]->getName());
    switch(in_edges[i]){ 
    case TRANSCRIPTION:
      s.append("->");break;
    case COACTIVATION:
      s.append("=>");break;    
    case DEGRADATION:
      s.append("-|");break;
   case SEQUESTRATION:
     s.append("|-|"); break;
  case HILL:
     s.append("->>"); break;
  case WELL:
     s.append("=|"); break;
    default:;
    }
    s.append(name);
    out<<s<<" "<<val<<" "<<log(getParameterValues(i,0))<<" "<<log(src())<<" "<<log(degrad())<<endl;
  }
}

void Node::print(ostream& out)const{
  out<<name<<": "<<state<<endl;
}

void Node::printParams(ostream& out)const{
  ParamList::const_iterator pit;
  out<<src()<<" "<<degrad()<<" ";
  for(pit=params.begin();pit<params.end();pit++){
    out<<*pit<<" ";
  }
}

int Node::removeInput(int i){
  int from = param_idx[i];
  int to = ((unsigned int) i) < param_idx.size()-1?param_idx[i+1]-1:params.size()-1;
  input.erase(input.begin()+i);
  in_edges.erase(in_edges.begin()+i);
  param_idx.erase(param_idx.begin()+i);
  params.erase(params.begin()+from,params.begin()+to);
  return input.size();
}

int Node::removeOutput(Node *o){
  for(uint i=0;i<children.size();i++){
    if(children[i]==o){
      children.erase(children.begin()+i);return 1;
    }
  }
  return 0;
}

//#define NODE_MAIN
#ifdef NODE_MAIN
int main(int argc, char *argv[]){
  Node n;
  srand(time(NULL));
  n.init();
  n.setState(1);
  while(!n.flat()){
    n.step(0.1);
    cout<<n.getState()<<endl;
  }
  return 1;
}

#endif
