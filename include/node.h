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

#ifndef __NODE_H__
#define __NODE_H__


#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "macros.h"
#include "edge.h"



class Edge;
class Node;


typedef std::vector<Node *> NodeList;
typedef std::vector<EdgeT> EdgeTypeList;
typedef std::vector<float> ParamList;
typedef std::vector<int> IndexList;


typedef std::vector<Edge> EdgeList; 
typedef std::vector<Edge*> EdgePtList;



class Trajectory{
 protected:
  float *buffer;
  int size;
  int nb;
  int top;
  int bot;//bottom

 private:
  int idx(int i)const{return (bot+i)%size;} 
 public:
  Trajectory(int si=5000){size=si;buffer = new float[size];clear();}
  ~Trajectory(){delete buffer;}
  void addPoint(float p)
  {if(nb==size){bot=(bot+1)%size;}else{nb++;}buffer[top]=p;top=(top+1)%size;}
  int getNbElements()const{return nb;};
  int getSize()const{return size;}
  float slope(int n)const{if(n>nb)return 100;return (buffer[idx(nb-1)]-buffer[idx(nb-n)])/n;}
  void clear(){nb=bot=top=0;}
};




class Node{

 protected:
  NodeList input;

  EdgeTypeList in_edges;
  ParamList params;
  IndexList param_idx;
  IndexList input_idx;
  NodeList children;

  EdgePtList edge_in;
  EdgePtList edge_out;

  ParamList *all_parameters;

  float state; // activity level



  int src_idx;
  int degrad_idx;

  char name[80];
  int idx;


  Trajectory traj;
  bool constrained;
  bool no_degrad;
  bool no_source;


 public:
  static float logmin;// for uniform param sampling
  static float logmax;// for uniform param sampling

  static float trans_mean_log;
  static float trans_var_log;

  static float src_mean_log;
  static float src_var_log;

  static float degr_mean_log;
  static float degr_var_log;

  static int nb_instances;

 public:
  Node(){state=1;no_source=no_degrad=constrained=false;idx=nb_instances++;};
  Node(const char *nname,int ind)
    {state=1;strcpy(name,nname);no_source=no_degrad=constrained=false;idx=ind;nb_instances++;}
  Node(const char *nname)
    {state=1;strcpy(name,nname);no_source=no_degrad=constrained=false;idx=nb_instances++;}
  //  Node(ParamList *params);
  ~Node(){}
  // void init();
  float getState()const{return state;}
  const char *getName()const{return name;}
  int getIndex()const{return idx;}
  void setState(float s){if(!constrained)state = s;}
  void setSource(float s){if(!no_source)(*all_parameters)[src_idx] = s;}
  //  float getSource()const{return no_source?0:(*all_parameters)[src_idx];}


  float f(float y=-1)const;
#ifdef DISCRETE
  void step(float dt){if(!constrained){state =f(state);}traj.addPoint(state);};
#else
  void step(float dt){if(!constrained){state+=dt*f(state);}traj.addPoint(state);};
#endif
  void print(std::ostream& out)const;
  void printParams(std::ostream& out)const;
  void addOutput(Node *o){children.push_back(o);}
  int addOutput(Edge *e){edge_out.push_back(e);return edge_out.size();}
  int addInput(Edge *e){edge_in.push_back(e);return edge_in.size();}
  int addInput(Node *i,Edge *e, float p);
  int addInput(Node *i,Edge *e, float p1,float p2);
  int removeInput(int i);
  int removeOutput(Node *o);
  int flat(int n=200)const{return fabs(traj.slope(n))<0.00001;}
  void randomParams();
  void randomParamsNormal();
  void noDegradation(){setDegrad(0.0);no_degrad=true;}
  bool isDegraded()const{return !no_degrad;}
  void noSource(){setSource(0.0);no_source=true;}
  void constrain(float s){state = s;constrained= true;}
  void unconstrain(){constrained=false;}
  bool isConstrained()const{return constrained;}
  void clearTraj(){traj.clear();}

  int getNbIncomingEdges()const{return edge_in.size();}
  int getNbParameters()const;
  float getParameterValues(int i,int j)const{return edge_in[i]->getValue(j);}
 
  int initParameters(ParamList& params);
  float src()const{return no_source?0:(*all_parameters)[src_idx];}
#ifndef OPTIMIZE_NODE_PARAMS 
  float degrad()const{return no_degrad?0:(*all_parameters)[degrad_idx];}
  void setDegrad(float s){if(!no_degrad)(*all_parameters)[degrad_idx] = s;}
#else
  float degrad()const{return no_degrad?0.0:1.0;}
#endif

  static float randParam();
  int getParameterNames(std::vector<std::string>& pnames)const;
  void printEdgeValues(std::ostream& out)const;
  int hasParamIdx(int i)const
  {if(src_idx==i)return 1;if(degrad_idx==i)return 2;return 0;}
  int setParameters(const float *p){
    uint s = edge_in.size(); int k=2;
    setSource(p[0]);setDegrad(p[1]);
    for(uint i=0;i<s;i++)k+=edge_in[i]->setValues(p+k);
    //   setParametersOS(p);
    return k;
  }
  int getParameters(float *p)const{
    uint s = edge_in.size(); p[0]=src();p[1]=degrad();int k=2;
    for(uint i=0;i<s;i++)
      k+= edge_in[i]->getValues(p+k);
    return k;
  }

  // for continuous time models
#ifndef DISCRETE
  void initJacobianRow(Vector& vec)const;
  void updateJacobianRow(float *row,std::vector<int>& nmap)const;
#endif

  // model specific functions

#ifdef GATED_GROWTH
  void gateState(float gating){state = gating/(1+exp(-gating*(state-0.5*gating)/gating));} 
#endif

#ifdef MENDOZA
  float getOmega(float& sum_p,float& scal_p,float& sum_n,float& scal_n,int& n)const;
#endif
#ifdef SIGMOID 
  float getPotential()const;
#endif

};


inline float Node::randParam(){
  return exp(logmin+RND(logmax-logmin));
}
#endif
