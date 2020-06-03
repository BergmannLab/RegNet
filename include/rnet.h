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

#ifndef __RNET_H__
#define __RNET_H__


#include <string.h>
#include <iostream>
#include <set>
#include "node.h"
#include "Matrix.h"

#define remove_last(a) (a).resize((a).size()-1)


#define MAX_STEPS 100000

#ifdef MENDOZA
#define FAC 0.08 //factor to have adequate units
#else
#ifdef SIGMOID
#define FAC 0.08 //factor to have adequate units
#else
#define FAC 1
#endif
#endif
/* #define insert(x) push_back((x)) */
/* #define erase(x) ((x)) */

typedef unsigned int uint;

typedef std::vector<Node *> NodeSet;
typedef std::vector<float> NodeStates;

typedef std::set<int> NodeIndexList;


typedef std::vector<Edge> NetEdgeList; 

class RegNet{
  
 protected:
  NodeSet nodes;
  EdgePtList edges;

  NodeSet shaded_nodes;
  NetEdgeList shaded_edges;

  int size;
  int nb_params;

  NodeIndexList input;
  NodeIndexList output;

  NodeStates offstates; 
  bool valid_offstates;

  NodeStates onstates;  
  bool valid_onstates;

  Matrix jacobian;
  Vector state;

  ParamList parameters;		/**< all the parameters of the network */
#ifdef GATED_GROWTH
  mutable float gating_beta;
  static float default_beta;
#endif


  std::vector<int> subspace; 

 public:
  RegNet();
  RegNet(const RegNet& net); 
  ~RegNet();
  const Node *getNode(int i)const{if(i<getNbNodes())return nodes[i];else return NULL;}
  const Edge *getEdge(int i)const{if(i<getNbEdges())return edges[i];else return NULL;}
  Node *findNode(const char *name)const;
  int findNodeIndex(const char *name)const;
  int readFromFile(const char *filename,float edge_default=-1.0);
  void addEdge2(Node *from, Node *to, EdgeT type,float p=-1);
  void addEdge3(Node *from, Node *from2,Node *to, EdgeT type,float p=-1);
  void addEdge(EdgePt& edge);
  // void addShadedEdge(Node *from, Node *to, EdgeT etype, float p);
  Node *addNode(const char *name);
  void addShadedNode(Node *n){shaded_nodes.push_back(n);}
  void print(std::ostream&  out=std::cout);
  int steadyState(int n)const;
  void step(float dt=0.1);
  void getUnconstrainedDerivatives(float *der);
  void printState(std::ostream& out=std::cout)const;
  void printExperiment(std::ostream& out=std::cout)const;
  int toSteadyState(int maxstep,int verbose=0);
  int findSteadyState(int maxstep,float tol);
  void randomTopology(int nb_nodes,float connectivity);
  void randomParams();
  void randomParamsNormal();
  void setAllStates(float s);
  void noInput(){input.clear();}
  //does not check if already in
  void asInput(int n){input.insert(n);} 
  //does not check if already in
  void asOutput(int n){output.insert(n);}
  int defineOnOffStates(const char *filename);
  int defineOnOffStates(std::ifstream& inf);
  void resetSteadyStates(){valid_onstates=valid_offstates=false;}
  void resetInputOutput(){removeAllInputs();output.clear();}
  const NodeIndexList& getInputIdx()const{return input;}
  void asInput(int n,float v_off,float v_on)
  {asInput(n),offstates[n]=v_off*FAC;onstates[n]=v_on*FAC;resetSteadyStates();}
  void removeAllInputs();
  int removeInput(int i)
  {nodes[i]->unconstrain();if(input.erase(i)){resetSteadyStates(); return 1;}return 0;}
  void initJacobian();
  void updateJacobian();
  Matrix& getJacobian(){return jacobian;}
  void resizeNewton();
  void updateStateVector()
  {for(uint i=0;i<nodes.size();i++){
      if(!nodes[i]->isConstrained())state[nmap(i)]=nodes[i]->getState();}}
  void updateStates()
  {for(uint i=0;i<nodes.size();i++){
      if(!nodes[i]->isConstrained())nodes[i]->setState(state[nmap(i)]);}}
  void updateStates(float *statevec)
  {for(uint i=0;i<nodes.size();i++){
      if(!nodes[i]->isConstrained())nodes[i]->setState(statevec[nmap(i)]);}}

  void getTimeDerivative(Vector& v) //v is assumde to be of the correct size
  {for(uint i=0;i<nodes.size();i++){
      if(!nodes[i]->isConstrained())v[nmap(i)]=nodes[i]->f();}}
  void turnOn();
  void turnOff();
  int fillOnStates(int trace=0);
  int fillOffStates(int trace=0);
  int setOffState();
  int setOnState();
  void print(NodeStates *s, std::ostream& out)const
  {for(uint i=0;i<s->size();i++)out<<(*s)[i]<<" ";out<<std::endl;}
  int fillStates(Vector& v, int cnt)const
  {for(int i=0;i<getNbNodes();i++)v[cnt++]=nodes[i]->getState();return cnt;}
  void printTopology(std::ostream& out);
  void printParams(std::ostream& out=std::cout)const;
  int nmap(int i){return subspace[i];}
  int inputMap(int i){int j=0;
    for(NodeIndexList::const_iterator it=input.begin();it!=input.end();it++){if(*it==i)return j;j++;}return -1;}
  float getState(int i)const{return nodes[i]->getState();}
#ifndef GATED_GROWTH
  float getOutput()const{return getState(*(output.begin()))/FAC;}//only the first output is returned
#else

  float gate(float state)const{return gating_beta/(1+exp(-gating_beta*(state-0.5*gating_beta)/gating_beta));}
  float getOutput()const{
#ifdef ADAPTIVE_GG
    gating_beta = log(parameters[nb_params-1])+10;// params have been exp, +10 because default values is more around 8 than -2. 
#endif
    return gate(getState(*(output.begin()))/FAC);}//only the first output is returned
#endif
  int getOutputIndex()const{return *(output.begin());}
  int getParameters(Vector& params);
  int setParameters(const Vector& params)
  {resetSteadyStates();return setParameters(params.GetArray());}
  int setParameters(const float *params);
  int setParameters(const double *params);
  int getNbEdges()const{int s=0;for(uint i=0;i<nodes.size();i++)//return edges.size()
			     {s+=nodes[i]->getNbIncomingEdges();}return s;}
  int getNbNodes()const{return nodes.size();};

  int setEdgeValue(uint i, float v)
  {resetSteadyStates();return edges[i]->setValue(v);}
  float getEdgeValue(uint i)const{return edges[i]->getValue();}
  int getNbParameters()const{return nb_params;}
  int getParameterName(int i, char *name)const;
  int getParameterNames(std::vector<std::string>& pnames)const;
  void printEdgeValues(std::ostream& out)const;
  int getEdgeParamIndex(const EdgePt *ep)const;
  int getSubspaceSize(){return subspace.size();}
 };



#endif
