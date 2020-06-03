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


#include "rnet.h"
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <stdio.h>

using namespace std;

#ifdef GATED_GROWTH
float RegNet::default_beta = 10;
#endif


RegNet::RegNet(){
  size=nb_params=0;
  subspace.clear();
}

RegNet::RegNet(const RegNet& net){
  for(int i=0;i<net.getNbNodes();i++){
    addNode(net.getNode(i)->getName());
  }
  cout<<net.getNbEdges()<<endl;
  for(int i=0;i<net.getNbEdges();i++){
    //    cout<<i<<" "<<net.getNbEdges()<<endl;
    const Edge *e = net.getEdge(i);
    EdgeT t = e->getType();
    switch(t){
    case TRANSCRIPTION:
    case DEGRADATION:
    case SEQUESTRATION:
    case HILL:
      addEdge2(findNode(e->origin()->getName()),
	       findNode(e->target()->getName()),t,e->getValue());
      break;
    case AGGREGATE:
    case DISAGGREGATE:
    case WELL:
    case COACTIVATION:
      addEdge3(findNode(e->origin()->getName()),findNode(e->auxil()->getName()),
	       findNode(e->target()->getName()),t,e->getValue());
      break;
    default:break;
    }
  }
  valid_offstates=false;
  valid_onstates=false;
  offstates.resize(getNbNodes());
  onstates.resize(getNbNodes());    
}


RegNet::~RegNet(){}

Node *RegNet::findNode(const char *name)const{
  NodeSet::const_iterator it;
  for(it=nodes.begin();it!=nodes.end();it++){
    if(!strcmp((*it)->getName(),name)){
      return *it;
    }
  }
  return NULL;
}

int RegNet::findNodeIndex(const char *name)const{
  for(uint i=0;i<nodes.size();i++){
     if(!strcmp(nodes[i]->getName(),name)){
       return i;
     }
  }
  return -1;
}

void RegNet::randomTopology(int nb_nodes, float connectivity){
  char nname[10];
  for (int i=0;i<nb_nodes;i++){
    sprintf(nname,"n%d",i);
    addNode(nname);
  }
  for (int i=0;i<nb_nodes-1;i++){
    int nedge = (int)ceil(RND(2*connectivity-1));
    vector<int> cand;
    for(int j=0;j<nb_nodes;j++){if(j!=i)cand.push_back(j);}
    for(int j=0;j<nedge;j++){
      int l = RND_INT(cand.size());
      int k = cand[l];
      cand[l] = cand[cand.size()-1];
      cand.pop_back();
      EdgeT et = RND_INT(2)?TRANSCRIPTION:DEGRADATION;
      addEdge2(nodes[i],nodes[k],et);
    }
  }
  offstates.resize(nodes.size());
  onstates.resize(nodes.size());
  state.Resize(nodes.size());
}


int RegNet::readFromFile(const char *fname,float edge_default){
  int cnt=0;
  Node *fn,*tn,*fn2;
  ifstream inf(fname);
  if(!inf.is_open()){cerr<<"cannot open file "<<fname<<endl;return 0;}
  char from[80];
  char to[80];
  char edge[80];
  char *f1,*f2;
  while(!inf.eof()){
    char c = inf.peek();
    if(c=='\n'){inf.get(c);c = inf.peek();}
    if(c == '#'){//comment
       inf.ignore(1000,'\n');
       continue;
    }
    inf>>from>>edge>>to;
    if(inf.fail()){continue;}
    f1 = strtok(from,",");
    f2 = strtok(NULL,",");
    fn = findNode(f1);
    tn = findNode(to);
    if(!fn){
      fn=addNode(f1);
    }
    if(!tn){
      tn=addNode(to);
    }
    if(f2){ //complex edge
      fn2 = findNode(f2);
      if(!fn2){
	fn2 =addNode(f2);
      }
      addEdge3(fn,fn2,tn,Edge::edgeType(edge),edge_default);
    }
    else{// simple edge
      addEdge2(fn,tn,Edge::edgeType(edge),edge_default);
    }
    cnt++;
  }
  offstates.resize(nodes.size());
  onstates.resize(nodes.size());
  state.Resize(nodes.size());
  inf.close();
#ifdef GATED_GROWTH
#ifdef ADAPTIVE_GG
   parameters.push_back(RegNet::default_beta);
   nb_params += 1;
#else
  gating_beta = RegNet::default_beta;
#endif

#endif

  return cnt;
}

Node *RegNet::addNode(const char *name){
  Node *no = new Node(name,size++);
  nb_params +=no->initParameters(parameters);; 
  nodes.push_back(no);
  return no;
}

void RegNet::addEdge2(Node *from, Node *to, EdgeT etype, float p){
  Edge *e= new Edge(to,from,etype);
  e->initParameters(parameters); 
  edges.push_back(e);
  switch(etype){
  case TRANSCRIPTION: //single input single output arrows
  case DEGRADATION:
  case SEQUESTRATION:
  case HILL:  
    to->addInput(e);
    from->addOutput(e);
    break;
  default:;
  }
  nb_params += e->getNbParameters();
  resetSteadyStates();
}


void RegNet::addEdge3(Node *from, Node *aux,Node *to, EdgeT etype, float p){
  Edge *e= new Edge(to,from,etype);
  e->initParameters(parameters);
  edges.push_back(e);
  e->setAux(aux);
  switch(etype){
  case AGGREGATE:
    to->addInput(e);
    from->addInput(e);
    aux->addInput(e);
    from->addOutput(e);
    aux->addOutput(e);
    to->noSource();
    to->noDegradation();
    break;
  case DISAGGREGATE:
    to->addInput(e);
    from->addInput(e);
    aux->addInput(e);
    from->addOutput(e);
    break;
  case WELL:
  case COACTIVATION:
    to->addInput(e);
    from->addOutput(e);
    aux->addOutput(e);
    break;
  default:;
  }
  nb_params += e->getNbParameters();
  resetSteadyStates();
}



int RegNet::getEdgeParamIndex(const EdgePt *ep)const{
  int i = 0;
  int idx = 0;
  while(nodes[i]!=ep->pto){
    idx+=2+nodes[i]->getNbIncomingEdges();
    if(i++==size){
      cerr<<"edge not found"<<endl;
      return -1;
    }
  }
  return idx+2+ep->nb;
}

void RegNet::print(ostream& out){
  NodeSet::const_iterator it;
  for(it=nodes.begin();it!=nodes.end();it++){
    (*it)->print(out);
    out<<endl;
  }
}

void RegNet::printTopology(ostream& out){
  for(uint i=0;i<edges.size();i++){
    edges[i]->print(out);
    out<<endl;
  }
}
 
int RegNet::steadyState(int n)const{
  NodeSet::const_iterator it;
  for(it=nodes.begin();it!=nodes.end();it++){
    if(!(*it)->flat(n)){
      return 0;
    }
  }
  return 1;
}

void RegNet::step(float dt){
  NodeSet::iterator it;
  for(it=nodes.begin();it!=nodes.end();it++){
    (*it)->step(dt);
  }
}

void RegNet::getUnconstrainedDerivatives(float *yd){
  for(uint i=0;i<nodes.size();i++){
    if(!nodes[i]->isConstrained()){
      yd[nmap(i)]=nodes[i]->f(-1);
    }
  }
}

void RegNet::printState(ostream& out)const{
  NodeSet::const_iterator it;
  for(it=nodes.begin();it!=nodes.end();it++){
    out<<(*it)->getState()<<" ";
  }
  out<<endl;
}

void RegNet::printParams(ostream& out)const{
  for(uint i=0;i<parameters.size();i++){
    out<<parameters[i]<<" ";
  }
  out<<endl;
}

void RegNet::printEdgeValues(ostream& out)const{
  for(uint i=0;i<nodes.size();i++){
    nodes[i]->printEdgeValues(out);
  }
}

int RegNet::toSteadyState(int maxstep, int verbose){
  int cnt = 0;
  step();
  while(!steadyState(200) && cnt++ < maxstep){
    step();
    if(verbose)if(!(cnt%verbose)){printState(cout);}
  }
  if(cnt>= maxstep){
    //   cerr<<"cannot reach steady state"<<endl;
    return 0;
  }
  else{
    return 1;
  }
}

#ifndef CVODE_INTEGRATION

#ifdef STATIC_NETWORK

int RegNet::findSteadyState(int maxsteps, float tol){
  resizeNewton();
  Vector x(input.size()+1);//+ cst term
  int n = nodes.size()-input.size();
  Matrix B(n,input.size()+1);
  Matrix eye(n,n);
  eye.Identity();
  jacobian.Zero();
  B.Zero();
  for(uint i=0;i<edges.size();i++){
    Node *tar = edges[i]->target();
    Node *ori = edges[i]->origin();
    int sign = (edges[i]->getType()==TRANSCRIPTION)?1:-1;
    if((! tar->isConstrained()) && (! ori->isConstrained())){
      jacobian(nmap(tar->getIndex()),nmap(ori->getIndex()))=sign*edges[i]->getValue(); 
    }
    if(! tar->isConstrained() && ori->isConstrained()){
      B(nmap(tar->getIndex()),inputMap(ori->getIndex())) = sign*edges[i]->getValue(); 
    }
  }
  for(uint i=0;i<nodes.size();i++){
    if(!nodes[i]->isConstrained()){
	float comb = nodes[i]->src()-nodes[i]->degrad();
	if(nodes[i]->isDegraded()){
	  B(nmap(nodes[i]->getIndex()),input.size())=comb;
	}else{
	   B(nmap(nodes[i]->getIndex()),input.size())=(comb<0)?0:comb;
	}
    }
  }
  NodeIndexList::const_iterator it;
  int i=0;
  for(it=input.begin();it!=input.end();it++){
    x(i++)=nodes[*it]->getState();
  }
  x(i)=1;
  state = ((eye -jacobian).Inverse())*B*x;
  if (!Matrix::IsInverseOk()){
    cerr<<"bad inversion"<<endl;
    return 0;// trial
  }
  else{
    updateStates();  
    return 1;
  }  
}

#else //STATIC_NETWORK
#ifdef DISCRETE
int RegNet::findSteadyState(int maxsteps, float tol){
  return toSteadyState(maxsteps,0);
}

#else //home-made integration using Newton, better to use CVODE
int RegNet::findSteadyState(int maxsteps, float tol){
  Vector f(nodes.size());
  Vector dx(nodes.size());
  int cnt=0;
  int notrunc = 1;
  float change = FLT_MAX;
  float step_size = 0.05;
  float tol2 = tol*tol;
  resizeNewton();
  while(cnt<maxsteps && change > tol2){
    notrunc = 1;
    updateStateVector();
      getTimeDerivative(f);
    updateJacobian();
    dx = jacobian.Inverse()*f;//dx/dt = df/dt*dx/df
    state = state - (dx*step_size);
    for(uint i=0;i<state.Size();i++){
      if(state[i]<0){
	state[i] =0;
	notrunc=0;
      }
    }
    change = dx.Norm2();
    if (!Matrix::IsInverseOk()){
      cerr<<"bad inversion"<<endl;
      return 0;// trial
    }
    else{
      updateStates();  
    }  
    cnt++;
  }
  return (cnt<maxsteps)*notrunc;
}
#endif
#endif
#endif
void RegNet::randomParams(){
   NodeSet::iterator it;
   for(it=nodes.begin();it!=nodes.end();it++){
     (*it)->randomParams();
   }
   resetSteadyStates();
}

void RegNet::randomParamsNormal(){
   NodeSet::iterator it;
   for(it=nodes.begin();it!=nodes.end();it++){
     (*it)->randomParamsNormal();
   }
   resetSteadyStates();
}

void RegNet::setAllStates(float s){
NodeSet::iterator it;
   for(it=nodes.begin();it!=nodes.end();it++){
     (*it)->setState(s);
   }
}

void RegNet::removeAllInputs(){
  NodeIndexList::iterator it;
  for(it=input.begin();it!=input.end();it++){
    nodes[*it]->unconstrain();
  }
  resetSteadyStates();
  input.clear();
}


int RegNet::getParameters(Vector& params){
  int nb =getNbParameters();
   memcpy(params.GetArray(),&(parameters[0]),nb*sizeof(float));
  return nb;
}

int RegNet::getParameterNames(vector<string>& pnames)const{
  int n=0;
  for(uint i=0;i<nodes.size();i++){
    n+=nodes[i]->getParameterNames(pnames);
  }
  return n;
}


int RegNet::getParameterName(int ind, char *name)const{
  for(uint i=0;i<nodes.size();i++){
    int a = nodes[i]->hasParamIdx(ind);
    if(a>0){
      strcpy(name, nodes[i]->getName());
      switch(a){
      case 1: strcat(name,"_s");return 1;
      case 2: strcat(name,"_d");return 1;
      default:strcat(name,"_?");return 2;
      }
      return 1;
    }
  }
  for(uint i=0;i<edges.size();i++){
    int a = edges[i]->hasParamIdx(ind);
    if(a>0){
      char atype[5];
      Edge::arrowType(edges[i]->getType(),atype);
      strcpy(name,edges[i]->origin()->getName());
      strcat(name,atype);
      strcat(name,edges[i]->target()->getName());
      if(a>1){
	sprintf(atype,"_%d",a);
	strcat(name,atype);
      }
    }
    return 1;
  }
  return 0;
}

int RegNet::setParameters(const float *par){
  int nb = parameters.size();
  memcpy(&(parameters[0]),par,nb*sizeof(float));
  return nb;
}

int RegNet::setParameters(const double *par){
  int nb = parameters.size();
  for(int i=0;i<nb;i++){
    parameters[i] = (float)par[i];
  }
  return nb;
}

// unkonwn nodes are ignored
int RegNet::defineOnOffStates(const char *filename){
  ifstream inf(filename);
  //  cout<<filename<<endl;
  if(!inf.is_open()){cerr<<"cannot open file "<<filename<<endl;return 0;}
  defineOnOffStates(inf);
  inf.close();
  return 1;
}

int RegNet::defineOnOffStates(ifstream& inf){
  char name[80];
  float onval, offval;
  int fn;
  int inp=0,outp=0;
  resetInputOutput();
  while(!inf.eof()){
    char c = inf.peek();
    if(c=='#'){inf.ignore(1000,'\n');continue;}//comment
    if(c=='='){inf.ignore(1000,'\n');break;}//end of network
    if(c=='+'){inf.get();}
    inf>>name;
    fn = findNodeIndex(name);
    if(!inf.fail() && fn>=0){
      if(c=='+'){//no degradation
	nodes[fn]->noDegradation();
      }
      else{
	inf>>offval>>onval;
	if(inf.fail()){
	  asOutput(fn);
	  outp=1;
	  if(!inf.eof()){inf.clear();
	  }else{continue;}
	}
	else{
	  inp=1;
	  asInput(fn,offval,onval);
	}
      }
    }
    else{
      continue;
    }
  }
  resizeNewton();
  return inp*outp;
}


void RegNet::printExperiment(ostream& out)const{
  NodeIndexList::const_iterator it;
  for(it=input.begin();it!=input.end();it++){
    int i= *it;
    out<<nodes[i]->getName()<<" "<<offstates[i]/FAC<<" "<<onstates[i]/FAC<<endl;
  }
  for(int i=0;i<getNbNodes();i++){
    if(!nodes[i]->isDegraded()){
      out<<"+"<<nodes[i]->getName()<<endl;
    }
  }
  for(it=output.begin();it!=output.end();it++){
    int i= *it;
    out<<nodes[i]->getName()<<endl;
  }
  out<<"=============================="<<endl;
}

void RegNet::turnOn(){
  NodeIndexList::iterator sit;
  for(sit=input.begin();sit!=input.end();sit++){
    if((*sit)<(int)nodes.size()){
      nodes[*sit]->constrain(onstates[*sit]);
    }
  }
}

void RegNet::turnOff(){
  NodeIndexList::iterator sit;
  for(sit=input.begin();sit!=input.end();sit++){
    if((*sit)<(int)nodes.size()){
      nodes[*sit]->constrain(offstates[*sit]);
    }
  }
}
void RegNet::resizeNewton(){
  subspace.resize(nodes.size());
  int n=nodes.size()-input.size();
  jacobian.Resize(n,n,false);
  jacobian.Zero();
  //  initJacobian();
  state.Resize(n);
  int j=0;
  for(uint i=0;i<nodes.size();i++){
    if(!nodes[i]->isConstrained()){
      subspace[i]=j++; 
    }
    else{
      subspace[i]=-1;
    }
  }
}
#ifndef DISCRETE
// not used anymore
void RegNet::initJacobian(){
  int n = nodes.size()-input.size();
  jacobian.Resize(n,n);
  Vector vec(nodes.size());
  for(uint i=0;i<nodes.size();i++){
    if(!nodes[i]->isConstrained()){
      nodes[i]->initJacobianRow(vec);
      for(uint j=0;j<nodes.size();j++){
	if(!nodes[j]->isConstrained()){
	  vec[nmap(j)] = vec[j];
	}
      }
      jacobian.SetRow(vec,nmap(i));
    }
  }
}

void RegNet::updateJacobian(){
  for(uint i=0;i<nodes.size();i++){
    if(!nodes[i]->isConstrained()){
      float *r =jacobian.GetPointerToRow(nmap(i)); 
      nodes[i]->updateJacobianRow(r,subspace);    
    }
  }
}
#endif

int RegNet::fillOnStates(int trace){
  turnOn();
  if(toSteadyState(MAX_STEPS,trace)){
    //  onstates.resize(nodes.size());
    for (uint i=0;i<nodes.size();i++){
      onstates[i]=nodes[i]->getState();
    }
    valid_onstates = true;
    return 1;
  }
  valid_onstates = false;
  return 0;
}


int RegNet::fillOffStates(int trace){
  turnOff();
  if(toSteadyState(MAX_STEPS,trace)){
    //   offstates.resize(nodes.size());
    for (uint i=0;i<nodes.size();i++){
      offstates[i]=nodes[i]->getState();
    }
    valid_offstates = true;
    return 1;
  }
    valid_offstates = false;
  return 0;
}

int RegNet::setOffState(){
  if(!valid_offstates)fillOffStates();
  if(valid_offstates){
    for(uint i=0;i<nodes.size();i++){
      nodes[i]->setState(offstates[i]);
    }
    return 1;
  }
  return 0;
}

int RegNet::setOnState(){
  if(!valid_onstates)fillOnStates();
  if(valid_onstates){
    for(uint i=0;i<nodes.size();i++){
      nodes[i]->setState(onstates[i]);
    }
    return 1;
  }
  return 0;
}



