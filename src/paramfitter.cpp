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

#include <fstream>
#include <stdlib.h>
#include "paramfitter.h"
#include <float.h>

#define NEWTON_METHOD

using namespace std;

ParamFitter * ParamFitter::instance  =NULL;

ParamFitter::ParamFitter(const char *netfname, const char *explistfname){
  ifstream inf(explistfname);
  nb_nets= 0;
  if(!inf.is_open())
    {cerr<<"cannot open file "<<explistfname<<endl;}
  else{
    char fname[200];
    while(!inf.eof()){
      char c = inf.peek();
      if(c=='#'){  inf.ignore(200,'\n');continue;}//ignore comments
      if(c=='\n' || inf.eof()){continue;}
      inf>>fname;
      inf.ignore(200,'\n');//ignore rest of line
      if(nb_nets>=MAX_NB_NETS){
       cerr<<"To many networks "<<nb_nets<<
	 "MAX_NB_NETS should be increased"<<endl;
     }   
     Node::nb_instances=0;//needed to have a coherent node indexing with idx ??
     nets[nb_nets] = new RegNet();   
     nets[nb_nets]->readFromFile(netfname);
     nets[nb_nets]->defineOnOffStates(fname);
     nets[nb_nets]->setAllStates(10);
     nb_nets++;
    }
    params.Resize(nets[0]->getNbParameters());
  }
  instance =this;
}

ParamFitter::~ParamFitter(){
}

void ParamFitter::init(const char *netfname, int nb_gen){
 for(int i=0;i<nb_gen;i++){
   nets[i] = new RegNet();   
   nets[i]->readFromFile(netfname);
 }
 params.Resize(nets[0]->getNbParameters());
 log_params.Resize(nets[0]->getNbParameters());
 nb_nets = nb_gen;
 instance =this;
}


void ParamFitter::setTargets(const char *fname, int wspan){
  ifstream inf(fname); 
  char line[200];
  if(!inf.is_open()){cerr<<"cannot open file "<<fname<<endl;return;}
  int i=0;
  while(!inf.eof()){
    char c = inf.peek();
    if(c=='#'){inf.getline(line,200);continue;}//comment
    inf>>targets[i][0]>>targets[i][1];
    if(wspan){inf>>span[i][0]>>span[i][1];}
    inf.ignore(200,'\n');
    i++;
    if(i>=nb_nets){
      break;
    }
  }
  inf.close();  
  if(i<nb_nets){cerr<<"expected more target values"<<endl;}
 
}

void ParamFitter::readExperimentList(const char *expfname){
  ifstream inf(expfname);
  if(!inf.is_open()){cerr<<"cannot open file "<<expfname<<endl;return;}
  for(int i=0;i<nb_nets;i++){
    if(!nets[i]->defineOnOffStates(inf)){
      cerr<<"cannot read experiment "<<i<<endl;
    }
  }
  inf.close();
}

int ParamFitter::eval(float& res,const Vector& par){
   for (int i=0;i<nb_nets;i++){
     nets[i]->setParameters(par);
   }
   return eval(res);
}

float ParamFitter::evalSA(const float *lpar){
  int n = getNbParameters();
  Vector v(n);
  for (int j=0;j<n;j++){
    v[j]=exp(lpar[j]);
  }
  for (int i=0;i<nb_nets;i++){
    nets[i]->setParameters(v);
  }
  float res=FLT_MAX;
  eval(res);
  return res;
}

/**
 * @ returns 0 if no steady state found, -1 if non matching steady-state, 1 if matching steady-state
 */
int ParamFitter::eval(float& res){
  float sum=0;
  int ok=1;
  int found = 0;
  //  cout<<"eval"<<endl;
  for (int i=0;i<nb_nets;i++){
    nets[i]->setAllStates(0);
    nets[i]->turnOff();
#ifdef NEWTON_METHOD
    if(nets[i]->findSteadyState(1000,0.001)){
      found = 1;
    }
    else{
      nets[i]->setAllStates(1);
      nets[i]->turnOff();
      if(nets[i]->findSteadyState(1000,0.001)){
	found = 1;
      }
    }
#else
    found = nets[i]->toSteadyState(100000,50);
#endif
    if(found){
      //     cout<<"found"<<" "<<nets[i]->getOutput()<<endl;
      float diff = fabs(nets[i]->getOutput()-targets[i][0])/span[i][0];
      sum+=diff*diff;
      if(diff>1){ok=-1;}
    }
    else{
      return 0;
    }
    found = 0;
    nets[i]->setAllStates(0);
    nets[i]->turnOn();
#ifdef NEWTON_METHOD
    if(nets[i]->findSteadyState(10000,0.001)){
      found = 1;
    }
    else{
      nets[i]->setAllStates(1);
      nets[i]->turnOn();
      if(nets[i]->findSteadyState(1000,0.001)){
	found  =1;
      }
    }
#else
    found = nets[i]->toSteadyState(100000,50);
#endif
    if(found){
      //     cout<<"found 2 "<<" "<<nets[i]->getOutput()<<endl;
      float diff = fabs(nets[i]->getOutput()-targets[i][1])/span[i][1];
      sum+=diff*diff;
      if(diff>1){ok=-1;}
    }
    else{
      return 0;
    }
  }
  //  cout<<res<<endl;
  res = sqrtf(sum);
  return ok;
}

float ParamFitter::eval(Vector& el_off, Vector& el_on,int leaveout){
  float sum = 0;
  float diff=0;
  for (int i=0;i<nb_nets;i++){
    if(i!= leaveout){
      diff= (el_off[i]-targets[i][0])/span[i][0];
      sum+= diff*diff;
      diff= (el_on[i]-targets[i][1])/span[i][1];
      sum+= diff*diff;
    }
  }
  return sqrtf(sum);
}

int ParamFitter::distinguishable(float& res){
  Vector el_off,el_on;
  getOutputs(el_off,el_on);
  // el_on /= FAC;
  // el_off /= FAC;
  res = fabs(10+el_on.Min()-el_on.Max()) +fabs(5+el_off.Min()-el_off.Max());

  int ok= res < 1?1:-1;
  float m =min(el_on.Min(),el_off.Min());
  if(m<0){
    res -= 100*m;  
  }
  return ok;
}


int ParamFitter::getOutputs(Vector& o_off,Vector& o_on,int allstates){
  int m = allstates?nets[0]->getNbNodes():1;
  o_on.Resize(nb_nets*m);
  o_off.Resize(nb_nets*m);
  int ok=1;
  int cnton=0;
  int cntoff=0;
  for (int i=0;i<nb_nets;i++){
    nets[i]->setAllStates(1);
    nets[i]->turnOff();
    if(nets[i]->findSteadyState(1000,0.001)){
      
      if(allstates){
	cntoff=nets[i]->fillStates(o_off,cntoff);	  
      }
      else{
	o_off[i] = nets[i]->getOutput();
      }
    }
    else{
	nets[i]->setAllStates(0);
	nets[i]->turnOff();
	if(nets[i]->findSteadyState(1000,0.001)){//second try
	  if(allstates){
	    cntoff=nets[i]->fillStates(o_off,cntoff);	  
	  }
	  else{
	    o_off[i] = nets[i]->getOutput();
	  }
      }
      else{
	o_off[i] = -1;
	ok=-1;
      }
    }
    nets[i]->setAllStates(1);
    nets[i]->turnOn();
    if(nets[i]->findSteadyState(1000,0.001)){
      if(allstates){
	cnton=nets[i]->fillStates(o_on,cnton);      }
      else{    
	o_on[i] = nets[i]->getOutput();
      }
    }
    else{
      nets[i]->setAllStates(0);
      nets[i]->turnOn();
      if(nets[i]->findSteadyState(1000,0.001)){//second try
	if(allstates){
	  cnton=nets[i]->fillStates(o_on,cnton);
	}
	else{
	  o_on[i] = nets[i]->getOutput();
	}
      }
      else{
	o_on[i] = -1;
	ok=-1;
      }
    }
  }
  return ok;
}


int ParamFitter::readLogParams(const char *fname){
 ifstream inf(fname);
 if(!inf.is_open()){cerr<<"cannot open file "<<fname<<endl;return 0;}
 int n = nets[0]->getNbParameters();
 log_params.Resize(n);
 for(int i=0;i<n;i++){
   inf>>log_params[i];
   if (inf.fail()){
     cerr<<"bad format in "<<fname
	 <<". Expecting "<<n<<" numbers"<<endl;
     return 0;
   }
 }
 inf.close();
 params = log_params.Exp();
 return n;
}


// early simulated annealing method, not used

float ParamFitter::climb(int (ParamFitter::*cost)(float&),float step,float temp,float cool){

  Vector p,diff,best_p,p1,vbest_p;
  float v,best_v,vbest_v;
  int val, best_val=0,vbest_val=0,cnt=1000;
  int niter  =50000;

  diff.Resize(params.Size());
  p.Resize(params.Size());
  p1.Resize(params.Size());
  vbest_p=best_p= log_params;
  v=vbest_v=best_v = FLT_MAX; 
  int i;
  for(int j=0;j<nb_nets;j++){
    nets[j]->setAllStates(1);
  }
  for(i=0;i<niter;i++){
    float st= step*RND(1);
    diff.Random();
    diff = diff*2-p1.One();
    diff /= (diff.Norm()/st);
     p = best_p+diff;
    p1 = p.Exp();
    setParameters(p1);
    val = (this->*cost)(v);

    //simulated annealing
    int mod = val?((v<best_v)?1:(RND(1)<exp((best_v-v)/temp))):0;  
    if(mod){
      if(v>best_v && vbest_v>best_v){ //keep the very best
	vbest_v = best_v;
	vbest_p = best_p;
	vbest_val = best_val;
      }
      best_v = v;
      best_p = p;
      best_val = val;
      temp *= cool;
      cnt = 1000;
    }
    if(val)if(!(cnt--)){break;}
  }
  if(!best_val && !vbest_val) return -FLT_MAX;//no steady-state found
  if(vbest_v >= best_v){
    log_params = best_p;
    return best_v*best_val;
  }
  else{
    log_params = vbest_p;
    return vbest_v*vbest_val;
  }
}

int ParamFitter::setNetworks(RegNet **networks, int n){
  if(n>=MAX_NB_NETS){
    cerr<<"can only accept "<<MAX_NB_NETS<<" networks"<<endl;
    return 0;
  }
  for(int i=0;i<n;i++){
    nets[i] = networks[i];
  }
  nb_nets=n;
  int p = nets[0]->getNbParameters();
  params.Resize(p);
  log_params.Resize(p);
  return 1;
}

void ParamFitter::generateData(int n, ostream& out){
  int np =getNbParameters();
  Vector el_off(np),el_on(np),sum_on(np),sum_off(np),sum2_on(np),sum2_off(np);
  Vector log_params_bak;
  randomParameters();
#ifdef GENERIC_CLIMB
  climb(&ParamFitter::distinguishable,0.4,50);
#endif 
 log_params_bak=log_params;
  Vector jiggle(np);
  for (int cnt=0;cnt<n;){
    jiggle.RandomNormal(0,0.03);
    log_params = log_params_bak+jiggle;
    params = log_params.Exp();
    setParameters();
    if(getOutputs(el_off,el_on)>0){
      //    cout<<el_off<<" "<<el_on<<endl;
      sum_on +=el_on;
      sum_off +=el_off;
      sum2_on += el_on^el_on;
      sum2_off += el_off^el_off;
      cnt++;
    }
  }
  //setting back optimal parameters
  log_params=log_params_bak;
  params = log_params.Exp();
  setParameters();

  sum_on /= n;
  sum_off /= n;
  sum2_on /= n;
  sum2_off /= n;

  sum2_on =  sum2_on - (sum_on ^sum_on);//now the variance
  sum2_off = sum2_off -(sum_off^sum_off);//now the variance
  for(int i=0;i<nb_nets;i++){
    out<<sum_off[i]<<" "<<sum_on[i]<<" "
	<<sqrtf(sum2_off[i])<<" "<<sqrtf(sum2_on[i])<<endl;
  }
}

void ParamFitter::printExperimentList(ostream& out){
  for(int i=0;i<nb_nets;i++){
    nets[i]->printExperiment(out);
  }
}

void ParamFitter::listEdges(const char *fname)const{
  ofstream out(fname);
  if(!out.is_open()){cerr<<"cannot open file "<<fname<<endl;return;}
  // edge relative values
  nets[0]->printEdgeValues(out);
  // param values
  out<<"log_params "<< log_params<<endl; 
  // param names
  vector<string> names;
  getParameterNames(names);
  for(uint i=0;i<names.size();i++){
    out<<names[i]<<" ";
  }
  cout<<endl;
  out.close();
}


void ParamFitter::checkEdgeRelevance(const char *fname){
  ifstream inf(fname);
  if(!inf.is_open()){cerr<<"cannot open file "<<fname<<endl;return;}
  int n = getNbParameters();
  int nbedge = getNbEdges();
  log_params.Resize(n);
  float ind;
  while(!inf.eof()){
    inf>>ind;
    cout<<ind;
    for(int i=0;i<n;i++){
      inf>>log_params[i];
    }
    if(!inf.fail()){
      Vector el_off,el_on;
      params = log_params.Exp();
      setParameters();
      getOutputs(el_off,el_on);
      // for all edges, set parameter to zero, get output and compare
        for(int i=0;i<nbedge;i++){
	float vbak = nets[0]->getEdgeValue(i);
	setEdgeValue(i,0);
	Vector s_off,s_on;
	getOutputs(s_off,s_on);
	setEdgeValue(i,vbak);
	float effect = sqrtf((s_off-el_off).Norm2()+(s_on-el_on).Norm2());
	cout<<" "<<effect;
      }
      cout<<endl;
    }
    else{
      cerr<<"failed to read parameters"<<endl;
    }
  }
  inf.close();
}

int ParamFitter::readAsciiLine(ifstream& inf,int n,int pos_score,vector<int>& zeros){
  char num[20];  
  float score;
  inf>>score;
    if(score<0 && pos_score){
      inf.ignore(1000,'\n');
      return 0;
    }
    for(int i=0;i<n;i++){
      inf>>num;
      if(!strcmp(num,"-inf")){
	log_params[i] = 0;
	zeros.push_back(i);
      }
      else{
	log_params[i]=atof(num);
      }
    }
    return 1;
}
int ParamFitter::readBinaryLine(ifstream& inf,int n,vector<int>& zeros){
  // if the format is right....
  //  inf.read((char* )log_params.GetArray(),n*sizeof(float));
  float val;
  for(int i=0;i<n;i++){
    inf.read((char* )&val,sizeof(float));
    log_params[i]=val;
  }
  return 1;
}

int ParamFitter::readMatHeader(ifstream& inf){
  char header[130];
  //  char c1,c2;
  int32_t dtype, dsize;
  inf.read(header,128);
  inf.read((char *)&dtype,4);
  cout<<"type is "<<dtype<<endl;
  if(dtype!=14){
    cout<<"wrong mat data type, expecting 14: "<<dtype<<endl;
  }
  inf.read((char *)&dsize,4);
  cout <<"data size is "<<dsize<<endl;
  return (int) dsize;
}


int ParamFitter::checkParameters(const char *fname, int pos_score,int allstates,int map_leaveout){
  vector<int> zeros;
  ifstream inf;
  int binary = 0;
  int ok;
  Vector opt_el_off(0),opt_el_on(0);
  if(!strcmp(fname+strlen(fname)-4,".par")){//binary file
    inf.open(fname, ios_base::in | ios_base::binary);
    binary =1;
    //    cout<<"binary file, expected float size: "<<sizeof(float)<<endl;
  }
  else{ //ascii file
    //  cout<<"ascii file"<<endl;
    inf.open(fname, ios_base::in);
  }
  if(!inf.is_open()){cerr<<"cannot open file "<<fname<<endl;return 0;}
  int n = getNbParameters();
  log_params.Resize(n);
  float min_score = FLT_MAX;
  while(!inf.eof()){
    if(binary){
      ok = readBinaryLine(inf,n,zeros);
    }
    else{
      ok = readAsciiLine(inf,n,pos_score,zeros);
    }
    if(!inf.fail()){
      if(ok){
	Vector el_off,el_on;
	params = log_params.Exp();
	for(uint i=0;i<zeros.size();i++){//putting all zeros
	  params[zeros[i]] = 0;
	}
	setParameters();
	getOutputs(el_off,el_on,allstates);

	if(map_leaveout>=0){ // only check output
	  float score = eval(el_off,el_on,map_leaveout);
	  if(score < min_score){
	    min_score=score;
	    opt_el_off = el_off;
	    opt_el_on = el_on;
	  }
	}
	else{ // writes down output
	  cout<<el_off<<" "<<el_on<<endl;
	}
      }
    }
    else{
      if(!inf.eof()){
      cerr<<"failed to read parameters"<<endl;
      }
    }
  }
  if(map_leaveout>=0){
	cout<<opt_el_off<<" "<<opt_el_on<<endl;
  }
  inf.close();
  return 1;
}

#ifndef WITH_ASA
void ParamFitter::fitParameters(int nbtrials){
  float v,vb=FLT_MAX;
  
  //  cout<<rand()<<endl; //for proper seeding check
  for(int i=0;i<nbtrials;i++){
    randomParameters();
    v=climb(&ParamFitter::eval,0.4,50);//was 0.4 10
    if(fabs(v)<vb || v> 0){
      cout<<v<<" ";
      printLogPar(cout);
      cout<<endl;
      if(v<0)vb=fabs(v);
    }
  }
}
#endif






#define SAMPLE_AND_CLIMB

#ifdef PARAMFITTER_MAIN

void usage(){

#ifdef CHECK_RESULTS
  // <
  cout<<"./checkres <network> <experimentslist> <parameters> [<pos_score> <allstates> <map>]"<<endl;
  cout<<"<pos_score> if one, discards negative score (not used)"<<endl;
  cout<<"<allstates> if one, outputs all network states (should be 0)"<<endl;
  cout<<"<map> only prints out the result of the optimal parameter vector "<<endl;

#else
#ifdef DATA_GENERATION
  cout<<"./gendata <network> <experimentslist> <outedgefile> "<<endl;
#else
  cout<<"./shademod <network> <experimentslist> <target>"<<endl;
#endif
#endif
  return;
}

int main(int argc,const char *argv[]){
  if(argc<=2){usage();return 1;}
  srand(time(NULL));
  ParamFitter pf(argv[1],argv[2]);
#ifdef DATA_GENERATION
  pf.generateData(100);
  pf.listEdges(argv[3]);
#else  

#ifdef CHECK_RESULTS
  int che = 1;
  if (argc>=4){
    che = atoi(argv[4]);
  }
  int allstates = 0;
  if (argc>=5){
    allstates = atoi(argv[5]);
  }
  int map = -1;
  if(argc >= 6){
    map = atoi(argv[6]);
  }
  if(argc >= 7){
    pf.setTargets(argv[7],1);
  }
  pf.checkParameters(argv[3],che,allstates,map);
 #else
#ifdef SAMPLE_AND_CLIMB
  //test
  // pf.init(argv[1],6);
  // pf.readExperimentList(argv[2]);
  //end test
 pf.setTargets(argv[3],1);
  float v,vb=FLT_MAX;
  for(int i=0;i<5;i++){
    pf.randomParameters();
    v=pf.climb(&ParamFitter::eval,0.4,50);//was 0.4 10
    cout<<"done"<<endl;
    //  if(v>0){
    cout<<v<<endl;
    if(fabs(v)<vb || v> 0){
      cout<<v<<" ";
      pf.printLogPar(cout);
      cout<<endl;
      if(v<0)vb=fabs(v);
    }
    if((i+1)%1000==0){cerr<<argv[1]<<": "<<i+1<<endl;}
  }

#else
  vector<string> names;
  Vector el_off,el_on;
  pf.readLogParams(argv[argc-1]);
  pf.setParameters();
  pf.getOutputs(el_off,el_on);
  pf.getParameterNames(names);
  for(uint i=0;i<names.size();i++){
    cout<<names[i]<<" ";
  }
  cout<<endl;
#endif
#endif
#endif
  return 0;
}

#endif
