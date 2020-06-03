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


#ifndef __PARAMFITTER_H__
#define __PARAMFITTER_H__

#include "rnet.h"
#include <iostream>
#include "Vector.h"


#define MAX_NB_NETS 30


class ParamFitter;


class ParamFitter{
 protected:
  
  RegNet *nets[MAX_NB_NETS];
  float targets[MAX_NB_NETS][2];
  float span[MAX_NB_NETS][2];
  int nb_nets;
  

  Vector params;
  Vector log_params;


  static ParamFitter* instance;

 protected:
  int readAsciiLine(std::ifstream& inf,int n,int pos_score,std::vector<int>& zeros);
  int readBinaryLine(std::ifstream& inf,int n,std::vector<int>& zeros);
  int readMatHeader(std::ifstream& inf);
 public:
  ParamFitter(){instance=NULL;};
  ParamFitter(const char *netfname, int argc, const char *expfiles[]);
  ParamFitter(const char *netfname, const char *explistfile);
  ~ParamFitter();
  void init(const char *netfname,int nb_gen);
  void setTargets(const char *fname,int wspan);
  int eval(float& res);
  float evalSA(const float *params);
  int eval(float& res, const Vector& par);
  float eval(Vector& el_off,Vector& el_on,int leaveout=-1);
  int getSteadyStates(Vector& o_off,Vector& o_on)
  {return getOutputs(o_off,o_on,1);}
  int getOutputs(Vector& o_off,Vector& o_on,int allstates=0);

  // @param cost: pointer to member function (the cost function)
  float climb(int (ParamFitter::*cost)(float&),float step = 0.05,
	      float temp = 0,float cool=0.99);
  int readLogParams(const char *filename);
  void setParameters()
  {for(int i=0;i<nb_nets;i++){nets[i]->setParameters(params);}}
  void setParameters(const Vector& par)
  {for(int i=0;i<nb_nets;i++){nets[i]->setParameters(par);}}
  void setLogParameters(const Vector& par)// to check
  {setLogParameters(par.GetArray());}
  void setLogParameters(const float *par) // to check
  {for(uint i=0;i<log_params.Size();i++)log_params[i]= par[i];
    params = log_params.Exp();setParameters();}
  void setLogParameters(const double *par)
  {for(uint i=0;i<log_params.Size();i++)log_params[i]=(float)par[i];
    params = log_params.Exp();setParameters();}
  void setEdgeValue(uint ind, float v)
  {for(int i=0;i<nb_nets;i++){nets[i]->setEdgeValue(ind,v);}}
  void printLogPar(std::ostream& out)
  {for(uint i=0;i<log_params.Size();i++){out<<log_params[i]<<" ";}}
  void getParameterNames(std::vector<std::string>& pnames)const
  {nets[0]->getParameterNames(pnames);}
  void getParameterName(int i, char *name)const
  {nets[0]->getParameterName(i,name);}
  void randomParameters(int gauss=0)
  {if(gauss){nets[0]->randomParamsNormal();}else{nets[0]->randomParams();}
    nets[0]->getParameters(params);for(int i=1;i<nb_nets;i++)
  { nets[i]->setParameters(params);}log_params = params.Log();   }
  int getNbParameters(){return nets[0]->getNbParameters();}
  int getNbEdges(){return nets[0]->getNbEdges();}
  int checkParameters(const char *filename,int pos_score=1,int allstates=0,int map_leaveout=-1);
  void fitParameters(int nbtrials=10);//number of seeds for simulated annealing
  void checkEdgeRelevance(const char *fname);
  void generateData(int n,std::ostream& out=std::cout);//n:number of data points
  int distinguishable(float& res);
  void listEdges(const char *fname)const;
  int setNetworks(RegNet **networks, int n);
  void printTopology(std::ostream& out)const{nets[0]->printTopology(out);}
  void printExperimentList(std::ostream& out);
  void readExperimentList(const char *expfname);
  static ParamFitter* getInstance(){return instance;}

#ifdef WITH_ASA
  static float asaCostFunction(double *par,int nb,int *cost_flag);
#endif

};

#endif
