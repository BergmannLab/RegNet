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

#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <iostream>
#include <stddef.h>
#include <stdio.h>
#include <sys/un.h>
#include <signal.h>
#include <string.h>
#include <errno.h>

#include "paramfitter.h"

using namespace std;

#define BUF_SIZE 100
#define NO_CONV_VAL 1000
#define NAME_LENGTH 30


int sock[2];
char socketname[100];

ParamFitter& loadNetworks(char *networkfile, char *explistfile, char *targetfile, int n=0){
  ParamFitter *pf;
  if(!n){
    pf= new ParamFitter(networkfile,explistfile);
  }
  else{
    pf = new ParamFitter();
    pf->init(networkfile,n);
    pf->readExperimentList(explistfile);
  }
  pf->setTargets(targetfile,1);
  return *pf;
}



void initServer(ParamFitter& pf){
  struct sockaddr_un name;
  sock[0] = socket(AF_LOCAL,SOCK_STREAM,0);
  if(sock[0]<0){cout<<"cannot create socket"<<endl;}
  name.sun_family = AF_FILE;
  strcpy(name.sun_path,socketname);
  size_t size = (offsetof (struct sockaddr_un, sun_path)
		 + strlen (name.sun_path) + 1);
  if (bind (sock[0], (struct sockaddr *) &name, size) < 0){
    cout<<"cannot bind"<<" "<<strerror(errno)<<endl;
  }
  listen(sock[0],2);
}


void runServer(ParamFitter& pf){
  struct sockaddr_un name_cl;
  int n = pf.getNbParameters();
  //  int cnt=0;
  Vector p1(n);
  float v;
  float ret;
  int neval=0,a=0;
  int vec_size = n*sizeof(float);
  sock[1]=-1;
  while(1){//main loop
    socklen_t size_cl;
    sock[1] = accept(sock[0],(struct sockaddr *) &name_cl,&size_cl);//should be blocking
    if(sock[1]>0){
      write(sock[1],&n,sizeof(int));//sends out the number of expected parameters
    }
    else{
      cerr<<"did not accept connection"<<" "<<strerror(errno)<<endl;
      return;// socket was deleted
    }
    if((a = read(sock[1],p1.GetArray(),vec_size))>0){//received a parameter vector
      if(a==vec_size){
	pf.setParameters(p1.Exp());
	int val = pf.eval(v);
	neval ++;
	if(!val){ret = NO_CONV_VAL;} // no convergence
	else{
	  if(val>0){
	    ret = v;
	  }
	  else{
	    ret = fabs(v);  
	  }
	}
	//send out results in sock[1]
	if(write(sock[1],&ret,sizeof(float))<0){
	  cerr<<"cannot write"<<endl;
	}
	close(sock[1]);
	sock[1]=-1;
      }
      else{//just a command
	if(*((int *)p1.GetArray())==0){ //close server
	  return;
	}
      }
    }
    else{
      cerr<<"cannot read parameters"<<endl;
      break;
    }
  }
}


void closeServer(){
  unlink(socketname);
  close(sock[1]);
  close(sock[0]);
}

void handle_sig(int s){
  closeServer();
  exit(0);
}

//#ifdef EVALSERVER_MAIN

int main(int argc, char *argv[]){
  if(argc<4){
    cerr<<"usage: evalserver <netfile> <expfile> <tarfile> <socketname> [<nbnets>]\n Beware: expfile format depends on <nbnets>: a list of files if n=0, a list of experiments (at least <nbnets>) if n>0)"<<endl;
    return 0;
  }
  int n=argc>4?atoi(argv[5]):0;
  signal (SIGINT,handle_sig);
  ParamFitter& pf = loadNetworks(argv[1],argv[2],argv[3],n); 
  strcpy(socketname,argv[4]);
  initServer(pf);
  sleep(5);
  runServer(pf);
  closeServer();
  delete &pf;
  return 0; 
}

//#endif
