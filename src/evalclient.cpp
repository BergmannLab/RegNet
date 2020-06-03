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

#include "Vector.h"
#include "evalclient.h"
#include "float.h"

using namespace std;
#define BUF_SIZE 100
#define NAME_LENGTH 30





/**
 * @returns the number of parameters
 */
int EvalClient::init(char *sockname){
 struct sockaddr_un name[2] ;
 sock[1] = socket(AF_LOCAL,SOCK_STREAM,0);
 
 name[0].sun_family = AF_FILE;
 strcpy(name[0].sun_path,sockname);
 size_t size = (offsetof (struct sockaddr_un, sun_path)
		+ strlen (name[0].sun_path) + 1);
 while(connect(sock[1],(struct sockaddr *) &(name[0]),size)<0){
 }
 int n=0;
 if(read(sock[1],&n,sizeof(int))>=0){
   //  cout<<"nb parameters "<<n<<endl;
 }
 else{
   cerr<<"cannot read number of parameters"<<endl;
 }
 return n;
}



void EvalClient::sendQuery(int n, float& val,Vector *par){
  if(n){
    Vector *pp;
    Vector p(n);
    if(!par){
      p.Random();
      pp = &p;
    }
    else{
      pp=par;
    }
    if(write(sock[1],pp->GetArray(),n*sizeof(float))<0){
      cerr<<"cannot write query"<<par<<endl;
    }
    else{
      if(read(sock[1],&val,sizeof(float))>=0){
      }
    }  
  }
  else{//close server
    int cmd=0;
    if(write(sock[1],&cmd,sizeof(int))<0){
      //      cout<<"cannot write"<<endl;
    }
  }
}


void EvalClient::sendNameQuery(int i, char *name){
  int cmd[2];
  cmd[0]= 1; //name query
  cmd[1]= i; //parameter index
  if(write(sock[1],cmd,2*sizeof(int))<0){
    cout<<"cannot send name query"<<endl;
  }
  else{
    int l=read(sock[1],name,sizeof(char)*NAME_LENGTH);
  } 
}
//obsolete
void EvalClient::sendOK(){
  int cmd=2;
  write(sock[1],&cmd,sizeof(int));
}

void EvalClient::closeClient(){
  close(sock[1]);
}

void EvalClient::closeServer(char *sockname){
  float v;
  init(sockname);
  sendQuery(0,v);
  close(sock[1]);
}

float EvalClient::eval(char *sockname, Vector *par){
  float v=FLT_MAX;
  int n=init(sockname);
  sendQuery(n,v,par);
  closeClient();
  return v;
}

void EvalClient::getParameterName(char *sockname,int i, char *name){
  int n=init(sockname);
  sendNameQuery(i,name);
  closeClient();
}


// int main(int argc, char *argv[]){
//   int n;
//   float val;
//   char c;
//   for(int i=0;i<1000;i++){
//     EvalClient client;
//     client.eval(argv[1],NULL);
//   }
//   EvalClient ec;
//   ec.init(argv[1]);
//   ec.sendQuery(0,val);
//   ec.closeClient();
// }

