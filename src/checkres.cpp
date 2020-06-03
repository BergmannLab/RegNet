#include "paramfitter.h"

using namespace std;

void usage(){
  cout<<"./checkres <network> <experimentslist> <parameters> [<pos_score> <allstates> [<map> <targetfile>]]"<<endl;
  cout<<"<pos_score> if one, discards negative score (not used)"<<endl;
  cout<<"<allstates> if one, outputs all network states (should be 0)"<<endl;
  cout<<"<map> only prints out the result of the optimal parameter vector (not used)"<<endl;
}

int main(int argc,const char *argv[]){
  if(argc<=2){usage();return 1;}
  srand(time(NULL));
  ParamFitter pf(argv[1],argv[2]);
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
  return 0;
}
