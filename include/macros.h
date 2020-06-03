/** 
    
    Copyright (C) 2013 Micha Hersch, University of Lausanne
    email: micha.hersch@unil.ch

    This program is free software: you can redistribute it and/or modify
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


#ifndef __MACROS_H__
#define __MACROS_H__


#define RND(x)  ((x)*((double)rand())/((double)(RAND_MAX+1.0)))
#define RND_INT(x) ((rand())%(x))

#define sign(a) (a)>0?1:((a)==0?0:-1)
//#define max(a,b) ((a)>=(b)?(a):(b))
#define printvec(out,v) for(uint i=0;i<v.size();i++){out<<v[i]<<" ";}

/**
 * samples from a Gaussian pdf. Uses method described p.648 of Zwillinger (2003)
 * @param mean mean of normal distribution
 * @param var variance of normal distribution
 * @return a sample
 */
inline float norm_sample(float mean, float var){
 float v1,v2,r,x;
  r=2;
  while(r>1 || r==0.0){
    v1 = rand()/((float)RAND_MAX);
    v1 = 2*v1-1;
    v2 = rand()/((float)RAND_MAX);
    v2 = 2*v2-1; 
   r= v1*v1+v2*v2;
  }
  x=v1*sqrtf(-2*log(r)/r);
  return x*sqrtf(var)+mean;
}


/* inline ostream& operator<<(std::ostream& out, */
/* 			    const vector<T>& v){ */
/*   for(int i=0;i<v.size();i++){out<<v[i]<<" ";} */
/*   return out; */
/* } */

#endif
