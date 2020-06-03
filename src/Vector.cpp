/* 
 * Copyright (C) 2008 Eric Sauser, Micha Hersch, EPFL
 * RobotCub Consortium, European Commission FP6 Project IST-004370
 * email:   micha.hersch@robotcub.org
 * website: www.robotcub.org
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * A copy of the license can be found at
 * http://www.robotcub.org/icub/license/gpl.txt
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */

#include "Vector.h"

float Vector::undef = 0.0f;


std::ostream& operator<<(std::ostream& out, const Vector& v){
  for(unsigned int i=0;i<v.Size()-1;i++){
    out<<v.At(i)<<" ";
  }
  out<<v.At(v.Size()-1);
  return out;
}
