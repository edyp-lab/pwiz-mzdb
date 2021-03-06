//
// $Id: Random.hpp 1191 2009-08-14 19:33:05Z chambm $
//
//
// Darren Kessner <darren@proteowizard.org>
//
// Copyright 2009 Spielberg Family Center for Applied Proteomics
//   Cedars Sinai Medical Center, Los Angeles, California  90048
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//


#ifndef _RANDOM_HPP_
#define _RANDOM_HPP_


namespace pwiz {
namespace math {


class Random
{
    public:

    static void initialize();               // initialize the randomizer
    static double real(double a, double b); // random double in [a,b)
    static int integer(int a, int b);       // random int in [a,b)
    static double gaussian(double sd);      // random gaussian w/deviation sd
};


} // namespace math
} // namespace pwiz


#endif // _RANDOM_HPP_

