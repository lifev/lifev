/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include<cmath>
#include<iostream>
#include<algorithm>

//! Processs markers in the form of a long int
/*! The marker flag is given by a long int of the form
\verbatim
     xxxyyyzzz

\endverbatim
where xxx is a group of n digits which may be subsequently reanalysed
une by one
*/
typedef long  MarkerFlag;


class process_marker_flag
{
public:
  explicit process_marker_flag(unsigned int ndg, unsigned int ng);
  process_marker_flag();

  void set_number_digits_per_group(unsigned int const ndg);
  void set_number_of_groups(unsigned int const ng);
  inline unsigned int number_digits_per_group() const;
  inline  unsigned int  number_of_groups() const;
  unsigned int extract_group(MarkerFlag const & mf, unsigned int const & group) const;
  inline  unsigned int
  extract_digit_in_group(unsigned int const &  group , unsigned int const & digit) const;
private:
  void set_factor();
  unsigned int my_ndg;
  unsigned int my_ng;
  unsigned int factor;
  unsigned int factor10;
};


process_marker_flag::process_marker_flag(unsigned int ndg, unsigned int ng):
  my_ndg(ndg),my_ng(ng){
  set_factor();
}

process_marker_flag::process_marker_flag():  my_ndg(0),my_ng(0){
  set_factor();
}

void
process_marker_flag::set_factor()
{
    factor = std::pow(10.0,std::max((int)my_ndg, (int)0));
    //cout << "FACTOR: "<<factor<<std::endl;
}

void process_marker_flag::set_number_digits_per_group(unsigned int const ndg){
  my_ndg=ndg;
  set_factor();
}

void
process_marker_flag::set_number_of_groups(unsigned int const ng){
  my_ng=ng;
}

unsigned int process_marker_flag::number_digits_per_group() const{
  return my_ndg;
}


unsigned int  process_marker_flag::number_of_groups() const{
  return my_ng;
}


unsigned int process_marker_flag::extract_group(MarkerFlag const &  mf, unsigned int const & group) const
{
    unsigned long ff = std::pow((double)factor,std::max(int(my_ng- group),int(0))); //group numbering starts from 1
    //std::cout << "FF_GROUP: "<<ff<<std::endl;
    return (mf/ff) % factor;
}

unsigned int
process_marker_flag::extract_digit_in_group(unsigned int const &  group , unsigned int  const &  digit) const{
  unsigned int ff= std::pow(10.0,std::max(int(my_ndg-digit),int(0)));
  //std::cout << "FF_DIGIT: "<<ff<<std::endl;
  return (group/ff) % 10;
}

int
main(){
    MarkerFlag flag(1);
    unsigned int groups=4;
    unsigned int group;
    unsigned int ndg=2;
    unsigned int digit;

    process_marker_flag proc(ndg,groups);

    do{
        std::cout<< " Gimme MarkerFlag (digits) "<<ndg*groups<<std::endl;
        std::cin >> flag;
        std::cout<<flag<<std::endl;
        do{
            std::cout<< " Gimme Group"<<std::endl;
            std::cin >> group;
            digit=group;
            std::cout<<proc.extract_group(flag,group)<<std::endl;
            do{
                std::cout<< " Gimme Digit"<<std::endl;
                std::cin >> digit;
                std::cout<<proc.extract_digit_in_group(proc.extract_group(flag,group),digit)<<std::endl;
            }while(digit!=0);
        }while(group !=0);
    }while(flag!=0);

}
