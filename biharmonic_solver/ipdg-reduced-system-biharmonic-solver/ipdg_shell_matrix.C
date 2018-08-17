// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local includes
#include "ipdg_shell_matrix.h"
#include "libmesh/numeric_vector.h"
// #include "libmesh/petsc_vector.h"

namespace libMesh
{

template <typename T>
void IPDGShellMatrix<T>::vector_mult (NumericVector<T> & dest,
                                      const NumericVector<T> & arg) const
{
  dest.zero();
  this->vector_mult_add(dest,arg);
}


template <typename T>
void IPDGShellMatrix<T>::vector_mult_add (NumericVector<T> & dest,
                                          const NumericVector<T> & arg) const
{
    // local variables
    UniquePtr <NumericVector<T> > dest1 = dest.clone();
    UniquePtr <NumericVector<T> > dest2 = dest.clone();
    UniquePtr <NumericVector<T> > arg_local = dest.clone();
   
    *arg_local = arg;
    (*dest1).zero();
    AT.vector_mult_add(*dest1, *arg_local); // compute A^T x
    *arg_local = *dest1; (*dest1).zero();
    Minv.vector_mult_add(*dest1, *arg_local); // compute Minv A^T x
    *arg_local = *dest1; (*dest1).zero();
    A.vector_mult_add(*dest1, *arg_local); // compute A Minv A^T x
    
    (*dest2).zero();
    J.vector_mult_add(*dest2, arg); // compute J x
    
    // compute dest + A Minv A^T x + J x
   *dest1 += dest;
    dest += *dest1;
    dest += *dest2;
    
}

template <typename T>
void IPDGShellMatrix<T>::get_diagonal(NumericVector<T> & dest) const
{
    
    UniquePtr<NumericVector<T> > unit_vector = dest.clone();
    UniquePtr<NumericVector<T> > foo = dest.clone();
    dest.zero();
    (*foo).zero();
    (*unit_vector).zero();
    
    for(int ii = 0; ii < A.m(); ii++)
    {
        (*unit_vector).set(ii,1.0);
        this->vector_mult(*foo, *unit_vector);         // pick out the iith column
        dest.set( ii, (*foo).el(ii) );                 // take the diagonal element from that column
        (*unit_vector).set(ii, 0.0);
        
    }
    
}



//------------------------------------------------------------------
// Explicit instantiations
template class IPDGShellMatrix<Number>;

} // namespace libMesh
