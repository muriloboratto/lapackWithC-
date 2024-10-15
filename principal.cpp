/*%****************************************************************************80
%  Code: 
%   principal.cpp 
%
%  Purpose:
%   Implements the LAPACK routine in C++ (OO).
%
%  Modified:
%   Oct 15 2024 10:15 
%
%  Author:
%    Murilo Do Carmo Boratto [muriloboratto@gmail.com]
%
%  How to Compile:
%    make 
%
%  How to Execute: 
%   ./principal 10           
%   
%****************************************************************************80*/

#include "libreria.h"

int main(int argc, char **argv)
{
    int rows = atoi(argv[1]); // size problem
  
    LAPACK Object(rows);
  
        Object.dgemm();   /*multiply matrix*/  
      
        Object.dsyevd();  /*eigenvalues Divide and Conquer type 1 - Dense Matrix */
      
        Object.dstevd();  /*eigenvalues Divide and Conquer type 2 - Tridiagonal Symmetric*/
     
    return 0;

}/*main*/

