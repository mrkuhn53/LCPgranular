%
% |+|+|+|+|+|+|+|+|  Revision 1   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function z = Cross_Product(a)
%
% this function creates a matrix that can be used to compute the cross-product
% of two vectors.  The following cross-product,
%   a X b
% where "a" is a vector (either a row vector or a column vector)
% and "b" is a column vector (note, "b" must be a column vector) is
% produced by the the product
%   Cross_Matrix(a)*b
%
% Here is the matrix:
  z = [  0   -a(3)  a(2); ...
       a(3)   0   -a(1); ...
      -a(2)  a(1)   0  ];
