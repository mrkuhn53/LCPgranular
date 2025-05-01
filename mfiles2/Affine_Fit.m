function [Affine_Matrix] = Affine_Fit(x);
%
% this function determines a matrix that finds the "best fit" of an
% affine field to the actual particle motions
%   dx \approx Affine_Matrix(1:3) + Affine_Matrix(4:7)*x
% or
%   dx_1_p = A_1 +  A_4*x_1_p +  A_5*x_2_p +  A_6*x_3_p
%   dx_2_p = A_2 +  A_7*x_1_p +  A_8*x_2_p +  A_9*x_3_p
%   dx_3_p = A_3 + A_10*x_1_p + A_11*x_2_p + A_12*x_3_p
%
  N = size(x,1);
%
% particle displacements
% Affine_Matrix = zeros(3*N,3*N);
  Affine_Matrix = zeros(3*N,12);
%
  Affine_Matrix(1:3:3*N,1) = 1;
  Affine_Matrix(2:3:3*N,2) = 1;
  Affine_Matrix(3:3:3*N,3) = 1;
%
  Affine_Matrix(1:3:3*N,4:6)   = x;
  Affine_Matrix(2:3:3*N,7:9)   = x;
  Affine_Matrix(3:3:3*N,10:12) = x;
