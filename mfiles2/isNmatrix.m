function [z] = isNmatrix(M)
%
% determine whether M is a P-matrix
%
  z = 1;
%
% these tolerances are used to determine whether a determinate is nearly 
% positive 
  EPS = 10*eps;
%
  n = size(M,1);
%
  N = dec2bin(0:(2^n-1),n);
  Ones = N=='1';
%
% search through all 2^n possible solutions
  for r = 1:2^n
%   construct the complementary matrix "Ar"
    Indices = find(Ones(r,:));
    A = M(Indices,Indices);
    z = z && det(A) < -EPS;
    if ~z; break; end
  end  % for r = 1:2^n
%
