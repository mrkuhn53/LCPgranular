function [X, DegSol, DegMatrix, IsolSol] = LCP_enumerate_old(M,q)
%
% solve the linear complementarity problem (LCP) using an enumerative
% approach that searches all 2^n possibilities
%
  X = [];
  DegSol = [];
  DegMatrix = [];
  IsolSol = [];
%
  EPS = 1e-40;
%
  n = length(q);
%
  N = dec2bin(0:(2^n-1),2);
  Ones = N=='1';
%
% search through all 2^n possible solutions
  for r = 1:2^n
%   construct the complementary matrix "Ar"
    Ar = eye(n);
    Indices = find(Ones(r,:));
    Ar(:,Indices) = -M(:,Indices);
%
%   check whether the principal submatrix has zero determinant
%   detAr = det(Ar);
    if det(Ar) ~= 0
%     y = Ar\q;
      [invAr,RCOND] = inv(Ar);
      y = invAr*q;
%
%     a solution exists when all elements of "y" are non-negative
      if all(y>=0)
%
%       the zero and non-negative "x" values
        x = zeros(n,1);
        x(Indices) = y(Indices);
        X = [X x];
%
%       the solution is degenerate if any of the "y" values are zero
        DegSol = [DegSol any(y==0)];
%
%       because the det(Ar) is not zero, the solution is isolated
        IsolSol = [IsolSol 1];
%
%       degenerate solutions occurs as pairs.  Remove redundant ones
        [Y,I,J] = unique(X','rows');
        X = X(:,I);
        DegSol = DegSol(I);
        IsolSol = IsolSol(I);
      end
    else
%     this principal submatrix has a zero determinant. Output the
%     indices of the principal submatrix in "DegMatrix" 
      Deg = zeros(n,1);
      Deg(Indices) = 1;
      DegMatrix = [DegMatrix Deg];
%
      if rank(Ar,EPS)==(n-1)
        y = Ar\q;
        if all(y>=0) || all(y<=0)
          x = zeros(n,1);
          x(Indices) = y(Indices);
          X = [X x];
          DegSol = [DegSol any(y==0)];
          IsolSol = [IsolSol 0];
          [Y,I,J] = unique(X','rows');
          X = X(:,I);
          DegSol = DegSol(I);
          IsolSol = IsolSol(I);
        end
      end
    end
  end
