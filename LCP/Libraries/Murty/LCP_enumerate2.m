function [X, DegenSoln, IsolatedSoln, NonIsolated, NonIsolatedGrp, ...
          ColCompetent, DegenMatrix, n_DegenMatrix] = ...
         LCP_enumerate2(M, q, EPS)
%
% solve the linear complementarity problem (LCP) using an enumerative
% approach that searches all 2^n possibilities to find all solutions
% of the LCP:
%   find "x" such that  M*x >= 0; x >=0; x'*(M*x) = 0
%
% Input: M = square matrix, nxn
%        q = column vector, nx1
%
% Output: X =            matrix with columns x that are solutions of the LCP
%                        Note that each column is unique, as redundant columns
%                        are removed.
%         DegenSoln =    row vector of 0s and 1s. A 0 means that the
%                        corresponding solution in X is not a degenerate
%                        solution. A 1 means that the corresponding solution 
%                        is degenerate.
%         IsolatedSoln = row vector of 0s and 1. A 1 means that the
%                        corresponding solution in X is an isolated
%                        solution.  A 2 means that the soltions is not
%                        isolated.
%         ColCompetent = 1 or 0, whether matrix M is column competent or not
%         DegenMatrix =  if the matrix is non-degenerate, then the empty
%                        matrix [] is assigned.  If M is degenerate, then
%                        DegenMatrix is a matrix with n rows.  
%                        The columns of DegenMatrix indicate the index
%                        sets for which the principal minor is zero.
%                        The index set is given by the rows of the
%                        1's in a column.  A matrix M can be degenerate
%                        with multiple prinicpal minors that are zero.
%                        Hence, possible multiple columns in DegenMatrix.
%         NonIsolated =  matrix of directions in which non-isolated solutions
%                        (like a polygonal edge) emanate from the solutions 
%                        in matrix X.  NonIsolated has n rows, and one column
%                        for each solution.  A column is filled with zeros, 
%                        if the solution is isolated.  Otherwise, the column
%                        contains a unit vector of the direction of the null
%                        space of the complemenary matrix, giving the
%                        direction "z" of solutions near to "x": x + eps*z.
%                        Note that the code only handles complemenary 
%                        matrices with a rank of n-1, not less.
%
% Limitation: for degenerate matrices, we currently are limited to cases
% in which complementary matrices are of rank n-1.
%
% Method: we use the method described my Murty and Yu, page 8:
%   Katta G. Murty and Feng-Tien Yu, Linear Complementarity: Linear and
%   Nonlinear Programming, Internet edition, 1997.

%   Murty, Katta G. "On the number of solutions to the complementarity 
%   problem and spanning properties of complementary cones." Linear 
%   Algebra and its Applications 5.1 (1972): 65-108.
%
% Limitation: the implementation time is of order 2^n, so large problems 
% (say, n>15) are intractable.  If you only need a single solution of a
% multi-solution LCP, or if you are sure that the LCP has only a single
% solution (for example, if M is a P-matrix), then there are several other
% Matlab/Octave codes that will solve the problem more efficiently than
% this function.  I suggest trying qp() or quadprog() in Octave, or 
% quadprog in Matlab.
%
  X = [];
  DegenSoln = [];
  DegenMatrix = [];
  IsolatedSoln = [];
  NonIsolated = [];
  ColCompetent = [];
%
% these tolerances are used to determine whether a matrix is nearly singular
% EPS1 = 1e-10; EPS1 = 10*eps;
%
% this tolerance is used to determine whether complementary values are
% positive
% EPS2 = 1e-12; EPS2 = 10*eps;
%
% we set a limit on how large the matrix DegenMatrix can grow
  Max_Recorded_DegenMatrix = 100;
%
% number of degenerate matrices embedded in M.  Note that the matrix
% DegenMatrix is limited in size to Max_Recorded_DegenMatrix columns
  n_DegenMatrix = 0;
%
  n = length(q);
%
% initialize DegenMatrix
  DegenMatrix = zeros(n, Max_Recorded_DegenMatrix);
%
  NonIsolatedGrp = [];
%
  N = dec2bin(0:(2^n-1),n);
  Ones = N=='1';
%
% matrix is column competent unless proven otherwise
  ColCompetent = 1;
%
% check whether M_lambda is zero, and M_hat and q_hat are empty
  if n==0
    disp('M_hat and q_hat are empty.  Invalid LCP in LCP_enumerate2.m')
  else
%
%   search through all 2^n possible solutions
    for r = 1:2^n
%     construct the complementary matrix "Ar"
      Ar = eye(n);
      Indices = find(Ones(r,:));
      Ar(:,Indices) = -M(:,Indices);
%
%     check whether the principal submatrix is ill-conditioned (nearly singular).
%     Here, we take advantage of Octave able to compute both the inverse
%     (nvAr) and the condition number of matrix Ar.  Matlab might not have this
%     ability.  If Matlab is used, you might need to use the cond() function
%     to find the condition number and the inv() function to find the inverse.
      [invAr, RCOND] = inv(Ar);
%
      if RCOND >= EPS
%       the RCOND value is the reciprocal condition number.  In this case,
%       matrix Ar is non-singular and the corresponding principal submatrix
%       has a non-zero determinant

%       solve matrix equation for a trial vector
        y = invAr*q;
%
%       a solution exists when all elements of "y" are non-negative
        if all(y >= -EPS)
%
%         the zero and non-negative "x" values
          x = zeros(n,1);
          x(Indices) = y(Indices);
%
%         append the solution "x" to the other solutions
          X = [X x];
%
%         the solution is degenerate if any of the "y" values are zero
          DegenSoln = [DegenSoln any(abs(y)<=EPS)];
%
%         because the det(Ar) is not zero, the solution is isolated
          IsolatedSoln = [IsolatedSoln 1];
%
%         an isolated solution
          NonIsolated = [NonIsolated, zeros(n,1)];
%
%         non-isolated group
          NonIsolatedGrp = [NonIsolatedGrp, 0];
        end
      else
%       the corresponding principal submatrix of M has a zero determinant. Output 
%       the indices of the principal submatrix in "DegenMatrix" 
        Degen = zeros(n,1);
        Degen(Indices) = 1;
        n_DegenMatrix = n_DegenMatrix + 1;
        if n_DegenMatrix <= Max_Recorded_DegenMatrix
          DegenMatrix(:,n_DegenMatrix) = Degen;
        end
%
%       determine whether matrix M is column competent using Theorem 2.5 of 
%       Song Xu, "On local w-uniqueness of solutions to linear complementarity
%       problem," Linear Algebra and its Applications, 1999, 290:23-29.
        ColCompetent = ColCompetent && rank(M(:,Indices),EPS) < length(Indices);
%
%       use the Moore-Penrose pseudo-inverse (employing singular value
%       decomposiiton).  This will give a projection of q onto the column
%       space of Ar.  The result "z" might not be a solution of Ar*z = q
        z = pinv(Ar) * q;
%
%       check whether the projection results actually satisfy Ar*z = q
        if all(abs(Ar*z - q) <= 10*EPS)
%         "z" is actually a solution of Ar*z = q, but it is not a unique
%         solution. We must now determine whether the vector "q" allows for 
%         a set of solutions "y" to satisfy y>=0
%
%         because Ar is singular, z is not a unique solution, and the set 
%         of all solutions is y plus the range of the null space of Ar
          NullAr = null(Ar, EPS);
%
%         rank of the null space of Ar
          RankNullAr = size(NullAr, 2);
%
%         search the hyperplane of z + NullSpace for intercepts with the 
%         the hyperplanes of the lambda axes
          Combinations = nchoosek(1:n, RankNullAr);
%
          for Combination = 1:size(Combinations, 1)
            iknown_y = Combinations(Combination,:);
            iunknown_y = setdiff(1:n,iknown_y);
%
            y = zeros(n,1);
            c = - NullAr(iknown_y,:) \ z(iknown_y);
            y(iunknown_y) = z(iunknown_y) + NullAr(iunknown_y,:)*c;
%
            if all(y  > -EPS)
%             we have found an intercept y = z + NullSpace that is a solution
%
%             the zero and non-negative "x" values
              x = zeros(n,1);
              x(Indices) = y(Indices);
%
%             append the solution "x" to the other solutions
              X = [X x];
%
%             the solution is degenerate if any of the "y" values are zero
              DegenSoln = [DegenSoln any(abs(y)<=EPS)];
%
%             because the det(Ar) is not zero, the solution is isolated
              IsolatedSoln = [IsolatedSoln 0];
%
%             a non-isolated solution
              NonIsolated = [NonIsolated, ismember([1:n]', iunknown_y')];
%
%             non-isolated group
              NonIsolatedGrp = [NonIsolatedGrp, r];
            end
          end
%
        else
%         q is not within the range of Ar, so this "r" does not result in a
%         solution
        end
      end  % if 1/RCOND <= COND
    end  % for r = 1:2^n
%
    DegenMatrix = DegenMatrix(:,1:min(Max_Recorded_DegenMatrix,n_DegenMatrix));
  end
%
