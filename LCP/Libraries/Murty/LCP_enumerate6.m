function [ X, Xr, DegenSoln, IsolatedSoln, NonIsolated, NonIsolatedGrp, ...
           NonIsolatedDim, ColCompetent, DegenMatrix, n_DegenMatrix] = ...
         LCP_enumerate6(M, q, EPS1, EPS2)
%
% solve the linear complementarity problem (LCP) using an enumerative
% approach that searches all 2^n possibilities to find all solutions
% of the LCP:
%   find "x" such that  M*x >= 0; x >=0; x'*(M*x) = 0
%
% Input: M = square matrix, nxn
%        q = column vector, nx1
%        EPS1 = tolerance applied to the rank() and RCOND values
%        EPS2 = tolerance applied to whether y() is non-negative
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
  Xr = [];
  DegenSoln = [];
  DegenMatrix = [];
  IsolatedSoln = [];
  NonIsolated = [];
  ColCompetent = [];
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
  NonIsolatedDim = [];
%
% enumerate the 2^n possible combinations of sliding/nonsliding contacts
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
%   maximum singular value of M
    Max_S = max(abs(svd(M)));
%
%   search through all 2^n possible solutions
    for r = 1:2^n
%     construct the complementary matrix "A"
      A = eye(n);
      J = find(Ones(r,:));  % the set J_C or (I_2 union I_3) \ I_1
      A(:,J) = -M(:,J);
%
%     check whether the principal submatrix is ill-conditioned (nearly singular).
%     Here, we take advantage of Octave able to compute both the inverse
%     (invA) and the condition number of matrix A.  Matlab might not have this
%     ability.  If Matlab is used, you might need to use the cond() function
%     to find the condition number and the inv() function to find the inverse.
      [invA, RCOND] = inv(A);
%
      if RCOND >= EPS1
%       the RCOND value is the reciprocal condition number.  In this case,
%       matrix A is non-singular and the corresponding principal submatrix
%       has a non-zero determinant

%       solve matrix equation for a trial vector
        y = invA*q;
%
%       a solution exists when all elements of "y" are non-negative
        if all(y >= -EPS2)
%
%         the zero and non-negative "x" values
          x = zeros(n,1);
          x(J) = y(J);
%
%         append the solution "x" to the other solutions
          X = [X x];
          Xr = [Xr r];
%
%         the solution is degenerate if any of the "y" values are zero
          DegenSoln = [DegenSoln any(abs(y)<=EPS2)];
%
%         because the det(A) is not zero, the solution is isolated
          IsolatedSoln = [IsolatedSoln 1];
%
%         an isolated solution
          NonIsolated = [NonIsolated, zeros(n,1)];
%
%         non-isolated group
          NonIsolatedGrp = [NonIsolatedGrp, 0];
%
%         dimension of the solution set (isolated point of dimension 0)
          NonIsolatedDim = [NonIsolatedDim, 0];
        end
      else
%       the corresponding principal submatrix of M has a zero determinant. Output 
%       the indices of the principal submatrix in "DegenMatrix" 
        Degen = zeros(n,1);
        Degen(J) = 1;
        n_DegenMatrix = n_DegenMatrix + 1;
        if n_DegenMatrix <= Max_Recorded_DegenMatrix
          DegenMatrix(:,n_DegenMatrix) = Degen;
        end
%
%       determine whether matrix M is column competent using Theorem 2.5 of 
%       Song Xu, "On local w-uniqueness of solutions to linear complementarity
%       problem," Linear Algebra and its Applications, 1999, 290:23-29.
        ColCompetent = ColCompetent && rank(M(:,J),EPS1*Max_S) < length(J);
%
        if rank(A, EPS1*Max_S) == rank([A, q], EPS1*Max_S)
%         although A is singular, the equation A*y=q has solutions.
%         The non-negative solutions will form a convex set (a convex
%         component) of non-isolated solutions.
%
%         rank of A
          rankA = rank(A, EPS1*Max_S);
%
%         dimension of the convex component of non-isolated solutions
          dimNonIsol = n - rankA;
%
%         the complement of J, with J_C representing the elements of x
%         that are assumed to be zero, with this complementary matrix A
          J_C = setdiff(1:n, J);
%
%         find the extreme-points of the convex set.  These extreme-points
%         are at the intersection of n + dimNonIsol planes:
%           1) the x(J_C) >= 0 planes, as presumed with the A matrix
%           2) the (M*x + q)(J_C) >=0 planes, as presumed with the A matrix
%           3) a combination of dimNonIsol of the following planes
%              a) x(i) = 0, for those i \in J of the combination
%              b) (M*x + q)(i) = 0, for those i \in J_C of the combination
%
%         create a set of dimNonIsol elements taken from the set 1:n. Each row
%         in "Ks" is a combination of dimNonIsol elements in 1:n
          Ks = nchoosek(1:n, dimNonIsol);
%
%         solutions with this complementary matrix A
          X_of_A = [];
%
          for k = 1:size(Ks,1)
%           a row in Ks that gives the combination of dimNonIsol elements
            K = Ks(k,:);
%
%           the complementary matrix A is augmented with additional
%           conditions, given in (3) above
            A_aug = zeros(dimNonIsol, n);
            q_aug = zeros(dimNonIsol, 1);
%
            for i = 1:dimNonIsol
              if ismember(K(i), J)
%               although previously assumed non-negative, x(i) is now zero.
%               Note that x(J_C) is assumed zero.
                A_aug(i,:) = eye(n)(K(i),:);
                q_aug(i) = 0;
              else
%               K(i) is an element of J_C, (M*x + q)(i) was previously assumed
%               non-negative, but is now zero
                A_aug(i,J) = -M(K(i),J);
                q_aug(i) = q(K(i));
              end
            end
%
%           append these conditions to A, with zeros to the q vector
            B = [A; A_aug];
            b = [q; q_aug];
%
%           use the Moore-Penrose pseudo-inverse (employing singular value
%           decomposiiton)
            y = B \ b;
%
%           check whether y is an actual solution
            if sum(abs(b - B*y)) < n*EPS2
%
%             a possible solution
              x = zeros(n,1);
              x(J) = y(J);
%
              w = M*x + q;
%
%             check the non-negative and complementarity conditions
              if   min(min(w),min(x)) > -EPS2 ...
                && isempty(intersect(find(w>EPS2),find(x>EPS2)))
%
%               append the solution "x" to the other solutions
                X = [X, x];
                Xr = [Xr, r];
                X_of_A = [X_of_A, x];
%
%               the solution is degenerate if any of the "y" values are zero
                DegenSoln = [DegenSoln any(abs(y)<=EPS2)];
%
%               because the det(A) is zero, the solution is isolated
                IsolatedSoln = [IsolatedSoln 0];
%
%               a non-isolated solution
                NonIsolated = [NonIsolated, ismember([1:n]', K')];
%
%               non-isolated group
                NonIsolatedGrp = [NonIsolatedGrp, r];
%
%               dimension of the solution set
%               NonIsolatedDim = [NonIsolatedDim, dimNonIsol];
              end
            end % if sum(abs(b - B*y)) < n*EPS2
          end % for k = 1:size(Ks,1)
%
%         find the dimension of the non-isolated solutions, dim(P). Begin
%         by finding the index types of the solutions that were found with
%         this complementary matrix A
          if ~isempty(X_of_A)
            Index_Type = zeros(size(X_of_A));
            for k = 1:size(X_of_A,2)
              x = X_of_A(:,k);
              w = M*x + q;
              Index_Type(find(x < EPS2 & w > EPS2),k) = 1;
              Index_Type(find(x > EPS2 & w < EPS2),k) = 2;
              Index_Type(find(x < EPS2 & w < EPS2),k) = 3;
            end
%
%           the index set C(P), the union of indices with Index_Type=2
            C_P = find(any(Index_Type==2, 2));
%
%           the index set E(P), the union of indices with Index_Type=2 or 3
            E_P = find(any(Index_Type==2 | Index_Type==3, 2));
%
%           the dimension of the non-isolated solutions
            dim_P = length(C_P) - rank(M(E_P, C_P), EPS1*Max_S);
%
%           place the same dim(P) into the vector NonIsolatedDim
            NonIsolatedDim = [NonIsolatedDim, dim_P*ones(1, size(X_of_A,2))];
          end
        else
%         q is not within the range of A, so this "r" does not result in a
%         solution
        end % if rank(A, EPS1*Max_S) == rank([A, q], EPS1*Max_S)
      end  % if 1/RCOND <= COND
    end  % for r = 1:2^n
%
    DegenMatrix = DegenMatrix(:,1:min(Max_Recorded_DegenMatrix,n_DegenMatrix));
    if isempty(DegenMatrix)
      DegenMatrix = 0;
    end
  end
%
