function [Lambda, NegEigenGrp, NegEigenValue, NegEigenVector] = ...
         LCP_Instability(P_L, H_x, H_lambda, Q_x, Q_lambda, EPS)
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
  Lambda = [];
  NegEigenGrp = [];
  NegEigenValue = [];
  NegEigenVector = [];
%
% number of active contacts, M_lambda
  n = size(H_lambda, 2);
%
% we set a limit on how large the matrix DegenMatrix can grow
  Max_Recorded_DegenMatrix = 100*n;
%
% the LCP matrix
  M = eye(n);
%
% inverse of the H_lambda_lambda matrix
  Q_lambda_inv = inv(Q_lambda);
%
  N = dec2bin(0:(2^n-1),n);
  Ones = N=='1';
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
      J = find(Ones(r,:));
      Ar(:,J) = -M(:,J);
%
      if mod(r,5000)==0; disp(sprintf('Combination %6i of %6i', r, 2^n)); end
%
%     intermediate matrix, for finding the "lambda" contribution to W_2
      B = zeros(n);
      B(J,:) = Ar(J,:) * Q_lambda_inv;
%
%     2nd-order work for this complementary matrix. Note that inv(Ar) = Ar
      W_2_matrix = ...
        P_L' ...
        * (H_x' +  H_lambda * B * Q_x)...
        * P_L;
%
%     symmetric part of W_2_matrix
      W_2_matrix_sym = 0.5 * (W_2_matrix + W_2_matrix');
%
%     find the negative eigenvalues and the corresponding eigenvectors
      [V, LAMBDA] = eig(balance(W_2_matrix_sym,'noperm'));
%  
%     convert LAMBDA from a diagonal matrix to a vector
      LAMBDA = diag(LAMBDA);
%  
%     we are concerned with finding negative eigenvalues, which correspond
%     to possible instability modes.  Ignore the positive eigenvalues and
%     their eigenvectors
      if any(LAMBDA < -10*EPS)
%
%       isolate the eigenvectors that have negative eigenvalues.
        Neg_Eigenvalues = find(LAMBDA < -10*EPS);
        Neg_Eigenvectors = V(:,Neg_Eigenvalues);
%
        Number_Neg_Eigenvalues = length(Neg_Eigenvalues);
%
        if Number_Neg_Eigenvalues > 0
%         eigenvectors are potential solutions.  However, they must 
%         also satisfy the LCP inequalities
%
%         note that inf(Ar) = Ar
          y = Ar * Q_lambda_inv * Q_x * P_L * V(:,Neg_Eigenvalues);
%
%         check the LCP inequalities
          for i =1:Number_Neg_Eigenvalues
            if all(y(:,i) >= -1*EPS) || all(y(:,i) <= 1*EPS)
%             all of the LCP criteria are satisfied for this eigenvector
%
              if all(y(:,i) >= -1*EPS)
                Sign = 1;
              else
                Sign = -1;
              end
%
%             the zero and non-negative "lambda" values
              lambda = zeros(n,1);
              lambda(J) = Sign*y(J,i);
%
%             append the solution "lambda" to the other solutions
              Lambda = [Lambda lambda];
              NegEigenGrp = [NegEigenGrp, r];
              NegEigenValue = [NegEigenValue, LAMBDA(Neg_Eigenvalues(i))];
              NegEigenVector = [NegEigenVector, Sign*V(:,Neg_Eigenvalues(i))];
            end
          end % i =1:Number_Neg_Eigenvalues
        end % if Number_Neg_Eigenvalues > 0
      end % if any(LAMBDA < -10*EPS)
    end % for r = 1:2^n
  end % if ~(n==0)
