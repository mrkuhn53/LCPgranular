function [ Lambda, NegEigenGrp, NegEigenValue, NegEigenVector] = ...
         LCP_Instability_alt( ...
           P_L, H_xx_BD_inv, ...
           H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, EPS);
%        LCP_Instability2(P_L, H_x, H_lambda, Q_x, Q_lambda, EPS)
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
% Limitation: the implementation time is of order 2^n, so large problems 
% (say, n>15) are intractable.  If you only need a single solution of a
% multi-solution LCP, or if you are sure that the LCP has only a single
% solution (for example, if M is a P-matrix), then there are several other
% Matlab/Octave codes that will solve the problem more efficiently than
% this function.  I suggest trying qp() or quadprog() in Octave, or 
% quadprog in Matlab.
%
  lambda = [];
  Lambda = [];
  NegEigenGrp = [];
  NegEigenValue = [];
  NegEigenVector = [];
%
% number of active contacts, M_lambda
  n = size(H_lambda_lambda, 2);
%
% we set a limit on how large the matrix DegenMatrix can grow
  Max_Recorded_DegenMatrix = 100*n;
%
  N = dec2bin(0:(2^n-1),n);
  Ones = N=='1';
%
% check whether M_lambda is zero, and M_hat and q_hat are empty
  if n==0
    disp('M_hat and q_hat are empty.  Invalid LCP in LCP_enumerate2.m')
  else
%   inverse of the H_lambda_lambda matrix
    H_lambda_lambda_inv = inv(H_lambda_lambda);
%
%   is multiplied by lambda to compute z
    Q = H_xx_BD_inv * H_x_lambda;
%
%   inner part of quadratic form for second-order work, W_2
    R = Q' * (H_x_x' * Q - H_x_lambda);
%
%   matrices for computing q_dot
    S = H_lambda_lambda_inv * H_lambda_x * Q;
%
%   quadratic form for second-order work, W_2
    T = S' * R * S;
%
%   list of indices
    I = [1:n];
%
%   XXXX the case of r=1 here XXXXX
%
%   search through all 2^n possible solutions
    for r = 2:2^n
%     list of sliding contacts
      J = find(Ones(r,:));
      Jc = setdiff(I, J);
%
      if mod(r,5000)==0; disp(sprintf('Combination %6i of %6i', r, 2^n)); end
%
%     eliminate the Jc rows and columns from R, since lambda(Jc) = 0
      T_eff = T(J,J);
%
%     symmetric part of W_2_matrix
      W_2_matrix_sym = balance(0.5 * (T_eff + T_eff'), 'noperm');
%
%     try Cholesky factorization to determine whether the matrix is positive
%     definite
      [Rchol, P] = chol(W_2_matrix_sym);
%
      if P
%       if P is not zero, then W_2_matrix_sym is not positive definite, and
%       we can proceed to find the negative eigenvalues
%
%       find the negative eigenvalues and the corresponding eigenvectors
        [V, LAMBDA] = eig(W_2_matrix_sym);
%  
%       convert LAMBDA from a diagonal matrix to a vector
        LAMBDA = diag(LAMBDA);
      else
%       P is zero, so W_2_matrix_sym is positive definite.  Force a positive
%       eigenvalue, so the the next statements are skipped
        LAMBDA = 1;
      end
%  
%     we are concerned with finding negative eigenvalues, which correspond
%     to possible instability modes.  Ignore the positive eigenvalues and
%     their eigenvectors
      if any(LAMBDA < -10*EPS)
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
%         the Jc rows are assumed zero
          n_Negs = zeros(n,Number_Neg_Eigenvalues);
          lambda_Jc = zeros(length(Jc),1);
%
%         check the LCP inequalities
          for i =1:n_Negs
            nu = [lambda_Jc; ...
                  Neg_Eigenvectors(J,i)];
            if    all(S(Jc,:)*nu <  EPS) ...
               && all(S(J ,:)*nu > -EPS)
              lambda = [lambda_Jc; ...
                        S(J ,:)*nu];
            elseif   all(S(Jc,:)*nu > -EPS) ...
                  && all(S(J ,:)*nu <  EPS)
              lambda = [ lambda_Jc; ...
                        -S(J ,:)*nu];
            end
%
            if ~isempty(lambda)
%             append the solution "lambda" to the other solutions
              Lambda = [Lambda lambda];
              NegEigenGrp = [NegEigenGrp, r];
              NegEigenValue = [NegEigenValue, LAMBDA(Neg_Eigenvalues(i))];
              NegEigenVector = [NegEigenVector, ...
                                -Q*lambda];
            else
              lambda = [];
            end
          end % i =1:Number_Neg_Eigenvalues
        end % if Number_Neg_Eigenvalues > 0
      end % if any(LAMBDA < -10*EPS)
    end % for r = 1:2^n
  end % if ~(n==0)
