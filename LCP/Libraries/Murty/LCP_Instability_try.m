function [ Lambda, NegEigenGrp, NegEigenValue, NegEigenVector] = ...
         LCP_Instability_try(P_L, H_x_x, H_x_lambda, ...
                             H_lambda_x, H_lambda_lambda, EPS);
%
%
  Lambda = [];
  NegEigenGrp = [];
  NegEigenValue = [];
  NegEigenVector = [];
%
% number of active contacts, M_lambda
  n = size(H_x_lambda, 2);
%
  N = dec2bin(0:(2^n-1),n);
  Ones = N=='1';
%
% check whether M_lambda is zero, and M_hat and q_hat are empty
  if n==0
    disp('M_hat and q_hat are empty.  Invalid LCP in LCP_enumerate2.m')
    W_2_matrix_sym = 0.5 * (H_x_x + H_x_x');
%
%   try Cholesky factorization to determine whether the matrix is positive
%   definite.  chol is much faster than the isdefinite function
    [R, Flag] = chol(W_2_matrix_sym);
%
    if Flag
%     if Flag is not zero, then W_2_matrix_sym is not positive definite, and
%     we can proceed to find the negative eigenvalues
%
%     find the negative eigenvalues and the corresponding eigenvectors
      [V, LAMBDA] = eig(W_2_matrix_sym);
%  
%     convert LAMBDA from a diagonal matrix to a vector
      LAMBDA = diag(LAMBDA);
%  
%     we are concerned with finding negative eigenvalues, which correspond
%     to possible instability modes.  Ignore the positive eigenvalues and
%     their eigenvectors
      if any(LAMBDA < -10*EPS)
%       isolate the eigenvectors that have negative eigenvalues.
        Neg_Eigenvalues = find(LAMBDA < -10*EPS);
%
        NegEigenValue = [NegEigenValue, LAMBDA(Neg_Eigenvalues)];
        NegEigenVector = [NegEigenVector, V(:,Neg_Eigenvalues)];
      end
    end
  else
    Q = inv(H_lambda_lambda) * H_lambda_x * P_L;
%
    I = [1:n];
    N = size(H_x_x, 1);
%
    S = P_L' * H_x_x' * P_L;
    T = P_L' * H_x_lambda;
%
%   search through all 2^n possible solutions
    for r = 1:2^n
%     these lambdas are possibly positive
      J = find(Ones(r,:));
      Jc = setdiff(I, J);
%
      if mod(r,5000)==0; disp(sprintf('Combination %6i of %6i', r, 2^n)); end
%
      B = zeros(n,N);
      B(J,:) = -Q(J,:);
%
%     2nd-order work for this complementary matrix. Note that inv(Ar) = Ar
      W_2_matrix = S +  T * B;
%
%     symmetric part of W_2_matrix
      W_2_matrix_sym = 0.5 * (W_2_matrix + W_2_matrix');
%     W_2_matrix_sym = balance(0.5 * (W_2_matrix + W_2_matrix'), 'noperm');
%
%     try Cholesky factorization to determine whether the matrix is positive
%     definite.  chol is much faster than the isdefinite function
      [R, Flag] = chol(W_2_matrix_sym);
%
      if Flag
%       if Flag is not zero, then W_2_matrix_sym is not positive definite, and
%       we can proceed to find the negative eigenvalues
%
%       find the negative eigenvalues and the corresponding eigenvectors
        [V, LAMBDA] = eig(W_2_matrix_sym);
%  
%       convert LAMBDA from a diagonal matrix to a vector
        LAMBDA = diag(LAMBDA);
%  
%       we are concerned with finding negative eigenvalues, which correspond
%       to possible instability modes.  Ignore the positive eigenvalues and
%       their eigenvectors
        if any(LAMBDA < -10*EPS)
%
%         isolate the eigenvectors that have negative eigenvalues.
          Neg_Eigenvalues = find(LAMBDA < -10*EPS);
          Number_Neg_Eigenvalues = length(Neg_Eigenvalues);
%
%         eigenvectors are potential solutions.  However, they must 
%         also satisfy the LCP inequalities
%
%         the "q" vector in the LCP
          q = Q * V(:,Neg_Eigenvalues);
%
%         check the LCP inequalities
          for i =1:Number_Neg_Eigenvalues
%           with each negative eigenvalue, check whether the q indicates
%           compliance with the complementarity condition
            if all(q(Jc,i) >= -1*EPS) && all(-q(J,i) >= -1*EPS)
%             the "y" vector is non-negative.  The zero and non-negative 
%             "lambda" values are contained in "y"
              lambda = zeros(n,1);
              lambda(J) = -q(J,i);
%
%             append the solution "lambda" to the other solutions
              Lambda = [Lambda lambda];
              NegEigenGrp = [NegEigenGrp, r];
              NegEigenValue = [NegEigenValue, LAMBDA(Neg_Eigenvalues(i))];
              NegEigenVector = [NegEigenVector, V(:,Neg_Eigenvalues(i))];
            elseif all(q(Jc,i) <= 1*EPS) && all(-q(J,i) <= 1*EPS)
%             the "y" vector is non-postive.  The zero and non-negative 
%             "lambda" values are contained in "-y"
              lambda = zeros(n,1);
              lambda(J) = q(J,i);
%
%             append the solution "lambda" to the other solutions
              Lambda = [Lambda lambda];
              NegEigenGrp = [NegEigenGrp, r];
              NegEigenValue = [NegEigenValue, LAMBDA(Neg_Eigenvalues(i))];
              NegEigenVector = [NegEigenVector, -V(:,Neg_Eigenvalues(i))];
            end
          end % i =1:Number_Neg_Eigenvalues
        end % if any(LAMBDA < -10*EPS)
      end % if Flag
    end % for r = 1:2^n
  end % if ~(n==0)
