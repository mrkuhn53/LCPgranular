function [ M_hat, q_hat, H_x_BD_inv, RCond_BD] = ...
         Setup_LCP_Three_Disks( ...
           H_x, H_lambda, dp, P_L, P_Lperp, Q_x, Q_lambda, C_MP_inv, dc)
%
%   with the LCP, find the matrix M_hat and vector q_hat
%
%   pre-Bott-Duffin inverse
    [Pre_H_x_BD_inv, RCond_BD] = inv(H_x*P_L + P_Lperp);
%
%   check whether the matrix (H_x*P_L + P_Lperp) is singular.
%   If so, we must use the generalize Bott-Duffin inverse.  Because the
%   matrix is singular, we must use the SVD to find the Bott-Duffin inverse
    if abs(RCond_BD) < 10*eps
      [U, S, V] = svd(H_x*P_L + P_Lperp);
%
      S_diag = diag(S);
      S_diag_zeros = find(abs(S_diag) < 10*eps);
      S_diag_nonzeros = setdiff([1:length(S_diag)], S_diag_zeros);
%
      S_diag_inv = zeros(size(S_diag));
      S_diag_inv(S_diag_zeros) = 0;
      S_diag_inv(S_diag_nonzeros) = 1 ./ S_diag(S_diag_nonzeros);
%
      S_inv = diag(S_diag_inv);
%
      Pre_H_x_BD_inv = V * S_inv * U';
    end
%
%   Bott-Duffin inverse
    H_x_BD_inv = P_L * Pre_H_x_BD_inv;
%
%   the M and q matrices of the LCP problem
    M_hat = Q_lambda - Q_x * H_x_BD_inv * H_lambda;
%
    q_hat =   Q_x * (C_MP_inv - H_x_BD_inv * H_x * C_MP_inv) * dc ...
            + Q_x * H_x_BD_inv * dp;
