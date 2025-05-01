function [ H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, ...
           dxcell_vec, dp_dstress] = ...
         Combine_H_alt( ...
           H_x, H_lambda, H_stress_x, H_stress_lambda, Q_x, Q_lambda, ...
           dp, dstress_vec, dxcell);
%
% arrange dxcell_vec as a column vector
  dxcell_vec = [dxcell(1,1);
                dxcell(2,2);
                dxcell(3,3);
                dxcell(1,2);
                dxcell(1,3);
                dxcell(2,3)];
%
  H_x_x =           [H_x; ...
                     H_stress_x];
%
  H_x_lambda =      [H_lambda; ...
                     H_stress_lambda];
%
  H_lambda_x =      [Q_x];
  H_lambda_lambda = [Q_lambda];
%
  dp_dstress =      [dp; ...
                     dstress_vec];
