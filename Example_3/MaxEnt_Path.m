function [ lambda] = ...
         MaxEnt_Path( ...
           Lambda, Choose_MaxEnt_Path, is_Soln, ...
           H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, ...
           H_xx_BD_inv, C_MP_inv, dxcell_vec, dp_dstress, dc_w_dxcell);
%
% if ~is_Soln
%   the LCP does not have a solution that will maintain equilibrium in the
%   deformed configuration.  Look for possible instabilities of equilibrium
%   in the undeformed configuration.  If such instability exists, we will
%   temporarily halt the loading (dxcell, dp, etc) and allow the system to
%   proceed to a non-equilibirum configuration, but from which equilibrium
%   can be reestablished upon resumption of loading
%   [ Lambda_Instab, NegEigenGrp, NegEigenValue, NegEigenVector] = ...
%     LCP_Instability2(P_L, H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, ...
%                      1000*eps);
% end
%
  if Choose_MaxEnt_Path && size(Lambda,2) > 1
%   2nd order entropy production
    I_2 = zeros(size(Lambda,2),1);
%
    for i = 1:size(Lambda,2)
      lambda_trial = Lambda(:,i);
%
%     since we have solved lambda, find the particle movements
      dx = (C_MP_inv - H_xx_BD_inv * H_x_x * C_MP_inv) * dc_w_dxcell ...
            + H_xx_BD_inv * dp_dstress ...
            - H_xx_BD_inv * H_x_lambda * lambda_trial;
%
%     reaction forces
      dy = - dp_dstress ...
           + H_x_x * dx ...
           + H_x_lambda * lambda_trial;
%
%     2nd order entropy production
      I_2(i) = 0.5*(dy'*dx - dp_dstress'*dx);
    end
%
    if Choose_MaxEnt_Path
      [I_2_min, iLambda] = min(I_2);
    else
      iLambda = 1;
    end
    lambda = Lambda(:,iLambda);
%
  elseif size(Lambda,2)==1
    lambda = Lambda(:,1);
  else
    lambda = [];
  end
