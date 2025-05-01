function [ Time, lambda, lambda_original, dp, dstress_vec, dc, dgamma, ...
           xcell, dxcell, dxcell_vec, def, stress_target] = ...
         Check_Lambda( ...
           lambda, lambda_threshold, ...
           Time, dTime, Time_prev, dp, dstress_vec, dc, dgamma, ...
           xcell, dxcell, xcello, dxcell_vec, def, defo, ddefm_apply, ...
           stress_target, stress_targeto, dstress_apply, nStress)
%
% this function checks whether max(lambda) exceeds a threshold value,
% and if so, the lambda results are scaled back to the threshold value
%
if max(abs(lambda)) > lambda_threshold
  if max(abs(lambda)) > 0
    rTime = lambda_threshold / max(abs(lambda));
  else
    rTime = 0
  end
%
  dTime = rTime * dTime;
%
% revise the current time
  Time = Time_prev + dTime;
%
% revise the lambda values
  lambda_original = lambda;
  lambda = rTime * lambda;
%
% revise the increments of force and constrained displacements
  dp = rTime * dp;
  dstress_vec = rTime * dstress_vec;
  dc = rTime * dc;
  dgamma = rTime * dgamma;
%
% advance the deformation gradient
  def = defo + rTime * ddefm_apply;
%
% Eulerian increment of deformation
  ddefme = ddefm_apply*inv(defo);
%
% increment in periodic cell dimensions
  dxcell = ddefme*xcello;
%
% advance the periodic cell dimensions
  xcell = xcello + dxcell;
%
% arrange dxcell_vec as a column vector
  dxcell_vec = [dxcell(1,1);
                dxcell(2,2);
                dxcell(3,3);
                dxcell(1,2);
                dxcell(1,3);
                dxcell(2,3)];
%
  if nStress > 0
%   advance the target stress
    stress_target = stress_targeto + rTime*dstress_apply;
  end
else
  lambda_original = lambda;
end
