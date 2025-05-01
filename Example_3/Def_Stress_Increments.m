function [ xcell, xcello, dxcell, def, defo, ddefm_apply, ...
           Vol, stress, stress_target, Reset_target_stress, Time, dTime, ...
           stress_targeto, dstress_apply, Repetitions] = ...
         Def_Stress_Increments( ...
           xcell, def, ddefm, r_pq_p, f_pq, stress, ...
           nStress, iStress, Reset_target_stress, stress_target, dstress, ...
           Time, Time_prev, b_p_Rel, is_Soln, I, J, idim, ...
           Max_stress_deviation, Max_b_p, Repetitions, Max_Repetitions)
%
% this function computes the increments of stress and deformation
% for the next time step
%
% the previous cell dimensions
  xcello = xcell;
%
% the previous deformation
  defo = def;
%
% time increment
  dTime = 1;
%
% current volume of assembly
  Vol = prod(diag(xcello)(1:idim));
%
% previous stress
  stresso = stress;
%
% current stress. 
  stress = (1/Vol) * r_pq_p' * f_pq;
%
  if nStress > 0 
    if Reset_target_stress
%     a new control period has been entered with som controlled stresses.  
%     Reset the target stress
      Reset_target_stress = 0;
      stress_target = stress;
    end
%
%   express the actual and target stresses as vectors, so that we can determine
%   whether the stress is deviating from the target stress
    stress_vec = zeros(1,6);
    stress_target_vec = zeros(1,6);
    for i = 1:6
      stress_vec(i) = stress(I(i),J(i));
      stress_target_vec(i) = stress_target(I(i),J(i));
    end
  end
%
% if stresses are being controlled...
  if nStress > 0 ...
    && ((max(abs(stress_vec(iStress) - stress_target_vec(iStress))) ...
          <= Max_stress_deviation ...
         && is_Soln ...
         && mean(mean(abs(b_p_Rel(:,[1:idim])))) <= Max_b_p ...
        ) ...
        || Repetitions >= Max_Repetitions ...
       )
%   the maximum deviation from the target stress is less than the tolerable
%   amount and a solution existed for the previous step.  Conntinue by using 
%   the prescribed changes in target stress and deformation
    dstress_apply = dstress;
    ddefm_apply = ddefm;
%
%   reset repetitions
    Repetitions = 0;
  elseif nStress > 0
%   stresses are being controlled, but the stress has deviated from the
%   target stress
    dstress_apply = zeros(3);
    ddefm_apply = zeros(3);
    Time = Time - 1;
    dTime = 0;
    Repetitions = Repetitions + 1;
  elseif nStress <= 0
%   stresses are not being controlled
    dstress_apply = zeros(3);
    ddefm_apply = ddefm;
    Repetitions = 0;
  else
%   stresses are not being controlled
    dstress_apply = zeros(3);
    ddefm_apply = zeros(3);
    Repetitions = 0;
  end
%
% advance the deformation gradient
  def = def + ddefm_apply;
%
% Eulerian increment of deformation
  ddefme = ddefm_apply*inv(def);
%
% increment in periodic cell dimensions
  dxcell = ddefme*xcell;
%
% advance the periodic cell dimensions
  xcell = xcell + dxcell;
%
  stress_targeto = stress_target;
  if nStress > 0
%   advance the target stress
    stress_target = stress_target + dstress_apply;
  end
