function [ a_k, x_k, v_k, w_k, z_k, ...
           Psi_value, Psi_value_0, x_next, next_x_near, ...
           in_minoriter, avg_minor, restarts, ...
           search_type_change, pathsearch_type, pg_search_type, ...
           do_crash, sigma_crash, max_m_alpha_w, ...
           take_pg_step, Delta0, memory, Ref, in_pivot, ...
           mark_check, bestpoint, checkpointx, checkpointz, ...
           checkpointxN, checkpointzN, ...
           sigma_NMS, deltachange, ...
           stop, restart] =  ...
         RestartCOMPASS_func( ...
           n, l, u, f, Df, infty, ...
           a_k, x_k, v_k, w_k, z_k, ...
           Psi_value, Psi_value_0, Psi_value_u, Psi_value_l, ...
           Delta0, x_next, next_x_near, memory, Ref, in_pivot, ...
           in_minoriter, avg_minor, restarts, tot_major, Var, bestpoint, ...
           save_avg_minor, save_Psi, cyclesize, maybe_cycle, mark_check, ...
           search_type_change, pathsearch_type, pg_search_type, ...
           do_crash, sigma_crash, max_m_alpha_w, max_m_alpha_pg, ...
           sigma_NMS, deltachange, max_d_count, max_no_solution_LCP_row, ...
           take_pg_step, checkpointxN, checkpointzN, ...
           checkpointx, checkpointz, stop, ...
           max_m_alpha_crash, min_dim_crash, definitive_cycling, ...
           max_no_change_row, min_active_change, max_crash_count, ...
           display_depth, max_restarts, max_major, max_minor, target)
%
% An adaptation of the RestartCOMPASS.m file of Stefan Schmelzer, in the
% form of an Octave function.  The header information of Stefan Schmelzer
% is given below.  This adaptation is free software: you can redistribute
% it and/or modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the License,
% or (at your option) any later version.
%
% Matthew R. Kuhn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restartCOMPASS: COMPASS is restarted with different options.  
% Copyright (C) 2012 Stefan Schmelzer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This file is part of COMPASS.

% COMPASS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% COMPASS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
  
% You should have received a copy of the GNU General Public License
% along with COMPASS.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%
% Comments on this file: 

% goes back to x0 
% changes crucial parameter settings, 
% calls crash.m and majoriter.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


if display_depth ~= 0
  disp('')
  disp('          RESTART')
  disp('')
end

restarts = restarts +1; %count the restarts. if they exceed max_restarts, COMPASS
             % terminates.   
restart = 0; % exit the loop in COMPASSmain if majoriter.m terminates
             % successfully. (if another restart is necessary, majoriter
             % sets restart=1 again.   

if restarts <= max_restarts

  % reset statistic counters:
  artificial_count = 0;
  avg_minor = 0;

  % restart parameters  
  rsp(1,:) = [1 2 1 0.05  16 0.05 6 0.8];            
  rsp(2,:) = [2 2 1 0.05  18 0.01 8 0.7];
  rsp(3,:) = [1 1 1 0.05  22 0.01 9 0.9];
  rsp(4,:) = [2 1 1 0.05  20 0.01 5 0.3];
  rsp(5,:) = [1 2 1 0.05  15 0.05 6 0.8];              
  rsp(6,:) = [2 2 1 0.05  16 0.01 8 0.7];
  rsp(7,:) = [1 1 1 0.05  19 0.01 9 0.9];
  rsp(8,:) = [2 1 1 0.05  17 0.01 5 0.3];
  rsp(9,:) = [1 2 1 0.05  18 0.01 6 0.8];              
  rsp(10,:)= [2 2 1 0.05  20 0.01 8 0.7];
  rsp(11,:)= [1 1 1 0.05  22 0.01 9 0.9];
  rsp(12,:)= [2 1 1 0.05  19 0.01 5 0.3];
  
  if restarts >= 12
    kk = 12;
  else
    kk = restarts;
  end

  % change parameters:
  pathsearch_type = rsp(kk,1); % chose: line(1) or arc(2) 
  pg_search_type = rsp(kk,2);  % chose: line(1) or arc(2)
  do_crash = rsp(kk,3);        % flag for invoking/stopping the crash procedure
  sigma_crash = rsp(kk,4);     % minimum descent in crash technique  
  max_m_alpha_w = rsp(kk,5);   % exponent of minimum step size (pg&wd searches) 
  sigma_NMS = rsp(kk,6);       % minimum descent in pg and watchdog searches
  %max_d_count = rsp(kk,7);    % maximum number of consecutive d-steps 
  deltachange = rsp(kk,8);     % Delta is reduced by deltachange after
                               % consecutive d-steps 

  % change initial point: 
  if next_x_near == 1
    x0 = x_next;
  else
    x0 = norm(u-l)*rand(1,n);
    if rsp(kk,1) == 1  % try astarting point inside B every second time. 
      x0 = piB(x0,l,u);
    end
  end
  next_x_near = 0;

  x_k = x0;                      % new current point
  z_k = piB(x_k,l,u);
  Psi_value_0 = Psi(z_k,l,u,f,infty); % new initial merit fct value.
  if display_depth ~= 0
    Psi_value_0
  end


  %%%%%%%%%%%%%%%%
% crash          
  % Computes a good approximation to the active set quickly with projected Newton
  % steps. Returns z_k, x_k, w_k and v_k
% Kuhn: change crash.m script to crash_func.m function
  [ f_z_k, a_k, x_k, v_k, w_k, z_k, ...
    Psi_value, ...
    no_change_row, crash_count, do_crash, ...
    mark_check, bestpoint, Ref, memory, checkpointx, checkpointz] = ...
  crash_func( ...
    n, x_k, z_k, l, u, f, Df, Psi_value_0, ...
    do_crash, min_dim_crash, display_depth, ...
    sigma_crash, max_m_alpha_crash, max_crash_count, ...
    min_active_change, max_no_change_row, infty);
  %%%%%%%%%%%%%%%%%%%%%%% 


  %%%%%%%%%%%%%
% majoriter
  % core of COMPASS algorithm - iteratively solves MCP
% Kuhn: change majoriter.m script to majoriter_func.m function
  [ x_k, v_k, w_k, z_k, a_k, f_z_k, ...
    stop, restart, sigma_NMS, take_pg_step, avg_minor, ...
    in_minoriter, pathsearch_type, ...
    Delta, Delta0, d_count, next_x_near, x_next, ...
    memory, Ref, in_pivot, search_type_change, ...
    mark_check, bestpoint, checkpointx, checkpointz, ...
    checkpointxN, checkpointzN] = ...
  majoriter_func( ...
    n, z_k, f, Df, l, u, ...
    f_z_k, a_k, x_k, v_k, w_k, ...
    max_major, target, in_minoriter, Var, ...
    Psi_value, Psi_value_0, Psi_value_u, Psi_value_l, ...
    Delta0, next_x_near, x_next, memory, Ref, in_pivot, ...
    restarts, tot_major, save_avg_minor, save_Psi, ...
    mark_check, bestpoint, checkpointx, checkpointz, ...
    checkpointxN, checkpointzN, ...
    max_d_count, max_no_solution_LCP_row, definitive_cycling, ...
    display_depth, infty, cyclesize, maybe_cycle, ...
    pathsearch_type, sigma_NMS, pg_search_type, search_type_change, ...
    max_minor, max_m_alpha_pg, max_m_alpha_w, deltachange);
  %%%%%%%%%%%%%%%%%%


else
  restarts = max_restarts;
  disp('          max number of restarts reached')
end

% PARENTFILE: COMPASSmain
