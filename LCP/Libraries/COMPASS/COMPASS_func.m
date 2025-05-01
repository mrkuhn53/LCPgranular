function [x_k] = COMPASS_func(Mbar, dQbar, x0, ...
                            max_major, max_minor, max_restarts, target)
%
% An adaptation of the COMPASSmain.m file of Stefan Schmelzer, in the
% form of an Octave function.  The header information of Stefan Schmelzer
% is given below.  This adaptation is free software: you can redistribute 
% it and/or modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the License, 
% or (at your option) any later version.
%
% Matthew R. Kuhn
% 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPASSmain: Main file of the COMPASS routine
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

% COMPASS is an algorithm to solve the mixed complementarity problem
% (MCP):

% given  f: R^n -> R^n, and (possibly infinite) bounds l,u 
% find   z in [l,u], w,v in R_+^n
% s.t.:  f(z) = w-v
%        <(z-l),w> = 0
%        <(u-z),v> = 0

% A crash technique is used at the beginning of the algorithm in order
% to determine an approximation to the active set at the solution. 

% Solving the normal equation, a reformulation of the MCP in terms of
% the nonsmooth normal map, as proposed by Robinson, is the core of the
% algorithm. A Newton like method is used:

% A general first order approximation to the normal map is
% computed. The problem of finding a zero of this approximation is
% recast as a linear MCP. This linear MCP is solved by a pivot
% technique similar to that of Lemke, or that described by Dantzig.  
% The pivot procedure (inner, or minor algorithm) yields a piecewise
% linear path connecting the current iterate and the Newton point. Each
% pivot step results in a linear piece of the path. 

% If the new iterate does not provide the necessary descent in the
% merit function, the constructed path (or the line segment connecting
% its endpoints) is searched for an appropriate point.  
% A Non monotone stabilization scheme (NMS) and a watchdog techniqe are
% used in order to reduce the number of function and gradient
% evaluations necessary (outer, or major algorithm).  

% This procedure is repeated until the normal equation, and hence the
% MCP, is solved. 

% An example for how to use the program is in COMPASS_driver.m

% The program was written by Stefan Schmelzer (University of Vienna) as
% part of his diploma thesis. The current Version 1.0 is dated from
% April 16, 2012. 
 
% Please inform the author at schmelzer@ihs.ac.at if you make 
% serious use of this code. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kuhn: customize for an LCP (not the MLCP of Schmelzer's code)
  n=length(dQbar);            % problem dimension
  l = zeros(1,n);             % upper and lower bounds
  u = 1e80*ones(1,n);
  f=@(x) (Mbar*x' + dQbar)';  % the objective function
  Df=@(x) Mbar;               % gradient of the objective function
%
% Kuhn: I have selected minimal displayed output
display_depth = 0;     % display details on solving process:
                         % 0 -solution only
                         % 1 -major iteration info
                         % 2 -major and minor iteration info

%%%%%%%%%%%%%%%%
% Options:            % (can be changed by user; some are automatically
                      % changed after a restart)   

% crash technique
do_crash = 1;           % flag for running/stopping the crash procedure 
min_dim_crash = 10;     % minimum problem dimension: n=10
sigma_crash = 0.05;     % minimum descent in crash technique  
max_m_alpha_crash = 14; % minimum step size in crash technique is
                        % (0.5)^12   
max_crash_count = 50;   % maximum number of crash iterateions
min_active_change = 1;  % minimum change in active index set 
max_no_change_row = 3;  % max number of consecutive steps with
						% insufficient change in active set.  

% pivot technique
in_minoriter = 0;      % flag: is the algorithm in the file minoriter? 
in_pivot = 0;		   % flag: is the algorithm in the pivot procedure?  
infty=1000;            % number used as infinity border
cyclesize = min(max(5, n/10), 20);   % size that the length of a cycle
               						 % (of pivotroutine) is assumed.   
maybe_cycle = 2;       % iteration at which cycling possibility is
					   % checked.   
definitive_cycling = 0.5*n;   % maximum number of the same variable
							  % entering without abort due to cycling   
max_no_solution_LCP_row = 10; % maximum number of unsuccessful pivot
							  % procedures in a row before a restart.   

% NMS
pathsearch_type = 1;   % chose: line(1) or arc(2)
pg_search_type = 2;    % chose: line(1) or arc(2)
sigma_NMS = 0.01;      % minimum descent in pg and watchdog searches
max_m_alpha_w = 20;    % minimum step size in watchdog searches
max_m_alpha_pg = 25;   % minimum step size in pg searches
max_d_count = 3;       % maximum number of consecutive d-steps 
Delta0 = 0.5;          % measure of closeness in every first of
                       % consecutive d-steps  
deltachange = 0.07;    % Delta is reduced by deltachange after
                       % consecutive d-steps 

% restarts: 
next_x_near = 0;       % flag to take a special x0 after restart 
x_next = norm(u-l)*rand(1,n); % default special x0 after restarts 
max_restarts = min(60,max(12,max_restarts));
					   % set max_restarts to min.12 and max.60 restarts. 

% statistics:
restarts = 0;          % number of restarts in the algorithm
tot_major = [];        % major iterations per model run (before a 
                       % restart is made)  
save_avg_minor = [];   % average pivot steps per model run
save_Psi =[];          % last Psi_value before restart or total stop.

% Kuhn: initialize Var and search_type_change
Var = [];
search_type_change = [];
checkpointxN = [];
checkpointzN = [];
take_pg_step = [];

% initial settings: 					
x_k = x0;                            % current point
z_k = piB(x_k,l,u);                  % projection onto box B 
Psi_value_0 = Psi(z_k,l,u,f,infty);  % merit fct value at x0.
Psi_value_l = Psi(l,l,u,f,infty);    % merit fct value at l.
Psi_value_u = Psi(u,l,u,f,infty);    % merit fct value at u.
Psi_value_init = Psi_value_0;
if display_depth ~= 0
  Psi_value_0
  Psi_value_init
end


%%%%%%%%%%%%%%%%%%
% call the crash technique:
% crash % Kuhn: replace the "crash.m" script with the "crash_kuhn.m" function
% Computes an approximation to the active set at the solution. 
% Returns z_k, x_k, w_k and v_k
%
% Kuhn: replace the "crash.m" script with the "crash_kuhn.m" function
[ f_z_k, a_k, x_k, v_k, w_k, z_k, ...
  Psi_value, ...
  no_change_row, crash_count, do_crash, ...
  mark_check, bestpoint, Ref, memory, checkpointx, checkpointz] = ...
crash_func( ...
  n, x_k, z_k, l, u, f, Df, Psi_value_0, ...
  do_crash, min_dim_crash, display_depth, ...
  sigma_crash, max_m_alpha_crash, max_crash_count, ...
  min_active_change, max_no_change_row, infty);
%
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%
% call the core of the algorithm
% majoriter  % Kuhn: replace the "majoriter.m" script with majoriter_func.m"
% iteratively solves MCP
%
% Kuhn: replace the "majoriter.m" script with majoriter_func.m" function
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
%%%%%%%%%%%%%%%%%%%

if restart == 1    
 % The algorithm could not solve the problem with the initial settings.
 % It is restarted with different options, multiple times if necessray.

  while restart == 1               

    %%%%%%%%%%%%%%%%
    % call the restart routine:
%   RestartCOMPASS
    % different options are defined
    % crash technique is called
    % majoriter is called. 
% Kuhn: replace the "RestartCOMPASS.m" script with 
% RestartCOMPASS_func.m" function
    [ a_k, x_k, v_k, w_k, z_k, ...
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
      display_depth, max_restarts, max_major, max_minor, target);
    %%%%%%%%%%%%%%%%%%%%%

  end
end
%
% Kuhn: return the solution as a column  vector
x = z_k';     % solution
