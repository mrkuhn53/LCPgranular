function [ a_k, v_k, w_k, z_k, z_b, ...
           stop, restart, next_x_near, x_next, ...
           search_type_change, memory, in_pivot, ...
           pathsearch_type] = ...
         watchdog_step_func( ...
           n, l, u, f, Df, infty, ...
           bestpoint, zwvt, x_next, ...
           a_k, v_k, w_k, z_k, z_b, memory, in_pivot, ...
           checkpointx, checkpointz, ...
           checkpointxN, checkpointzN, search_type_change, ...
           Ref, sigma_NMS, pg_search_type, ...
           pathsearch_type, display_depth, max_m_alpha_w, ...
           max_m_alpha_pg, stop, restart, ...
           next_x_near)
%
% An adaptation of the watchdog_step.m file of Stefan Schmelzer, in the
% form of an Octave function.  The header information of Stefan Schmelzer
% is given below.  This adaptation is free software: you can redistribute
% it and/or modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the License,
% or (at your option) any later version.
%
% Matthew R. Kuhn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% watchdog_step: A step kind in the NMS procedure in COMPASS.  
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

% A watchdog-step is taken if descent criterion is violated in m-steps,
% or if Psi(x_k_N) is too high in d-steps.   

% The line segment between checkpointz and checkpointzN (or
% alternatively the arc that results from projecting the line between
% checkpointx and checkpointxN onto B) is searched for a point that
% satisfies the nonmonotone descent criterion. If the stepsize becomes
% too small, a projected gradient step is taken.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


grad_Psi_chz = DPsi(checkpointz,l,u,f,Df,infty)';
								% value of gradient of Psi at checkpointz

m = 0;
search = 1;
search_type_change = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while search == 1  
  alpha_w = (0.5)^m;

  if pathsearch_type == 1       % z_alpha can be defined in two ways
    z_alpha = (checkpointz + alpha_w*(checkpointzN-checkpointz));
  elseif pathsearch_type == 2
    z_alpha = piB(checkpointx + alpha_w*(checkpointxN-checkpointx),l,u);
  end

  if ( grad_Psi_chz*(checkpointz - z_alpha)' < 0 )   
    % criterion which criterion to apply:
    % descent criterion 1
    if (Psi(z_alpha,l,u,f,infty)  >  Ref - sigma_NMS* grad_Psi_chz*(checkpointz  - z_alpha)' )
      m = m+1;
    else 
      search = 0;
      z_k = z_alpha;
      if display_depth ~= 0
        disp(' watchdog-step taken (case_1)')
      end
    end 
  else
    % descent criterion 2
    if Psi(z_alpha,l,u,f,infty)  >  (1-alpha_w*sigma_NMS)*Ref
      m = m + 1;
    else 
      search = 0;
      z_k = z_alpha;
      if display_depth ~= 0
        disp(' watchdog-step taken (case_2)')
      end 
    end
  end	      
  
  if m == max_m_alpha_w
    if search_type_change == 0  
      m = 0;
      if pathsearch_type == 1  % change searchtype, try again.    
        pathsearch_type = 2;
      elseif pathsearch_type == 2
        pathsearch_type = 1;
      end
      search_type_change = 1;
    else
      search = 0;
      if display_depth ~= 0
        disp('alpha_w   too small in watchdog-step')
      end

      %%%%%%%%%%%%%%%%%%
%     pg_step
      % Returns z_k, new checkpoint
% Kuhn: change pg_step.m script to pg_step_func.m function
      [ a_k, v_k, w_k, z_k, z_b, ...
        stop, restart, next_x_near, x_next, ...
        search_type_change, memory, in_pivot] = ...
      pg_step_func( ...
        n, bestpoint, l, u, f, Df, zwvt, infty, ...
        a_k, v_k, w_k, z_k, z_b, ...
        sigma_NMS, pg_search_type, max_m_alpha_pg, ...
        display_depth, stop, restart, next_x_near, ...
        search_type_change, memory, in_pivot);
      %%%%%%%%%%%%%%%%%%%%%
    end
  end
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if search_type_change == 1
  search_type_change = 0;
  if pathsearch_type == 1    % change back
    pathsearch_type = 2;
  elseif pathsearch_type == 2
    pathsearch_type = 1;
  end
end

% PARENTFILE: NMS.m
