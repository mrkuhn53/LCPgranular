function [ lambda, Lambda, fval, ...
           DegenSoln, IsolatedSoln, NonIsolated, ...
           ColCompetent, DegenMatrix, n_DegenMatrix, ...
           is_PD, is_NS, Method_used, is_Soln] = ...
         Solve_LCP_alt( ...
           M_hat, q_hat, x0, Method, PHI, Hc, System, LibrariesPath, ...
           M_lambda_max_enumerate, Only_Fast_Method, ...
           lambda_try, n_intersect_lambda, EPS);
%
  is_PD = -1;
  is_NS = -1;
%
  M_lambda = length(q_hat);
%
  EPS_is_Soln = 1e-10;
%
% List of Octave's Methods that can be used to solve the LCP, and the folder
% locations that contain the library
% qpOctave = 1; % Octave's qp() solver
% Quadprog = 2; % Octave's quadprog() solver
% COMPASS =  3; % COMPASS solver
% Murty =    4; % LCP_enumerate() solver
% Yuval =    5; % Mathworks code by Yuval
% Almqvist = 6; % Mathworks code by Andreas Almqvist
% CompEcon = 7; % Paul L. Fackler & Mario J. Miranda
% sqpOctave =8; % Octave's sqp() solver
%
  Method_Fast =     1;
  Method_Fast_Alt = 8;
  Method_Slow =     4;
%
  is_Soln = [];
%
  if M_lambda ~= 0
%   whether M_hat is non-singular
    is_NS = isNonSingular(M_hat,EPS);
%
%   whether M_hat is positive definite
    is_PD = isdefinite(0.5*(M_hat+M_hat'),EPS);
%
%   older versions of octave gave a value of "isdefinite" = -1 if its 
%   argument was not positive semi-definite.  We want a value of 0 if
%   the argument is not positive seme-definite
    is_PD = is_PD && ~(is_PD==-1);
%
    if isempty(Method)
      if is_PD || M_lambda > M_lambda_max_enumerate || Only_Fast_Method
        Method = Method_Fast; % octave's built-in qp() solver
      else
        Method = Method_Slow; % a full enumerative calculation of all solutions
      end
    end
  end
%
% solve the problem as an LCP
  [ lambda, Lambda, fval, ...
    DegenSoln, IsolatedSoln, NonIsolated, ...
    ColCompetent, DegenMatrix, n_DegenMatrix] = ...
  LCP3_Solve_Problem_alt( ...
    M_hat, q_hat, x0, Method, PHI, Hc, System, LibrariesPath);
%
% is 'lambda' a solution
  is_Soln = isSoln(M_hat, q_hat, lambda, EPS_is_Soln);
%
  if Method==Method_Slow && isempty(lambda)
%
    Method = Method_Fast_Alt;
%
    if n_intersect_lambda > 3
      x0 = lambda_try;
    end
%
%   solve the problem as an LCP
    [ lambda, Lambda, fval, ...
      DegenSoln, IsolatedSoln, NonIsolated, ...
      ColCompetent, DegenMatrix, n_DegenMatrix] = ...
    LCP3_Solve_Problem_alt( ...
      M_hat, q_hat, x0, Method, PHI, Hc, System, LibrariesPath);
%
%   is 'lambda' a legitimate solution
    is_Soln = isSoln(M_hat, q_hat, lambda, EPS_is_Soln);
%
  elseif M_lambda~=0 ...
         && Method==Method_Fast ...
         && ~is_Soln
%
%   we tried using the fast solution method, but it did not converge onto
%   a valid solution
%   Method = Method_Fast_Alt;
%
    if is_PD && M_lambda <= M_lambda_max_enumerate
      Method = Method_Slow; % a full enumerative calculation of all solutions
    else
      Method = Method_Fast_Alt;
    end
%
%   solve the problem as an LCP using the slow algorithm
    [ lambda, Lambda, fval, ...
      DegenSoln, IsolatedSoln, NonIsolated, ...
      ColCompetent, DegenMatrix, n_DegenMatrix] = ...
    LCP3_Solve_Problem_alt( ...
      M_hat, q_hat, x0, Method, PHI, Hc, System, LibrariesPath);
%
%   is 'lambda' a legitimate solution?
    is_Soln = isSoln(M_hat, q_hat, lambda, EPS_is_Soln);
  end
%
if 0
  if M_lambda ~= 0 ...
     && (Method==Method_Fast || Method==Method_Fast_Alt) ...
     && ~is_Soln ...
     && M_lambda<=M_lambda_max_enumerate ...
     && ~Only_Fast_Method
%
%   we tried using the fast solution method, but it did not converge onto
%   a valid solution
    Method = Method_Slow;
%
%   solve the problem as an LCP using the slow algorithm
    [ lambda, Lambda, fval, ...
      DegenSoln, IsolatedSoln, NonIsolated, ...
      ColCompetent, DegenMatrix, n_DegenMatrix] = ...
    LCP3_Solve_Problem_alt( ...
      M_hat, q_hat, x0, Method, PHI, Hc, System, LibrariesPath);
%
%   is 'lambda' a legitimate solution?
    is_Soln = isSoln(M_hat, q_hat, lambda, EPS_is_Soln);
%
    if isempty(lambda)
      lambda = zeros(size(q_hat));
    end
  end
end
%
%
% is 'lambda' a legitimate solution?
  is_Soln = isSoln(M_hat, q_hat, lambda, EPS_is_Soln);
  if is_Soln && Method==Method_Slow
    is_Soln = size(Lambda,2);
  end
%
  Method_used = Method;
