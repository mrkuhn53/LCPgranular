function [Combinations, H_g4, c, f, dx_c, C, dc, ...
          R, Prr, Pnrr, other_parameters] = ...
         Misc_Initialization(N)
%
% list of the combinations of contact stiffness branches
% that should be investigated.  A zero "0" or empty vector
% [] means to investigate all possible branches.  A negative
% value means to investigate no branches.
  Combinations = [];
%
% external follower forces
  H_g4 = [];
%
% initialize some vectors associated with constraints
  c = []; f = []; dx_c = []; C = []; dc = []; R = []; Prr = []; Pnrr = [];
%
% initialize the increments of force and moment on the particles.  The
% vector [dp] is the stacked vector of force increments [db] and moment
% increments [dw], as in KPD Eq. 2_2.  Alternatively, with external follower
% forces, [dp] can represent the lower-dimensional increments in KPD Eq. 5_2
% dp = zeros(6*N,1);
%
% future contact parameters
  other_parameters = [];
