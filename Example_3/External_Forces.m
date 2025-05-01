function dp = External_Forces(N, b_p, w_p)
%
% with the purpose of equilibrating the contact forces, apply the net
% force and moment imbalances as external forces
%
%-------------- INPUT -------------------------
% N = number of particles, excluding non-relevant particles, and including
%     ghost particles
% b_p = Nx3 force imbalance on particles
% w_p = Nx3 moment imbalance on particles
%
%------------- OUTPUT -------------------------
% dp = 6Nx1 vector of forces, db_p and dw_p
%
% in the event of small net force and moment imbalances, apply these
% imbalances as external forces, so that the system will equilibrate
% during this time step
%
% the force increment vector
  dp = zeros(6*N,1);
  dp(1:3:3*N,1) = b_p(:,1);
  dp(2:3:3*N,1) = b_p(:,2);
  dp(3:3:3*N,1) = b_p(:,3);
%
  dp(3*N + [1:3:3*N],1) = w_p(:,1);
  dp(3*N + [2:3:3*N],1) = w_p(:,2);
  dp(3*N + [3:3:3*N],1) = w_p(:,3);
%
% dp = 0.2*dp;
% dp = 0.9 * dp;
