%
function ...
  [ H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, dp_dstress] = ...
  LCP3_Convert_2D_or_3D_alt( ...
    N_Rel, H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, dp_dstress, idim);
%
% this function converts matrices for a 3D problem into a 2D problem.
% These matrices will be arranged for the following ordering of the
% displacements and rotations (see Eqs. 1 & 2, with movements stacked
% above rotations.  In each part, particle movements (or rotations) are
% arranged in primary order of particle number, and secondary order of 
% movement component.
%
% [du_1_1; du_1_2; du_2_1; du_2_2; du_3_1; du_3_2; etc. dtheta_1_3; dtheta_2_3;
% etc.]
%
% where the first subscript is the particle number and the second subscript
% is the movement component.
%
%---------------------  INPUT & OUTPUT ---------------------------
% idim =    dimension of the problem.  Either 2 or 3.
% Rows =    rows (and columns) to be extracted for problem solution
%
%
  if idim==2
%   number of particles.  The matrix
    N = N_Rel;
%
%   select the rows of the stacked vector [du / dtheta] that apply to 
%   the 2D case
    Rows2D = [ 1:3:3*N, 2:3:3*N, (3*N+3):3:6*N, 6*N+[1,2,4]];
    Rows2D = sort(Rows2D);
%
%   select the 2D rows and columns of these matrices
    H_x_x =           H_x_x(Rows2D, Rows2D);
    H_x_lambda =      H_x_lambda(Rows2D, :);
    H_lambda_x =      H_lambda_x(:, Rows2D);
    H_lambda_lambda = H_lambda_lambda;
%
    dp_dstress =      dp_dstress(Rows2D);
  end
%
