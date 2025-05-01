%
% |+|+|+|+|+|+|+|+|  Revision 2   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function ...
  [H_x, H_m_x, H_g1, H_g2, H_g3, H_g4, H_g1234, H_m_lambda, Q_x, ...
   dp, c, f, C, R, Prr, Pnrr, P_L, P_Lperp, C_MP_inv, Rows, HRows, ...
   Affine_Matrix] = ...
  LCP_Convert_2D_or_3D( ...
    dp, H_x, H_m_x, H_g1, H_g2, H_g3, H_g4, H_g1234, H_m_lambda, Q_x, ...
    c, f, C, R, Prr, Pnrr, idim, P_L, P_Lperp, C_MP_inv, ...
    Affine_Matrix, Type_Constraint)
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
% number of particles.  The matrix
  N = size(dp, 1) / 6;
%
% the number of rows and columns in the global stiffness matrix
  HRows = 6*N;
%
% this vector is a list of all degrees of freedom that must be considered
% for the problem.  For 3D problem, it is 6Nx1.  For 2D problems it is 3Nx1
% Rows = ones(6*N,1);
  Rows = [1:6*N]';
%
  if idim==2
%
    if mod(N,1)~=0
      disp('Oops! The length of input vector [dp] should be integer.');
    end
%
%   select the rows of the stacked vector [du / dtheta] that apply to 
%   the 2D case
    Rows2D = [ 1:3:3*N, 2:3:3*N, (3*N+3):3:6*N];
    Rows2D = sort(Rows2D);
%
%   leave the function with this value of Rows
    Rows = Rows2D;
%
%   the number of rows and columns in the global stiffness matrix
    HRows = 3*N;
%
%   select the 2D rows and columns of these matrices
    dp = dp(Rows2D);
    H_g1 =    H_g1(Rows2D,Rows2D);
    H_g2 =    H_g2(Rows2D,Rows2D);
    H_g3 =    H_g3(Rows2D,Rows2D);
    H_x =     H_x(Rows2D,Rows2D);
    H_m_x =   H_m_x(Rows2D,Rows2D);
    H_m_lambda = H_m_lambda(Rows2D,:);
    Q_x        = Q_x(:,Rows2D);
    if ~isempty(H_g4)
      H_g4 =  H_g4(Rows2D,Rows2D);
    end
    H_g1234 = H_g1234(Rows2D,Rows2D);
%
%
%   initialize the free degrees
    f = [];
%
    if Type_Constraint==1
%     the vector "c" listed the elements of [du /dtheta] that were constrained.
%     This list was for a 3D system.  Convert "c" to a list for the 2D system.
%     See Eqs. 54-56.  Note that the order of these "c" constraints must 
%     correspond to the order of the dx_c constrained displacements.
%
      disp('Oops! Constraint Type II not currently supported.')
      ERROR
%
      for i = 1:length(c)
        if c(i) <= 3*N          % a displacement
          if mod(c(i),3)==1     %   an x_1 displacement
            c(i) = 2*floor(c(i)/3) + 1;
          elseif mod(c(i),3)==2 %   an x_2 displacement
            c(i) = 2*floor(c(i)/3) + 2;
          else                  %   an x_3 displacement
            disp('Error! An x_3 constraint in a 2D system.');
          end
        else                    % a rotation
          if mod(c(i),3)==0     %   a theta_3 rotation
            c(i) = 2*N + (c(i) - 3*N) / 3;
          else
            disp('Error! An theta_1 or theta_2 constraint in a 2D system.');
          end
        end
      end
%
%     the list of unconstrained movements
      f = setdiff([1:3*N], c);
%
    elseif Type_Constraint==2 
      disp('Oops! Constraint Type II not currently supported.')
      ERROR
%
      Type_Constraint==3
%     constraint matrix in Eqs. 57 and 64
      C = C(:, Rows2D);
      P_L = P_L(Rows2D,Rows2D);
      P_Lperp = P_Lperp(Rows2D,Rows2D);
      C_MP_inv = C_MP_inv(Rows2D,:);
%
    elseif Type_Constraint==3
%     constraint matrix in Eqs. 57 and 64
      C = C(:, Rows2D);
      P_L = P_L(Rows2D,Rows2D);
      P_Lperp = P_Lperp(Rows2D,Rows2D);
      C_MP_inv = C_MP_inv(Rows2D,:);
%
    elseif Type_Constraint==4
      disp('Oops! Constraint Type II not currently supported.')
      ERROR
%
%     matrices in Eqs. 72-74
      R = R(Rows2D,3);
      Prr = Prr(Rows2D,Rows2D);
      Pnrr = Pnrr(Rows2D,Rows2D);
    end
%
    Affine_Matrix = Affine_Matrix(Rows2D(1:2*N),[1,2,4,5,7,8]);
  end
