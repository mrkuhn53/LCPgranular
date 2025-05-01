%
% |+|+|+|+|+|+|+|+|+|  Revision 3   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
  clear
% Example with a 3D Body Centered Cubic system of spheres
%
% In comments below, reference is made to equations in KPD, which is the
% Kuhn, Prunier, Daouadji paper: "Stiffness pathologies in discrete 
% granular systems: Bifurcation, neutral equilibrium, and instability 
% in the presence of kinematic constraints," Int. J. Num. Anal. Methods Geomech.,
% 2019
%
% connect with libraries
  addpath('../mfiles2/');
%  
% Path to the directory that holds all of the libraries
  LibrariesPath = '../LCP/Libraries/';
  addpath(cstrcat(LibrariesPath, 'Murty'));
%
% this code works with both Octave and Matlab.  Create the following 
% two constants. The two softwares have subtle differences.
  Octave = 1; Matlab = 2;
%
% select the software system here
  System = Octave;
%
  EPS1 = 1e-11;
%
% problem types, as described in the paper XXX.  Choose 1 or 2. (Problem 3
% is solved in the same manner as Problem 2.)

  Problem = 1; % Try solving the matrix equations for all stiffness matrices
  Problem = 2; % Treat the problem as a linear complementarity problem (LCP)
 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START providing data for the 3-disk problem
%
% number of unit cells: 
%   1 cell  =  9 particles
%   2 cells = 14 particles
%   3 cells = 19 particles
  Cells = 2;
%
% ellipse radius
  Radius = 1.00;  % radius of all spheres
%
% friction coefficient "mu", the normal stiffness "k", and the 
% tangential stiffness coefficient "alpha".  These are described in
% KPD Section 2.5
  mu = 0.50;
  k = 1.0;     % normal stiffness = k
  alpha = 1.0; % tangential stiffness = k * alpha
%
% normal contact force between spheres
  Fn = 0.001;
%
% the curvatures at the contacts will be multiplied by this factor
  Curvature = 0;
% Curvature = 2.00;
%
% vertical movement OR load of the cap
  Cap_Movement = -0.001; % active when LoadControl==0
  Cap_Load =     -0.001; % active when LoadControl==1
%
%
  Beta = 65;  %
  AspectRatio = sqrt(2) * tand(Beta);
%
% for AspectRatio = [0.75, 0.90, 1.00, 1.10, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00]
% for AspectRatio = [2.00 3.00]
  for AspectRatio = [2.00]
% for AspectRatio = [2*sqrt(2)]
%
% tangential contact force between spheres
  Ft = 0;
  Ft = mu*Fn;
%
% constraint conditions
  Turn = 0;      % whether base spheres can turn in unison about the x3 axis
  Tilt = 0;      % whether the cap is allowed to tilt
  Shift = 0;     % whether the cap is allowed to shift
  Twist = 0;     % whether the cap is allowed to twist
  Rotations = 3; % 0) all rotations disallowed, 
                 % 1) all can rotate, 
                 % 2) rotations of all exterior particles disallowed
                 % 3) rotations of exterior particles of cap & base disallowed
%
% whether a membrane is present, so that exterior forces depend on dx, and
% a H_g4 stiffness applies
  Membrane = 1;
%
% whether load control (1) or displacement control (0)
  LoadControl = 0;
%
% load direcction: positive 'P' (upward) or negative 'N' (downward)
  LD = 'P';
  if LoadControl==1 && Cap_Load<0
    LD = 'N';
  elseif LoadControl==0 && Cap_Movement<0
    LD = 'N';
  end
%
  Condition = 1;
%
% future contact parameters
  other_parameters = [];
%
  disp('Setting up the tetragonal BCC system...');
%
% create the input information for the system
  [ N, M, Unit, Pressure, Plate, q_Plate, ...
    x, V_pq, r_pq_p, r_pq_q, n_pq, f_pq, m_pq, rho_pq, K_pq, K_inv, ...
    mu_pq, k_pq, alpha_pq, FM_model_pq, idim, Beta] = ...
  Input_for_BCC( ...
    Radius, Fn, Ft, mu, k, alpha, Cells, AspectRatio, Condition, ...
    Curvature);
%
% report the Beta angle
  Beta
%
% initialize the increments of force and moment on the particles.  The
% vector [dp] is the stacked vector of force increments [db] and moment
% increments [dw], as in KPD Eq. 2_2.  Alternatively, with external follower
% forces, [dp] can represent the lower-dimensional increments in KPD Eq. 5_2
  dp = zeros(6*N,1);
%
  if Membrane==1
%   a triaxial chamber memberane is simulated.  dPdX is \partial p / \partial x
    [dPdX] = dpdx(N, Pressure, q_Plate, x);
    H_g4 = -dPdX;
  else
%   no external follower forces
    H_g4 = [];
  end
%
% type of constraints, as in KPD Section 2.4.  Note that in the current
% version, only Type 3 is accomodated
  Type_Constraint = 3;
%
% initialize some vectors associated with constraints
  c = []; f = []; dx_c = []; C = []; dc = []; R = []; Prr = []; Pnrr = [];
%
% label the base spheres
  Base = [1:4];
%
% label the cap spheres
  if Cells==1
    Cap = [6:9];
  elseif Cells==2
    Cap = [6:9] + 5;
  elseif Cells==3
    Cap = [6:9] + 10;
  end
%
  if Type_Constraint==3
%   displacment constraints
%
%   initialize C
    C = []; dc = [];
%
%   constraints at the base
    [C_sub, dc_sub] = Constraints_Base(N, Base, x, Turn);
    C = [C; C_sub];
    dc = [dc; dc_sub];
%
%   constraints at the cap
    if 1
      [ C_sub, dc_sub, dp] = ...
      Constraints_Cap2( ...
        N, Cap, x, Tilt, Shift, Twist, dp, ...
        LoadControl, Cap_Movement, Cap_Load);
    else
      [C_sub, dc_sub] = Constraints_Cap(N, Cap, Tilt, Shift, Cap_Movement);
    end
    C = [C; C_sub];
    dc = [dc; dc_sub];
%
%   prevent particle rotations, if Rotations=0
    [C_sub, dc_sub] = Constraints_NoRotations(N, Rotations, Cells);
    C = [C; C_sub];
    dc = [dc; dc_sub];
%
  else
    disp('NO CONSTRAINTS HAVE BEEN SPECIFIED')
  end
%
  disp('Finished setting up the tetragonal BCC system.');
%
% END providing data for the stacked-ellipse problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% list of the combinations of contact stiffness branches
% that should be investigated.  A zero "0" or empty vector
% [] means to investigate all possible branches.  A negative
% value means to investigate no branches.
  Combinations = [];
%
  Strain = [];
%
% do not wait until the end of run to report screen information
  more off; page_output_immediately(1)
%
% lambda_tolerance = an kx1 vector that gives the tolerances associated
% with the activation of possible slip conditions (lambdas).  For a 
% simple linear-frictional contact, k=1 and the lambda_tolerance is a scalar.
% For contact models with multiple slip mechanisms (e.g., both sliding and
% rotational friction), lambda_tolerance will be a vector
  lambda_tolerance = [0.0001];
%
% XXXXX - What are these numbers?
  Approx_Singular_Matrix = 0;
  Approx_Singular_Matrix = 1e7;
  Cond_numbers = 0;
%
  disp('Running LCP2_Stiffness_Pathologies_BCC. This can take a while...');
%
   [Details, Headings, ...
    M_hat, q_hat, H_x_BD_inv, dp, RCond_BD, ...
    C, P_L, P_Lperp, C_MP_inv, ...
    Pathologies, Available_Branch_Combinations, Condition_numbers, ...
    ...
    H_x, H_m_x, H_g1, H_g2, H_g3, H_g1234, B, ...
    H_lambda, Q_x, Q_lambda, M_lambda, ...
    map_lambda_to_pq_and_qp, map_pq_to_qp, ...
    ...
    H_m_pq_x, H_g1_pq, H_g2_pq, H_g3_pq, B_pq, H_m_max, ...
    H_m_pq_lambda, Q_pq_x_all, Q_pq_lambda_all, ...
    k_lambdas_pq_all, lambdas_pq_all ...
   ] = ...
  LCP2_Stiffness_Pathologies_BCC( ...
    idim, Problem, N, x, V_pq, r_pq_p, r_pq_q, ...
    n_pq, f_pq, m_pq, rho_pq, K_pq, K_inv, ...
    mu_pq, k_pq, alpha_pq, FM_model_pq, ...
    other_parameters, ...
    dp, Type_Constraint, c, f, dx_c, C, dc, ...
    H_g4, lambda_tolerance, Combinations, ...
    Approx_Singular_Matrix, Cond_numbers);
%
  disp('Finished running LCP2_Stiffness_Pathologies_BCC.');
%
  if Problem==2
%   solve the linear complementarity problem (LCP) for the lambda
%   values.  These lambda values will later be used for finding the dx
%   displacements.
%
%   initial guess, or enter the empty "[]" value, for zeros
    x0 = [];
%
%   bundles used in the spq function.  Not needed here
%   (see Instability/Octave/Simulation/sqp_bundles.m)
    PHI = [];
    Hc = [];
%
%   to reduce roundoff errors, we enforce symmetries in the M_hat matrix
%   (rotation about the cenral x3 axis)
    if 1 && Cells==2
      Permute = [2 3 4 1 6 7 8 5 10 11 12 9 14 15 16 13];
      Permute0 = [1:16];
      Permute1 = Permute0(Permute);
      Permute2 = Permute1(Permute);
      Permute3 = Permute2(Permute);
      M_hat = 0.25*(  M_hat(Permute0, Permute0) ...
                    + M_hat(Permute1, Permute1) ...
                    + M_hat(Permute2, Permute2) ...
                    + M_hat(Permute3, Permute3) ...
                   );
%
      if mean(abs(q_hat - mean(q_hat))) < 1e-16
        q_hat = mean(q_hat) * ones(size(q_hat));
      end
    end
%
    if 1
      disp('Running LCP_enumerate6.  This can take a while...');
      [ Lambda, Lambda_r, DegenSoln, IsolatedSoln, NonIsolated, ...
        NonIsolatedGrp, NonIsolatedDim, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate6(M_hat, q_hat, EPS1, 10*eps);
    elseif 0
      [ Lambda, Lambda_r, DegenSoln, IsolatedSoln, NonIsolated, ...
        NonIsolatedGrp, NonIsolatedDim, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate5(M_hat, q_hat, EPS1, 10*eps);
    elseif 0
      [ Lambda, Lambda_r, DegenSoln, IsolatedSoln, NonIsolated, ...
        NonIsolatedGrp, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate4(M_hat, q_hat, EPS1, 10*eps);
    elseif 0
      [ Lambda, Lambda_r, DegenSoln, IsolatedSoln, NonIsolated, ...
        NonIsolatedGrp, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate3(M_hat, q_hat, 100*eps);
    elseif 0
      [ Lambda, DegenSoln, IsolatedSoln, NonIsolated, NonIsolatedGrp, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate2(M_hat, q_hat, 10*eps);
    elseif 0
      [ Lambda, DegenSoln, IsolatedSoln, NonIsolated, ...
        ColCompetent, DegenMatrix, n_DegenMatrix] = ...
      LCP_enumerate(M_hat,q_hat);
    end
%  
%   number of active contacts
    M_lambda = length(q_hat);
%
%   number of solutions
    nLambda = size(Lambda, 2);
%
%   process the nLambda solutions, as many are duplicates.  Some of the
%   solutions are isolated, others are associated with multi-solution
%   non-isolated solutions.  We remove all of the duplicates, to arrive
%   at the list "Bif_all_list" of the solutions in "Lambda" that are 
%   legitimate
    [ Bif_I, Bif_II, Bif_III, Bif_II_Dim, Bif_III_Dim, ...
      Bif_I_Grp, Bif_II_Grp, Bif_III_Grp, Bif_I_X, ...
      Bif_II_Corners, Bif_III_Corners, ...
      Bif_I_list, Bif_II_III_list, Bif_all_list, Index_Type] = ...
    Count_Bifurcations( ...
      M_hat, q_hat, EPS1, ...
      Lambda, Lambda_r, DegenSoln, IsolatedSoln, NonIsolated, NonIsolatedGrp, ...
      NonIsolatedDim, ColCompetent, DegenMatrix, n_DegenMatrix);
%
%   the number of legitimate solutions
    nBif_all_list = length(Bif_all_list);
%
%   now that we have solved lambda, find the particle movements for the
%   subset of "legitimate" solutions
    if nLambda > 0
      dx = (  (C_MP_inv - H_x_BD_inv * H_x * C_MP_inv) * dc ...
            + H_x_BD_inv * dp ...
           ) * ones(1,nLambda)  ...
           - H_x_BD_inv * H_lambda * Lambda;
%  
%     reaction forces associated with constrained forces
      dy = - dp * ones(1,nLambda) ...
           + H_x*dx ...
           + H_lambda*Lambda;
%
%     the types of complementarity condition for each element if lambda
      if 0
        Limit1 = 1e-9;
        w = M_hat*Lambda + q_hat*ones(1,nLambda);
        Index_Type = 1000*ones(M_lambda,nLambda);
%
        Index_Type(find((Lambda > Limit1 & w > Limit1) ...
                        | ...
                        (Lambda< -Limit1& w< -Limit1))) = 0; % not a solution
        Index_Type(find(Lambda < Limit1 & w > Limit1)) =  1;
        Index_Type(find(Lambda > Limit1 & w < Limit1)) =  2;
        Index_Type(find(Lambda < Limit1 & w < Limit1)) =  3; % degen. solution
      end
%
%     du and dtheta
      du =     zeros(nLambda, N, 3);
      dtheta = zeros(nLambda, N, 3);
      db =     zeros(nLambda, N, 3);
      dw =     zeros(nLambda, N, 3);
%
      for i = 1:nLambda
        du(i,:,:) = [dx(1:3:3*N,i), dx(2:3:3*N,i), dx(3:3:3*N,i)];
        dtheta(i,:,:) = [dx(3*N+1:3:6*N,i), dx(3*N+2:3:6*N,i), dx(3*N+3:3:6*N,i)];
%
        db(i,:,:) = [dy(1:3:3*N,i),     dy(2:3:3*N,i),     dy(3:3:3*N,i)];
        dw(i,:,:) = [dy(3*N+1:3:6*N,i), dy(3*N+2:3:6*N,i), dy(3*N+3:3:6*N,i)];
      end
%
%     compute a parameter to characterize the stability of the current path.
%     The paramter is the second-order work in the direction of the current path
      if ~isempty(Lambda)
%       Path_Stability = dx' * H_x * dx ...
%                        + dx' * H_lambda * Lambda;
%       Path_Stability = 0.5 * sum(dy.*dx - dp.*dx, 1);
        Path_Stability = 0.5 * sum(dy.*dx - (-H_g4*dx).*dx, 1);
      else
        Path_Stability = [];
      end
%
%     determine the matrix for finding the "best fit" of the particle
%     displacements to an affine displacement field (strain field)
      [Affine_Matrix] = Affine_Fit(x);
%  
%     a consistent solution
      Affine = Affine_Matrix \ dx(1:3*N,:);
      Strain = Affine(4:12,:);  % arranged as 11, 12, 13, 21, 22, 23, 31, 32, 33
%
%     dilation (postive) values
      Dilate = sum(Strain([1,5,9],:),1);
    else
%     no solutions
      dx = [];
      dy = [];
      du = [];
      db = [];
      dw = [];
      Index_Type = [];
      dtheta = [];
      Path_Stability = [];
      Strain = [];
      Dilate = [];
    end
%
    disp('Finished running LCP_enumerate6');
%
H_xx_BD_inv = H_x_BD_inv;
H_x_x = H_x;
H_x_lambda = H_lambda;
H_lambda_x = Q_x;
H_lambda_lambda = Q_lambda;
    if 1
      disp('Running LCP_Instability2. This can take a while...');
%
      [ Lambda_Instab, NegEigenGrp, NegEigenValue, NegEigenVector] = ...
        LCP_Instability2(P_L, H_x, H_lambda, Q_x, Q_lambda, 10*eps);
%
      nLambda_Instab = size(Lambda_Instab,2);
%
%     du and dtheta
      du_Instab =     zeros(nLambda_Instab, N, 3);
      dtheta_Instab = zeros(nLambda_Instab, N, 3);
%
      for i = 1:nLambda_Instab
        du_Instab(i,:,:) = [NegEigenVector(1:3:3*N,i), ...
                            NegEigenVector(2:3:3*N,i), ...
                            NegEigenVector(3:3:3*N,i)];
        dtheta_Instab(i,:,:) = [NegEigenVector(3*N+1:3:6*N,i), ...
                                NegEigenVector(3*N+2:3:6*N,i), ...
                                NegEigenVector(3*N+3:3:6*N,i)];
      end
%
      disp('Finished running LCP_Instability2');
    else
      nLambda_Instab = [];
      Lambda_Instab = [];
      NegEigenGrp = [];
      NegEigenValue = [];
      NegEigenVector = [];
      du_Instab = [];
      dtheta_Instab = [];
    end
%
    if 1
      disp('Running LCP_Neutral2. This can take a while...');
      disp('Ignore errors: unable to recover undefined or non-optimal ...');
      [NeutralVec, NeutralGrp, NeutralRank] = ...
        LCP_Neutral2(M_hat, EPS1, 10*eps);
      nNeutralVec = size(NeutralVec,2);
%
      disp('Finished running LCP_Neutral2');
    elseif 0
      [ X_, Xr_, DegenSoln_, IsolatedSoln_, NonIsolated_, NonIsolatedGrp_, ...
        NonIsolatedDim_, ColCompetent_, DegenMatrix_, n_DegenMatrix_] = ...
      LCP_enumerate6(M_hat, zeros(size(M_hat,1), 1), EPS1, 10*eps);
%
      nNeutralVec = size(X_, 2) - 1;
      NeutralVec = X_;
      NeutralGrp = NonIsolatedGrp_;
      NeutralRank = NonIsolatedDim_;
    else
      nNeutralVec = [];
      NeutralVec = [];
      NeutralGrp = [];
      NeutralRank = [];
    end
%
    if 1
      disp('Running LCP_Sensitive. This can take a while...');
%
      [Sensitive] = LCP_Sensitive(M_hat, q_hat, Lambda, 10*eps);
%
      disp('Finished running LCP_Sensitive');
    else
      Sensitive = [];
    end
%
    if 1
      disp('Running LCP_enumerate6. This can take a while...');
%
      [ Lambda_, Lambda_r_, DegenSoln_, IsolatedSoln_, ...
        NonIsolated_, NonIsolatedGrp_, NonIsolatedDim_, ...
        ColCompetent_, DegenMatrix_, n_DegenMatrix_] = ...
      LCP_enumerate6(M_hat, ones(size(M_hat,1),1), EPS1, 10*eps);
%
      is_R_matrix = isempty(Lambda_) || all(all(Lambda_ < 1000*eps));
      is_P_matrix = isPmatrix(M_hat);
      is_PD_matrix = all(real(eig(M_hat)) > 10*eps);
      is_R0_matrix = ~nNeutralVec;
%
      disp('Finished running LCP_enumerate6.');
    end
%
    Details =  cell(nLambda,66);
    Headings = cell(66,1);
%
    Details{ 1} = M_lambda;        %  1) M_lambda
    Details{ 2} = M_hat;           %  2) M_hat
    Details{ 3} = q_hat;           %  3) q_hat
    Details{ 4} = nLambda;         %  4) nLambda = number of solutions
    Details{ 5} = Index_Type;      %  5) Index_Type
    Details{ 6} = RCond_BD;        %  6) condition no. of Bott-Duffin inverse
    Details{ 7} = DegenSoln;       %  7) whether a degenerate solution
    Details{ 8} = IsolatedSoln;    %  8) whether an isolated solution
    Details{ 9} = NonIsolated;     %  9) whether a nonisolated solution
    Details(10) = NonIsolatedGrp;  % 10) r for NonIsolatedGrp
    Details{11} = ColCompetent;    % 11) whether a column competent matrix
    Details{12} = DegenMatrix;     % 12) whether a degenerate matrix
    Details{13} = n_DegenMatrix;   % 13) type of degenerate matrix
    Details{14} = dx;              % 14) vectors of dx of solutions
    Details{15} = du;              % 15) vectors of du, 3N x 3 of solutions
    Details{16} = dtheta;          % 16) vectors of dtheta, 3N x 3 of solutions
    Details{17} = dy;              % 17) vectors of dy of solutions
    Details{18} = db;              % 18) vectors of db, 3N x 3 of solutions
    Details{19} = dw;              % 19) vectors of dw, 3N x 3 of solutions
    Details{20} = Path_Stability;  % 20) stability parameter of solutions
    Details{21} = C_MP_inv;        % 21) C_MP_inv
    Details{22} = H_x_BD_inv;      % 22) H_x_BD_inv
    Details{23} = H_x;             % 23) H_x
    Details{24} = H_lambda;        % 24) H_lambda
    Details{25} = nLambda_Instab;  % 25) number of instability modes
    Details{26} = Lambda_Instab;   % 26) Lambda of instability modes
    Details{27} = NegEigenGrp;     % 27) enumeration group, r, of inst. mode
    Details{28} = NegEigenValue;   % 28) neg. eigenvalues of inst. modes
    Details{29} = NegEigenVector;  % 29) neg. eigenvectors of inst. modes
    Details{30} = du_Instab;       % 30) k x N x 3 du's of inst. modes
    Details{31} = dtheta_Instab;   % 31) k x N x 3 dtheta's of inst. modes
    Details{32} = nNeutralVec;     % 32) no. of neutral equilibrium (NE) modes
    Details{33} = NeutralVec;      % 33) neutral equilibrium (NE) vectors
    Details{34} = NeutralGrp;      % 34) enumeration group, r, of (NE) vector
    Details{35} = NeutralRank;     % 35) rank of null space in this r group
    Details{36} = Sensitive;       % 36) whether soln. lambda is path-sensitive
    Details{37} = Strain;          % 37) Strains: 11,12,13,21,22,23,31,32,33
    Details{38} = Dilate;          % 38) Dilation: 11+22+33
    Details{39} = Lambda;          % 39) all lambda solutions
    Details{40} = Lambda_r;        % 40) r combinations of lambda solutions
    Details{41} = DegenSoln;       % from LCP_enumerate6
    Details{42} = IsolatedSoln;
    Details{43} = NonIsolated;
    Details{44} = NonIsolatedGrp;
    Details{45} = NonIsolatedDim;
    Details{46} = ColCompetent;
    Details{47} = DegenMatrix;
    Details{48} = n_DegenMatrix;
    Details{49} = Bif_I;           % from Count_Bifurcations
    Details{50} = Bif_II;
    Details{51} = Bif_III;
    Details{52} = Bif_II_Dim;
    Details{53} = Bif_III_Dim;
    Details{54} = Bif_I_Grp;
    Details{55} = Bif_II_Grp;
    Details{56} = Bif_III_Grp;
    Details{57} = Bif_I_X;
    Details{58} = Bif_II_Corners;
    Details{59} = Bif_III_Corners;
    Details{60} = Bif_I_list;
    Details{61} = Bif_II_III_list;
    Details{62} = Bif_all_list;
    Details{63} = is_R_matrix;     % from herein
    Details{64} = is_P_matrix;
    Details{65} = is_PD_matrix;
    Details{66} = is_R0_matrix;
%
    Headings(:) = ...
      {'M_lambda', ...
       'M_hat', ...
       'q_hat', ...
       'nLambda', ...
       'Index_Type', ...
       'RCond_BD', ...
       'DegenSoln', ...
       'IsolatedSoln', ...
       'NonIsolated', ...
       'NonIsolatedGrp', ...
       'ColCompetent', ...
       'DegenMatrix', ...
       'n_DegenMatrix', ...
       'dx', ...
       'du', ...
       'dtheta', ...
       'dy', ...
       'db', ...
       'dw', ...
       'Path_Stability', ...
       'C_MP_inv', ...
       'H_x_BD_inv', ...
       'H_x', ...
       'H_lambda', ...
       'nLambda_Instab', ...
       'Lambda_Instab', ...
       'NegEigenGrp', ...
       'NegEigenValue', ...
       'NegEigenVector', ...
       'du_Instab', ...
       'dtheta_Instab', ...
       'nNeutralVec', ...
       'NeutralVec', ...
       'NeutralGrp', ...
       'NeutralRank', ...
       'Sensitive', ...
       'Strain', ...
       'Dilate', ...
       'Lambda', ...
       'Lambda_r', ...
       'DegenSoln', ...
       'IsolatedSoln', ...
       'NonIsolated', ...
       'NonIsolatedGrp', ...
       'NonIsolatedDim', ...
       'ColCompetent', ...
       'DegenMatrix', ...
       'n_DegenMatrix', ...
       'Bif_I', ...
       'Bif_II', ...
       'Bif_III', ...
       'Bif_II_Dim', ...
       'Bif_III_Dim', ...
       'Bif_I_Grp', ...
       'Bif_II_Grp', ...
       'Bif_III_Grp', ...
       'Bif_I_X', ...
       'Bif_II_Corners', ...
       'Bif_III_Corners', ...
       'Bif_I_list', ...
       'Bif_II_III_list', ...
       'Bif_all_list', ...
       'is_R_matrix', ...
       'is_P_matrix', ...
       'is_PD_matrix', ...
       'is_R0_matrix'};
%
  elseif Problem==1
%
%
%   if the list of pathologies is small, print to the screen
%   if size(Pathologies,1) < 40
      printf(cstrcat(' %0.1f %6i %+12.3e %2i  ', ...
                     '%+0.6f %+0.6f %+0.6f %+0.6f\n'), ...
                     Pathologies(:,[1:4,[1,5,9,3]+4])')
%   end
%
  end  % Problem==1
%
% archive the results
%
  if Curvature <= 0
    File_Name = sprintf(cstrcat( ...
             'Prob_%1i_Cells_%1i_Asp_%1.3f_Turn_%1i_Tilt_%1i', 
             '_Shift_%1i_Twist_%1i_Rot_%1i_Memb_%1i_LC_%1i_LD_%1c_Fn_%7.1e', ...
             '_FtFn_%1.2f'), ...
             Problem, Cells, AspectRatio, Turn, Tilt, Shift, Twist, ...
             Rotations, Membrane, LoadControl, LD, ...
             Fn, Ft/Fn);
  else
    File_Name = sprintf(cstrcat( ...
             'Prob_%1i_Cells_%1i_Asp_%1.3f_Turn_%1i_Tilt_%1i', 
             '_Shift_%1i_Twist_%1i_Rot_%1i_Memb_%1i_LC_%1i_LD_%1c_Fn_%7.1e', ...
             '_FtFn_%1.2f_Curv_%1.2f'), ...
             Problem, Cells, AspectRatio, Turn, Tilt, Shift, Twist, ...
             Rotations, Membrane, LoadControl, LD, ...
             Fn, Ft/Fn, Curvature);
  end
%
  File_Path = cstrcat('Results/', File_Name, '.mat');
  save -binary Temp_Details ...
       Problem Cells AspectRatio Beta Turn Tilt Shift Rotations ...
       Membrane LoadControl Radius mu k alpha Ft Fn Cap_Movement ...
       Cap_Load P_L H_x H_lambda Q_x Q_lambda ...
       Details Headings
  system(cstrcat('mv Temp_Details ',File_Path));
  printf(cstrcat('Results saved as ', File_Path, '\n'));
%
  end
%
% [Min, Likely] = min(Path_Stability);
% Likely_Dilate = Dilate(Likely)
% Force_Likely = sum(dy(3*([11:14] - 1)+3, Likely))
