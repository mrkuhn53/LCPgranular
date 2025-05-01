function [Shape, idim, xcell, xcell_init, N, Size, u, Qp, ...
          satrad, cenrad, cirrad, njoin, s_rad, b_rad, xlocal, ...
          rmax, ravg, periodic_directions, Contact_List, ...
          def, istep, V_pq, C_pq, f_pq, m_pq, df_pq, dm_pq, n_pq, r_pq_p, ...
          Overlap_pq, k_pq, Map_p_Rel_to_All, iupdat, I, J, ...
          Time, istepo, icont, iconto, init_output, Output_Cell_Array, ...
          Output_Instab_Array, ...
          M_lambda, M_hat, q_hat, lambda, last_Save_Time, last_Write_Time, ...
          Method_used, is_Soln, Lambda, DegenSoln, DegenMatrix, ...
          du, b_p_Rel, RCond_BD, iOutput, File_Path, File_Path_Instab, ...
          Headings, Headings_Instab, ...
          Repetitions, Time_to_Stop, Path_Stability, ...
          R, Prr, Pnrr, P_L, P_Lperp, C_MP_inv, dTime, stress, ...
          H_x, H_lambda, Q_x, Q_lambda, H_stress_x, H_stress_lambda, ...
          C_wo_dxcell, ...
          H_x_x, H_x_lambda, H_lambda_x, H_lambda_lambda, ...
          map_lambda_to_pq_and_qp, u_update, xcell_update,  ...
          Threshold_separation, nbins, bin_list, surround_bins, dgamma, ...
          File_Path_alt, Output_Cell_Array_alt, Headings_alt, M_alt,q_alt] = ...
         Create_Assembly_D(Directory, File, scale_sizes)
%
% initialize Contact_List, since no contacts are yet known
  Contact_List = [];
%
  [Shape, xcell, np, ...
   HalfWidth, Aspect, x, theta, gamma, Beta_rad, Beta, ...
   nobs, nbumps, satrad, cenrad, cirrad, Qp] =  ...
  Read_D_file3(Directory, File);
%
% sizes of the particle (e.g., radii)
  Size = HalfWidth;
%
% basic particle shapes
  [Circle, Oval, Ellipse, Sphere, Ovoid, Nobby, Bumpy] = Shapes3();
  if Shape==Circle || Shape==Nobby
    idim = 2;
  elseif Shape==Sphere || Shape==Bumpy
    idim = 3;
  else
    disp('Only disks, spheres, nobbies, or bumpies allowed.')
    ERROR_Create_Assembly_D_1
  end
%
% number of particles, using "N" instead of "np"
  N = np;
%
% particle positions, using "u" instead of "x"
  u = x;
  if size(u,2) < 3
    u = [u, zeros(size(u,1),3-size(u,2))];
  end
%
% the code currently handles four shapes: circles, spheres, nobbies,
% and bumpies
%
% put the particles back inside of the periodic boundaries
  du = zeros(size(u));           % zero initial displacement
  dxcell = zeros(size(xcell));   % zero initial cell deformation
%
% periodic directions
  if idim==2
    periodic_directions = [1 1 0];
  else
    periodic_directions = [1 1 1];
  end
%
% put those particles back inside the periodic cell
  [u, du] = PutBack(u, du, xcell, dxcell, periodic_directions);
%
% numbers of circles/spheres joined together with each particle
  if Shape==Circle || Shape==Sphere
    njoin = 1;
%
    s_rad(1) = 1;
    b_rad(1) = 0;
%
    xlocal = zeros(3,njoin);
%
    satrad = 0;
    cenrad = 1;
    cirrad = 0;
%
%   quaternion orientation of the particle
    Qp = ones(N,1)*[cos(0/2), sin(0/2)*[0 0 1]];
%
%   maximum and mean particle radius, as used to find nearby particle-pairs
    rmax = max(Size);
    ravg = mean(Size);
%
  elseif Shape==Nobby
%   nobby particles - a composite of circles
%
%   A nobby particle is a non-convex shape that is the union of a set
%   of circles.  The set includes a central circle and several
%   satellite circles.  The number of satellite circles is the input
%   variable "nobs", so that the total number of circles is nobs+1.
%   The satellite circles are evenly spaced, with the centers lying
%   on a circle.  This circle has input radius "rad".  The satellite
%   circles have radius satrad*rad.  The central circle has radius
%   cenrad*rad.
%
    njoin = nobs + 1;
%
%   for nobbies, there is not an independent cirrad
    cirrad = 1;
%
    s_rad = zeros(njoin,1);
    b_rad = zeros(njoin,1);
%
%   the central sphere
    s_rad(1) = cenrad;
    b_rad(1) = 0;
%
%   the satellite spheres
    s_rad(2:njoin) = satrad;
    b_rad(2:njoin) = cirrad;
%
%   initialize relative positions of all surrounding (satellite) circles
    xlocal = zeros(3,njoin);
%
%   relative locations of satellite circles
    for i = 1:nobs
      xlocal(1,i+1) = cos((i-1)*2*pi/nobs);
      xlocal(2,i+1) = sin((i-1)*2*pi/nobs);
    end
%
%   convert degrees to radians.  "theta" is the angle, measured counterclockwise
%   from the x_1 axis to the center of the first circular "nob"
    theta = [zeros(N,2), (pi/180)*theta];
%
%   quaternion orientation of the particle, using the "theta" angle
%   for i = 1:N
%     Qp(i,:) = [cos(theta(i,3)/2), sin(theta(i,3)/2)*[0, 0, 1]];
%   end
    Qp = [cos(theta(:,3)/2), sin(theta(:,3)/2)*[0, 0, 1]];
%
%   maximum and mean particle radius, as used to find nearby particle-pairs
    rmax = max(cenrad, cirrad+satrad)*max(Size);
    ravg = max(cenrad, cirrad+satrad)*mean(Size);
%
  elseif Shape==Bumpy
%   A bumpy particle is an composite of multiple component spheres.
%   These component spheres include a central sphere and several
%   satellite spheres, having the following characteristics
%     cenrad  =  The radius of the central sphere, for a particle of
%                unit rad(igrain).  That is the radius of the
%                central sphere is cenrad*rad(igrain)
%     satrad  =  The radius of the satellite spheres, for a particle of
%                unit rad(igrain).  That is, the radius of the
%                satellite spheres is satrad*rad(igrain)
%     cirrad  =  The radius of the sphere on which the satellite
%                spheres are centered, for a particle of
%                unit rad(igrain).  That is, cirrad*rad(igrain) is 
%                the radius of the circumsphere of satellite centers.
%     nbumps  =  The number of satellite spheres, clustered around
%                the central sphere
%
%   All dimensions in this section of code are for a particle
%   of unit input radius (that is, assuming that input "rad" 
%   values are 1 in the StartFile).  Other "rad" values will later
%   be scaled accordingly.
%
%   The component spheres will be numbered as follows: the central
%   sphere will have index "0", and the satellite spheres will be
%   numbered from "1" through "nbumps"
%
%   Create lists of the characteristics of the component spheres:
%     s_rad(0:nbumps) is the radius of each component sphere.
%                     = cenrad for the central partcles
%                     = cirrad for the satellite particles
%     b_rad(0:nbumps) is the distance (radius) from the center of the
%                        particle to the center of the component sphere.
%                     = 0 for the central particle
%                     = cirrad for the satellite spheres
%   Later, these values will be multiplied by the "rad" values obtained
%   from the input StartFile
%
%   the number of spheres in the cluster
    njoin = nbumps + 1;
%
    s_rad = zeros(njoin,1);
    b_rad = zeros(njoin,1);
%
%   the central sphere
    s_rad(1) = cenrad;
    b_rad(1) = 0;
%
%   the satellite spheres
    s_rad(2:njoin) = satrad;
    b_rad(2:njoin) = cirrad;
%
%   initialize relative positions of all surrounding (satellite) spheres
    xlocal = zeros(3,njoin);
%
%   maximum and mean particle radius, as used to find nearby particle-pairs
    rmax = max(cenrad, cirrad+satrad)*max(Size);
    ravg = max(cenrad, cirrad+satrad)*mean(Size);
%
    if nbumps==0
%     A single (central) sphere with no satellite spheres.  We include
%     this case as a means of testing the bumpy code
    elseif nbumps==2
%     A single (central) sphere with two satellite spheres symmetrically
%     place on the x-x axis of the central sphere
%
%     The first satellite sphere is centered on the x-axis.
      xlocal(1,2) = 1;
%     The second satellite sphere is centered on the x-axis.
      xlocal(1,3) = -1;
%
    elseif nbumps==3
%     Three satellite spheres on the y-z plane, centered around 
%     the x-x axis.
%
%     Locations of the vertices on a circumsphere of unit radius
      xlocal(2,2) = 1;
%
      xlocal(2,3) = -0.5;
      xlocal(3,3) = sqrt(3)/2;
%
      xlocal(2,4) = -0.5;
      xlocal(3,4) = -sqrt(3)/2;
%
    elseif nbumps==4
%     Four satellite spheres centered on the vertices of a tetrahedron 
%
%     Locations of the vertices on a circumsphere of unit radius
%
%     The first satellite sphere is centered on the z-axis.  All others
%     are located below the x-y plane.
      xlocal(3,2) = 1;
%
%     The other three satellite spheres are centered on the z = -1/3 
%     plane
      xlocal(2,3) = sqrt(8 / 9);
      xlocal(3,3) = -1 / 3;
%
      xlocal(1,4) = -sqrt(2 / 3);
      xlocal(2,4) = -sqrt(2 / 9);
      xlocal(3,4) = -1 / 3;
%
      xlocal(1,5) =  sqrt(2 / 3);
      xlocal(2,5) = -sqrt(2 / 9);
      xlocal(3,5) = -1 / 3;
%
    elseif nbumps==6
%     Six satellite spheres centered on the vertices of an octohedron 
%
%     Locations of the vertices on a circumsphere of unit radius
%
%     The first and second satellite spheres are centered on the z-axis.
%     All others are centered on the x-y plane.
      xlocal(3,2) = 1;
%
      xlocal(3,3) = -1;
%
%     The other four satellite spheres are centered on the z = 0
%     plane, arranged counterclockwise around the z-axis
      xlocal(1,4) = 1;
      xlocal(2,5) = 1;
      xlocal(1,6) = -1;
      xlocal(2,7) = -1;
    elseif nbumps==8
%
%     Eight satellite spheres centered on the vertices of a cube 
%
%     Locations of the vertices on a circumsphere of unit radius
%
%     The first through fourth satellite spheres are centered 
%     centered above the x-y plane, arranged counterclockwise 
%     around the z-axis.  All others are centered below the x-y plane.
%
      sqrt13 = sqrt(1 / 3);
      xlocal(1,2) = sqrt13;
      xlocal(2,2) = sqrt13;
      xlocal(3,2) = sqrt13;
%
      xlocal(1,3) = -sqrt13;
      xlocal(2,3) = sqrt13;
      xlocal(3,3) = sqrt13;
%
      xlocal(1,4) = -sqrt13;
      xlocal(2,4) = -sqrt13;
      xlocal(3,4) = sqrt13;
%
      xlocal(1,5) = sqrt13;
      xlocal(2,5) = -sqrt13;
      xlocal(3,5) = sqrt13;
%
%     The fifth through eighth satellite spheres are centered 
%     centered below the x-y plane, arranged counterclockwise 
%     around the z-axis.
%
      xlocal(1,6) = sqrt13;
      xlocal(2,6) = sqrt13;
      xlocal(3,6) = -sqrt13;
%
      xlocal(1,7) = -sqrt13;
      xlocal(2,7) = sqrt13;
      xlocal(3,7) = -sqrt13;
%
      xlocal(1,8) = -sqrt13;
      xlocal(2,8) = -sqrt13;
      xlocal(3,8) = -sqrt13;
%
      xlocal(1,9) = sqrt13;
      xlocal(2,9) = -sqrt13;
      xlocal(3,9) = -sqrt13;
%
    else
      disp('Error in Create_Assembly_D.m: nbumps must be 2, 3, 4, 6, 8.')
      ERROR_Create_Assembly_D_2
    end
  end
%
% particle orientation.  Note that quaternions will be used to handle the
% particles' rotations from the initial orientation "theta"
  theta = zeros(size(u));
%
  if isempty(scale_sizes) || scale_sizes==0
%   do not scale the assembly size
%
  elseif scale_sizes > 0
    Factor = scale_sizes;
%
%   apply the scaling factor
    xcell = Factor * xcell;
    Size =  Factor * Size;
    u =     Factor * u;
    rmax =  Factor * rmax;
    ravg =  Factor * ravg;
%
  elseif scale_sizes < 0
    Factor = abs(scale_sizes) / mean(Size);
%
%   apply the scaling factor
    xcell = Factor * xcell;
    Size =  Factor * Size;
    u =     Factor * u;
    rmax =  Factor * rmax;
    ravg =  Factor * ravg;
  end
%
% initial cell dimensions
  xcell_init = xcell;
%
% initialize the control step
  istep = 0;
%
% initialize the deformation gradient
  def = eye(3);
%
% initialize V_pq, etc.
  V_pq = [];
  C_pq = [];
  f_pq = [];
  m_pq = [];
  df_pq = [];
  dm_pq = [];
  n_pq = [];
  r_pq_p = [];
  Overlap_pq = [];
  k_pq = [];
  V_pq_old = [];
  C_pq_old = [];
%
  iupdat = 0;
%
  Time = 0;
  istepo = 0;
  init_output = 1;
  iconto = zeros(6,1);
  Output_Cell_Array = [];
  Output_Instab_Array = [];
  Headings = [];
  Headings_Instab = [];
  M_lambda = 0;
  M_hat = [];
  q_hat = [];
  lambda = [];
  icont = [];
  last_Save_Time = 0;
  last_Write_Time = 0;
  iOutput = 0;
  File_Path = '';
  File_Path_Instab = '';
  RCond_BD = 0;
  du = [0 0 0];
  Map_p_Rel_to_All = [];
  Method_used = 0
  is_Soln = 1;
  Lambda = [];
  DegenSoln = [];
  DegenMatrix = [];
  b_p_Rel = [];
  Repetitions = 0;
  Time_to_Stop = 0;
  Path_Stability = 0;
  dTime = 1;
  stress = zeros(3);
  map_lambda_to_pq_and_qp = [];
  u_update = [];
  xcell_update = [];
  Threshold_separation = [];
  nbins = [];
  bin_list = [];
  surround_bins = [];
  H_x = [];
  H_lambda = [];
  Q_x = [];
  Q_lambda = [];
  H_stress_x = [];
  H_stress_lambda = [];
  dgamma = [];
  H_x_x = [];
  H_x_lambda = [];
  H_lambda_x = [];
  H_lambda_lambda = [];
  C_wo_dxcell = [];
  File_Path_alt = [];
  Output_Cell_Array_alt = [];
  Headings_alt = [];
  M_alt = [];
  q_alt = [];
%
% projection matrices for constraint Types II, II, and IV.
% First, initialize them
  R = []; Prr = []; Pnrr = []; P_L = []; P_Lperp = []; C_MP_inv = [];
%
% the I and J components of stress
  I = zeros(1,6);
  J = zeros(1,6);
  I(1) = 1; I(2) = 2; I(3) = 3; I(4) = 1; I(5) = 1; I(6) = 2;
  J(1) = 1; J(2) = 2; J(3) = 3; J(4) = 2; J(5) = 3; J(6) = 3;
