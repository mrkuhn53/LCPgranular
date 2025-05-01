%
% |+|+|+|+|+|+|+|+|+|  Revision 1   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function [ N, M, Unit, Pressure, Plate, q_Plate, ...
           x, V_pq, r_pq_p, r_pq_q, n_pq, f_pq, m_pq, rho_pq, K_pq, K_inv, ...
           mu_pq, k_pq, alpha_pq, FM_model_pq, idim, Beta] = ...
         Input_for_BCC( ...
           Radius, Fn, Ft, mu, k, alpha, Cells, AspectRatio, Condition, ...
           Curvature);
%
% This function creates the input for a stack of ellipses.  The stack lies 
% within the x_1 - x_2 plane
%
% ------------ INPUT ----------------------------------------------------
% Radius1 =   radius of top and bottom disks
% Radius2 =   radius of middle disk
% Fn =        current normal contact force at the disk-disk contacts
% mu =        friction coefficient, disk-disk contacts
% k =         normal stiffness, disk-disk contacts
% alpha =     coefficient for tangential stiffness, disk-disk contacts
% Beta =      angle of the lower pair from horizontal
% CurveFactor=increase the cuvatures of the disk-disk contacts by this
%             factor.  If factor is [], 0, 1, or negative, use 1
%
% -------------- OUTPUT -------------------------------------------------
% N =         number of particles
% M =         number of contacts.  Note that a pair of touching particles has
%             two contacts(!): a pq variant and a qp variant.  Therefore,
%             M is twice the number of actual contacts (i.e. twice the number
%             of touching pairs of particles)
% x =         locations of the centers of the particles, x(N,3)
% V_pq =      M x 2 matrix that gives the assembly topology, where M is twice
%             the number of contacts.  V_pq(pq,1) is the "p" particle of 
%             "pq" contact; whereas, V_pq(c,2) is the "q" particle of the 
%             "pq" contact.
% r_pq_p =    M x 3 array. Row r_pq_p(pq,:) is the vector from particle p
%             to contact pq
% r_pq_q =    M x 3 array. Row r_pq_q(pq,:) is the vector from particle q
%             to contact pq
% n_pq =      unit normal vector of contact pq, directed outward from particle 
%             p.  That is, outward from particle V_pq(pq,1).  An M x 3 array.
% f_pq =      force at contact pq, acting upon particle p.  That is, acting
%             upon particle V_pq(pq,1).  An M x 3 array, where M is twice
%             the number of contacts.
% rho_pq =    the radius of curvature of particles p and q at their contact pq.
%             An M x 2 array
% K_pq =      curvature tensors of the two particles, p and q, at their 
%             contact pq. It is an M x 2 x 3 x 3 array (a 4D array), where 
%             K_pq(pq,1,:,:) is the 3x3 tensor for particle p of contact pq,
%             and K_pq(pq,2,:,:) is the 3x3 tensor for particle q of 
%             contact pq
% K_inv =     the pseudo-inverse of [K_p + K_q].  An Mx3x3 array.
% mu_pq =     M x 1 array of the friction coefficients of the contacts
% k_pq =      M x 1 array of the normal stiffness of the contacts
% alpha_pq =  M x 1 array of the tangential stiffness coefficients
% FM_model_pq =M x 1 array.  The contact model for the force of each contact.
%             Only one model is coded in function F_M_pq.m: 
%               = 1, the linear-frictional model of Section 2.5, with zero
%                 contact moments
%
% a 3D problem
  idim = 3;
%
%---- create the stack of N_Disks disks with centers in matrix x(N,3) -----
%
  if Condition==1
%   a stack of equal-size ellipses on a half-ellipse base under 
%   self-weight conditions.  The first particle (n=1) is the base:
%   a half-ellipse, and the particles are numbered upward
%
%   number of unit cells: 
%     1 cell  =  9 particles
%     2 cells = 14 particles
%     3 cells = 19 particles
%
    if Cells==1
      N = 9;
    elseif Cells==2
      N = 14;
    elseif Cells==3
      N = 19;
    else
      disp('ERROR: invalid value of Cells')
      ERROR
    end
%
%   check for a height-to-width ratio that is within permissible range
    if AspectRatio > sqrt(2) || AspectRatio < 1/sqrt(2)
      disp('Warning: AspectRatio requires truncated spheres')
    end
%
%   unit cell size in 3 directions
    Unit = (4*Radius / sqrt(2 + AspectRatio^2)) * [1, 1, AspectRatio];
    hUnit = 0.5 * Unit;
    Unit2 = 2 * Unit;
    Unit3 = 3 * Unit;
    Unit4 = 4 * Unit;
    Unit5 = 5 * Unit;
    Unit6 = 6 * Unit;
%
%   positions of the 14 particles
    x( 1,:) = [ hUnit(1),  hUnit(2),     0];
    x( 2,:) = [-hUnit(1),  hUnit(2),     0];
    x( 3,:) = [-hUnit(1), -hUnit(2),     0];
    x( 4,:) = [ hUnit(1), -hUnit(2),     0];
%
    x( 5,:) = [    0,     0, hUnit(3)];
%
    x( 6,:) = [ hUnit(1),  hUnit(2), Unit(3)];
    x( 7,:) = [-hUnit(1),  hUnit(2), Unit(3)];
    x( 8,:) = [-hUnit(1), -hUnit(2), Unit(3)];
    x( 9,:) = [ hUnit(1), -hUnit(2), Unit(3)];
%
    if Cells >=2
      x(10,:) = [    0,     0, Unit(3)+hUnit(3)];
%
      x(11,:) = [ hUnit(1),  hUnit(2), Unit2(3)];
      x(12,:) = [-hUnit(1),  hUnit(2), Unit2(3)];
      x(13,:) = [-hUnit(1), -hUnit(2), Unit2(3)];
      x(14,:) = [ hUnit(1), -hUnit(2), Unit2(3)];
%
      if Cells ==3
        x(15,:) = [    0,     0, Unit2(3)+hUnit(3)];
%
        x(16,:) = [ hUnit(1),  hUnit(2), Unit3(3)];
        x(17,:) = [-hUnit(1),  hUnit(2), Unit3(3)];
        x(18,:) = [-hUnit(1), -hUnit(2), Unit3(3)];
        x(19,:) = [ hUnit(1), -hUnit(2), Unit3(3)];
      end
    end
%
%   the angle from the horizontal of diagonally opposed particles of
%   levels 1 and 2
    Beta = atand(x(8,3)/sqrt(sum((x(1,1:2)-x(8,1:2)).^2,2)));
%
%-----identify the contacts among the disks and build the "v" matrix --------
%
%   number of contacts. There are 2 contacts per pair of touching particles,
%   since each contact has a pq and qp variant
    M = 16 * Cells;
%
%   the contacts are as follows:
%     contact 1 = between the half-ellipse base (1) and first ellipse (2)
%     contact 2 = between first stacked ellipse (2) and the half-ellipse base (1)
%     contact 3 = between first stacked ellipse (2) and second ellipse (3)
%     contact 4 = between second ellipse (3) and first stacked ellipse (2)
%     etc.
%
%   initialize matrix "V_pq", which is a collapsed incidence matrix that gives
%   the topology of the particle assembly.  V_pq(pq,1) is the "p" particle of
%   "pq" contact; whereas, V_pq(c,2) is the "q" particle of the "pq" contact.
%   If M is the total number of contacts, then "V_pq" will have size 2Mx2,
%   since pq and qp are treated as separate entries.  Now, initialize the
%   matrix "V_pq"
    V_pq = zeros(M,2);
    V_pq( 1,:) = [ 1,  5];
    V_pq( 2,:) = [ 5,  1];
    V_pq( 3,:) = [ 2,  5];
    V_pq( 4,:) = [ 5,  2];
    V_pq( 5,:) = [ 3,  5];
    V_pq( 6,:) = [ 5,  3];
    V_pq( 7,:) = [ 4,  5];
    V_pq( 8,:) = [ 5,  4];
%
    V_pq( 9,:) = [ 5,  6];
    V_pq(10,:) = [ 6,  5];
    V_pq(11,:) = [ 5,  7];
    V_pq(12,:) = [ 7,  5];
    V_pq(13,:) = [ 5,  8];
    V_pq(14,:) = [ 8,  5];
    V_pq(15,:) = [ 5,  9];
    V_pq(16,:) = [ 9,  5];
%
    if Cells >=2
      V_pq(17,:) = [ 6, 10];
      V_pq(18,:) = [10,  6];
      V_pq(19,:) = [ 7, 10];
      V_pq(20,:) = [10,  7];
      V_pq(21,:) = [ 8, 10];
      V_pq(22,:) = [10,  8];
      V_pq(23,:) = [ 9, 10];
      V_pq(24,:) = [10,  9];
%
      V_pq(25,:) = [10, 11];
      V_pq(26,:) = [11, 10];
      V_pq(27,:) = [10, 12];
      V_pq(28,:) = [12, 10];
      V_pq(29,:) = [10, 13];
      V_pq(30,:) = [13, 10];
      V_pq(31,:) = [10, 14];
      V_pq(32,:) = [14, 10];
%
      if Cells ==3
        V_pq(33,:) = [11, 15];
        V_pq(34,:) = [15, 11];
        V_pq(35,:) = [12, 15];
        V_pq(36,:) = [15, 12];
        V_pq(37,:) = [13, 15];
        V_pq(38,:) = [15, 13];
        V_pq(39,:) = [14, 15];
        V_pq(40,:) = [15, 14];
%
        V_pq(41,:) = [15, 16];
        V_pq(42,:) = [16, 15];
        V_pq(43,:) = [15, 17];
        V_pq(44,:) = [17, 15];
        V_pq(45,:) = [15, 18];
        V_pq(46,:) = [18, 15];
        V_pq(47,:) = [15, 19];
        V_pq(48,:) = [19, 15];
      end
    end
%
%----- contact properties of the disks: contact vectors, --------
%      normal vectors, and contact force vectors.
%
%   the branch vectors for the contacts, from p to q
    rx = zeros(M,3);
    rx = x(V_pq(:,2),:) - x(V_pq(:,1),:);
%
%   impose uniformity among the rx values, to improve accuracy
    rx(:,1) = sign(rx(:,1)) .* hUnit(1);
    rx(:,2) = sign(rx(:,2)) .* hUnit(2);
    rx(:,3) = sign(rx(:,3)) .* hUnit(3);
%
%   the normal vectors of the contacts, outward from p
    n_pq = rx ./ (sqrt(sum(rx.^2, 2)) * [1 1 1]);
%
%   the normal contact forces at the contacts
    fn_pq = zeros(M,3);
    fn_pq = -Fn * n_pq;
%
%   directions of the unit tangential vecor in the direction of
%   tangential force
%   a = sqrt(1/6);
    a = 1 / sqrt(2 + 4*(Unit(1)/Unit(3))^2);
    b = 2*Unit(1)/Unit(3) * a;
    t1 = [-a, -a, -b]';
    t2 = [ a, -a, -b]';
    t3 = [ a,  a, -b]';
    t4 = [-a,  a, -b]';
%
    t_pq = zeros(M,3);
    t_pq( 1,:) = t1;
    t_pq( 3,:) = t2;
    t_pq( 5,:) = t3;
    t_pq( 7,:) = t4;
%
    t_pq( 9,:) = t3;
    t_pq(11,:) = t4;
    t_pq(13,:) = t1;
    t_pq(15,:) = t2;
%
    if Cells >=2
      t_pq(17,:) = t1;
      t_pq(19,:) = t2;
      t_pq(21,:) = t3;
      t_pq(23,:) = t4;
%
      t_pq(25,:) = t3;
      t_pq(27,:) = t4;
      t_pq(29,:) = t1;
      t_pq(31,:) = t2;
%
      if Cells ==3
        t_pq(33,:) = t1;
        t_pq(35,:) = t2;
        t_pq(37,:) = t3;
        t_pq(39,:) = t4;
%
        t_pq(41,:) = t3;
        t_pq(43,:) = t4;
        t_pq(45,:) = t1;
        t_pq(47,:) = t2;
      end
    end
%
    t_pq(2:2:end,:) = -t_pq(1:2:end,:);
%
%   the tangential contact forces at the contacts
    ft_pq = zeros(M,3);
    ft_pq = Ft * t_pq;
%
%   combined contact forces
    f_pq = fn_pq + ft_pq;
%
%   no contact moments at the contacts
    m_pq = zeros(size(f_pq));
%
%   lateral chamber pressure
    Pressure = 4 * f_pq(1,1) / (Unit(2) * Unit(3));
%
%------curvatures and curvature tensors of particles p and q ---------
%      at their contact pq
%
%   radii of curvature of particles p and q.  These radii will be stored in
%   the 2M x 2 matrix rho.  rho_pq(pq,1) is the radius of p; rho_pq(pq,2) is the
%   radius of q. Thus far, we have only identified contacts between disks.
%   The radius of curvature is the disk radius.
%   rho_pq = Radius*ones(M,2);
%
%   the radii of curvature of the particles at each contact. We will modify
%   the curvature of the base surface later
    if Curvature <= 0
      rho_pq = Radius * ones(M,2);
    else
      rho_pq = Curvature * Radius * ones(M,2);
    end
%
%   the contact vectors for each contact: from p to the contact pq, and
%   from q to the contact pq
    r_pq_p =  0.5 * rx;
    r_pq_q = -0.5 * rx;
%
%   we store the curvature tensors in a 2M x 2 x 3 x 3 array (a 4D array), K_pq.
%   Note that M is twice the number of contacts.  See Eqs. 27 and 32
    K_pq = zeros(M,2,3,3);
%
%   the pseud0-inverse of [k_p + K_q], a 2M x 3 x 3 array.  See Eq. 27.
    K_inv = zeros(M,3,3);
%
%   Each pq contact has a K_p curvature tensor and a K_q curvature tensor.
%   Each of these two tensors is 3 x 3.  Refer to Kuhn & Bagi (2004),
%   J. Engrg. Mech., vol. 130, no 7, 826-835.
    p = 1; q = 2;
    for pq = 1:M
      if idim==2
%       we will use this matrix in our calculations:
        Matrix = [ n_pq(pq,2)*n_pq(pq,2), -n_pq(pq,1)*n_pq(pq,2), 0; ...
                  -n_pq(pq,1)*n_pq(pq,2),  n_pq(pq,1)*n_pq(pq,1), 0; ...
                             0,                      0,           0];
%
%       the 3x3 K_p tensor for contact pq
        K_pq(pq,p,:,:) = -(1/rho_pq(pq,p)) * Matrix;
%
%       the 3x3 K_q tensor for contact pq
        K_pq(pq,q,:,:) = -(1/rho_pq(pq,q)) * Matrix;
%
%       also compute the 3x3 pseudo-inverse of [K_p + K_q], which involves
%       the quotient (rho_p * rho_q) / (rho_p + rho_q), as in eq. 41 of
%       Kuhn and Bagi (2004).
        if      isinf(rho_pq(pq,1)) && ~isinf(rho_pq(pq,2))
          Quotient = rho_pq(pq,2);
        elseif ~isinf(rho_pq(pq,1)) &&  isinf(rho_pq(pq,2))
          Quotient = rho_pq(pq,1);
        elseif ~isinf(rho_pq(pq,1)) && ~isinf(rho_pq(pq,2))
          Quotient = rho_pq(pq,1)*rho_pq(pq,2) / (rho_pq(pq,1) + rho_pq(pq,2));
        else
          Quotient = Inf;
        end
%
%       now, compute the 3x3 pseudo-inverse of [K_p + K_q]
        K_inv(pq,:,:) = - Quotient * Matrix;
%
      else
        if Curvature <= 0
          [K_pq, K_inv] = Curvature_Spheres( ...
                           pq,  n_pq(pq,:), Radius, Radius, K_pq, K_inv);
        else
          [K_pq, K_inv] = Curvature_Spheres( ...
                           pq,  n_pq(pq,:), ...
                           Curvature*Radius, Curvature*Radius, K_pq, K_inv);
        end
      end
    end
%
  end % if Condition==1
%
%-----contact stiffness and friction properties ---------
%
% contact frictional and stiffness properties. We allow for different
% properties for the disk-disk contacts
  mu_pq = mu * ones(M,1);
  k_pq = k * ones(M,1);
  alpha_pq = alpha * ones(M,1);
%
% the models for the contact force; and contact moments;
  FM_model_pq = 1*ones(M,1);
%
% "plates" that are used for computing the force produced on the particle
% by a chamber pressure on a membrane
%
% 1) for each particle, i
% 2) there are as many as 4 plates for exterior particles, and
% 3) two vectors (boundary vectors) are associated with each plate, and
% 4) each vector has three components
% 
  Plate = zeros(N,4,2,3);
  q_Plate = zeros(N,4,2);
%
% Base particles, Level 0
  Plate( 1,1,1,:) =  Unit(3)*[0; 0; 1]; q_Plate( 1,1,1) = 6;
  Plate( 1,1,2,:) = -Unit(2)*[0; 1; 0]; q_Plate( 1,1,2) = 4;
  Plate( 1,2,1,:) = -Unit(1)*[1; 0; 0]; q_Plate( 1,2,1) = 2;
  Plate( 1,2,2,:) =  Unit(3)*[0; 0; 1]; q_Plate( 1,2,2) = 6;
%
  Plate( 2,1,1,:) =  Unit(3)*[0; 0; 1]; q_Plate( 2,1,1) = 7;
  Plate( 2,1,2,:) =  Unit(1)*[1; 0; 0]; q_Plate( 2,1,2) = 1;
  Plate( 2,2,1,:) = -Unit(2)*[0; 1; 0]; q_Plate( 2,2,1) = 3;
  Plate( 2,2,2,:) =  Unit(3)*[0; 0; 1]; q_Plate( 2,2,2) = 7;
%
  Plate( 3,1,1,:) =  Unit(3)*[0; 0; 1]; q_Plate( 3,1,1) = 8;
  Plate( 3,1,2,:) =  Unit(2)*[0; 1; 0]; q_Plate( 3,1,2) = 2;
  Plate( 3,2,1,:) =  Unit(1)*[1; 0; 0]; q_Plate( 3,2,1) = 4;
  Plate( 3,2,2,:) =  Unit(3)*[0; 0; 1]; q_Plate( 3,2,2) = 8;
%
  Plate( 4,1,1,:) =  Unit(3)*[0; 0; 1]; q_Plate( 4,1,1) = 9;
  Plate( 4,1,2,:) = -Unit(1)*[1; 0; 0]; q_Plate( 4,1,2) = 3;
  Plate( 4,2,1,:) =  Unit(2)*[0; 1; 0]; q_Plate( 4,2,1) = 1;
  Plate( 4,2,2,:) =  Unit(3)*[0; 0; 1]; q_Plate( 4,2,2) = 9;
%
% Level 1 particles, top of first cell
  Plate( 6,1,1,:) = -Unit(2)*[0; 1; 0]; q_Plate( 6,1,1) = 9;
  Plate( 6,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate( 6,1,2) = 1;
  Plate( 6,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate( 6,2,1) = 1;
  Plate( 6,2,2,:) = -Unit(1)*[1; 0; 0]; q_Plate( 6,2,2) = 7;
%
  Plate( 7,1,1,:) =  Unit(1)*[1; 0; 0]; q_Plate( 7,1,1) = 6;
  Plate( 7,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate( 7,1,2) = 2;
  Plate( 7,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate( 7,2,1) = 2;
  Plate( 7,2,2,:) = -Unit(2)*[0; 1; 0]; q_Plate( 7,2,2) = 8;
% 
  Plate( 8,1,1,:) =  Unit(2)*[0; 1; 0]; q_Plate( 8,1,1) = 7;
  Plate( 8,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate( 8,1,2) = 3;
  Plate( 8,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate( 8,2,1) = 3;
  Plate( 8,2,2,:) =  Unit(1)*[1; 0; 0]; q_Plate( 8,2,2) = 9;
%     
  Plate( 9,1,1,:) = -Unit(1)*[1; 0; 0]; q_Plate( 9,1,1) = 8;
  Plate( 9,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate( 9,1,2) = 4;
  Plate( 9,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate( 9,2,1) = 4;
  Plate( 9,2,2,:) =  Unit(2)*[0; 1; 0]; q_Plate( 9,2,2) = 6;
%
  if Cells>=2
%   Level 1 particles, bottom of 2nd cell
    Plate( 6,3,1,:) =  Unit(3)*[0; 0; 1]; q_Plate( 6,3,1) = 11;
    Plate( 6,3,2,:) = -Unit(2)*[0; 1; 0]; q_Plate( 6,3,2) =  9;
    Plate( 6,4,1,:) = -Unit(1)*[1; 0; 0]; q_Plate( 6,4,1) =  7;
    Plate( 6,4,2,:) =  Unit(3)*[0; 0; 1]; q_Plate( 6,4,2) = 11;
%
    Plate( 7,3,1,:) =  Unit(3)*[0; 0; 1]; q_Plate( 7,3,1) = 12;
    Plate( 7,3,2,:) =  Unit(1)*[1; 0; 0]; q_Plate( 7,3,2) =  6;
    Plate( 7,4,1,:) = -Unit(2)*[0; 1; 0]; q_Plate( 7,4,1) =  8;
    Plate( 7,4,2,:) =  Unit(3)*[0; 0; 1]; q_Plate( 7,4,2) = 12;
%
    Plate( 8,3,1,:) =  Unit(3)*[0; 0; 1]; q_Plate( 8,3,1) = 13;
    Plate( 8,3,2,:) =  Unit(2)*[0; 1; 0]; q_Plate( 8,3,2) =  7;
    Plate( 8,4,1,:) =  Unit(1)*[1; 0; 0]; q_Plate( 8,4,1) =  9;
    Plate( 8,4,2,:) =  Unit(3)*[0; 0; 1]; q_Plate( 8,4,2) = 13;
%
    Plate( 9,3,1,:) =  Unit(3)*[0; 0; 1]; q_Plate( 9,3,1) = 14;
    Plate( 9,3,2,:) = -Unit(1)*[1; 0; 0]; q_Plate( 9,3,2) =  8;
    Plate( 9,4,1,:) =  Unit(2)*[0; 1; 0]; q_Plate( 9,4,1) =  6;
    Plate( 9,4,2,:) =  Unit(3)*[0; 0; 1]; q_Plate( 9,4,2) = 14;
%
%   Level 2 particles, top of 2nd cell
    Plate(11,1,1,:) = -Unit(2)*[0; 1; 0]; q_Plate(11,1,1) = 14;
    Plate(11,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate(11,1,2) =  6;
    Plate(11,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate(11,2,1) =  6;
    Plate(11,2,2,:) = -Unit(1)*[1; 0; 0]; q_Plate(11,2,2) = 12;
%
    Plate(12,1,1,:) =  Unit(1)*[1; 0; 0]; q_Plate(12,1,1) = 11;
    Plate(12,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate(12,1,2) =  7;
    Plate(12,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate(12,2,1) =  7;
    Plate(12,2,2,:) = -Unit(2)*[0; 1; 0]; q_Plate(12,2,2) = 13;
% 
    Plate(13,1,1,:) =  Unit(2)*[0; 1; 0]; q_Plate(13,1,1) = 12;
    Plate(13,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate(13,1,2) =  8;
    Plate(13,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate(13,2,1) =  8;
    Plate(13,2,2,:) =  Unit(1)*[1; 0; 0]; q_Plate(13,2,2) = 14;
%     
    Plate(14,1,1,:) = -Unit(1)*[1; 0; 0]; q_Plate(14,1,1) = 13;
    Plate(14,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate(14,1,2) =  9;
    Plate(14,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate(14,2,1) =  9;
    Plate(14,2,2,:) =  Unit(2)*[0; 1; 0]; q_Plate(14,2,2) = 11;
%
    if Cells>=3
%     Level 2 particles, bottom of 3rd cell
      Plate(11,3,1,:) =  Unit(3)*[0; 0; 1]; q_Plate(11,3,1) = 16;
      Plate(11,3,2,:) = -Unit(2)*[0; 1; 0]; q_Plate(11,3,2) = 14;
      Plate(11,4,1,:) = -Unit(1)*[1; 0; 0]; q_Plate(11,4,1) = 12;
      Plate(11,4,2,:) =  Unit(3)*[0; 0; 1]; q_Plate(11,4,2) = 16;
%
      Plate(12,3,1,:) =  Unit(3)*[0; 0; 1]; q_Plate(12,3,1) = 17;
      Plate(12,3,2,:) =  Unit(1)*[1; 0; 0]; q_Plate(12,3,2) = 11;
      Plate(12,4,1,:) = -Unit(2)*[0; 1; 0]; q_Plate(12,4,1) = 13;
      Plate(12,4,2,:) =  Unit(3)*[0; 0; 1]; q_Plate(12,4,2) = 17;
%
      Plate(13,3,1,:) =  Unit(3)*[0; 0; 1]; q_Plate(13,3,1) = 18;
      Plate(13,3,2,:) =  Unit(2)*[0; 1; 0]; q_Plate(13,3,2) = 12;
      Plate(13,4,1,:) =  Unit(1)*[1; 0; 0]; q_Plate(13,4,1) = 14;
      Plate(13,4,2,:) =  Unit(3)*[0; 0; 1]; q_Plate(13,4,2) = 18;
%
      Plate(14,3,1,:) =  Unit(3)*[0; 0; 1]; q_Plate(14,3,1) = 19;
      Plate(14,3,2,:) = -Unit(1)*[1; 0; 0]; q_Plate(14,3,2) = 13;
      Plate(14,4,1,:) =  Unit(2)*[0; 1; 0]; q_Plate(14,4,1) = 11;
      Plate(14,4,2,:) =  Unit(3)*[0; 0; 1]; q_Plate(14,4,2) = 19;
%
%     Level 3 particles, top of 3rd cell
      Plate(16,1,1,:) = -Unit(2)*[0; 1; 0]; q_Plate(16,1,1) = 19;
      Plate(16,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate(16,1,2) = 11;
      Plate(16,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate(16,2,1) = 11;
      Plate(16,2,2,:) = -Unit(1)*[1; 0; 0]; q_Plate(16,2,2) = 17;
%
      Plate(17,1,1,:) =  Unit(1)*[1; 0; 0]; q_Plate(17,1,1) = 16;
      Plate(17,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate(17,1,2) = 12;
      Plate(17,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate(17,2,1) = 12;
      Plate(17,2,2,:) = -Unit(2)*[0; 1; 0]; q_Plate(17,2,2) = 18;
% 
      Plate(18,1,1,:) =  Unit(2)*[0; 1; 0]; q_Plate(18,1,1) = 17;
      Plate(18,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate(18,1,2) = 13;
      Plate(18,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate(18,2,1) = 13;
      Plate(18,2,2,:) =  Unit(1)*[1; 0; 0]; q_Plate(18,2,2) = 19;
%     
      Plate(19,1,1,:) = -Unit(1)*[1; 0; 0]; q_Plate(19,1,1) = 18;
      Plate(19,1,2,:) = -Unit(3)*[0; 0; 1]; q_Plate(19,1,2) = 14;
      Plate(19,2,1,:) = -Unit(3)*[0; 0; 1]; q_Plate(19,2,1) = 14;
      Plate(19,2,2,:) =  Unit(2)*[0; 1; 0]; q_Plate(19,2,2) = 16;
    end
  end
