%
% |+|+|+|+|+|+|+|+|  Revision 2   |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
function [ P_L, P_Lperp, C_MP_inv, ...
           P_L_wo_constr, P_Lperp_wo_constr, C_MP_inv_wo_constr] = ...
         LCP3_Projections_Types_II_III(C, C_stress, C_wo_constr);
%
% -------------- INPUT ------------------------
% C =        r x 6N constraint matrix
%
% ----------- OUTPUT ------------------------------------------
% P_L =      6N x 6N matrix that projects onto the null space of C', KPD Eq. 58
% P_Lperp  = 6N x 6N matrix that projects onto the column space of C', KPD Eq. 59
% C_MP_inv = 6N x r Moore-Penrose inverse of C, near KPD Eqs. 67 & 68
%
%
% this function computes the matrices [P_L] and [P_Lperp], [Pnrr] in 
% KPD Eqs. 58-59.  These matrices are used with Types II and III constraints.
%
% append the stress constraints
  C = [C, C_stress];
%
% Moore-Penrose inverse near KPD Eq. 68
  C_MP_inv = C' * inv(C*C');
%
% the matrix P_Lperp in KPD Eq. 59
  P_Lperp = C_MP_inv * C;
%
% the matrix P_L in KPD Eq. 58
  P_L = eye(size(C,2)) - P_Lperp;
%
% Moore-Penrose inverse near KPD Eq. 68
  C_MP_inv_wo_constr = C_wo_constr' * inv(C_wo_constr*C_wo_constr');
%
% the matrix P_Lperp in KPD Eq. 59
  P_Lperp_wo_constr = C_MP_inv_wo_constr * C_wo_constr;
%
% the matrix P_L in KPD Eq. 58
  P_L_wo_constr = eye(size(C_wo_constr,2)) - P_Lperp_wo_constr;
