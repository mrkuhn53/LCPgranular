function [ K_pq, K_inv] = ...
         Curvature_Spheres(pq, n, Radius_p, Radius_q, K_pq, K_inv);
%
  Darboux1_p = zeros(1,3);
  Darboux2_p = zeros(1,3);
  Darboux3_p = zeros(1,3);
%
  Darboux1_q = zeros(1,3);
  Darboux2_q = zeros(1,3);
  Darboux3_q = zeros(1,3);
%
% Darboux unit vectors for the 1st (i) particle
%
  Darboux3_p = n;
%
  Darboux1_p(1) = -Darboux3_p(2);
  Darboux1_p(2) =  Darboux3_p(1);
%
% normalize
  Darboux1_p = (1/sqrt(sum(Darboux1_p.^2))) * Darboux1_p;
%
  Darboux2_p(1) =   Darboux3_p(2).*Darboux1_p(3) ...
                    - Darboux3_p(3).*Darboux1_p(2);
  Darboux2_p(2) =   Darboux3_p(3).*Darboux1_p(1) ...
                    - Darboux3_p(1).*Darboux1_p(3);
  Darboux2_p(3) =   Darboux3_p(1).*Darboux1_p(2) ...
                    - Darboux3_p(2).*Darboux1_p(1);
%
% Darboux unit vectors for the 2nd (j) particle
%
  Darboux3_q = -n;
%
  Darboux1_q(1) = -Darboux3_q(2);
  Darboux1_q(2) =  Darboux3_q(1);
%
% normalize
  Darboux1_q = (1/sqrt(sum(Darboux1_q.^2))) * Darboux1_q;
%
  Darboux2_q(1) =   Darboux3_q(2).*Darboux1_q(3) ...
                    - Darboux3_q(3).*Darboux1_q(2);
  Darboux2_q(2) =   Darboux3_q(3).*Darboux1_q(1) ...
                    - Darboux3_q(1).*Darboux1_q(3);
  Darboux2_q(3) =   Darboux3_q(1).*Darboux1_q(2) ...
                    - Darboux3_q(2).*Darboux1_q(1);
%
% Principal radii of curvature, 1st particle
  Rad1_p = Radius_p;
  Rad2_p = Radius_p;
%
% Principal radii of curvature, 2nd particle
  Rad1_q = Radius_q;
  Rad2_q = Radius_q;
%
  Kp_ap = zeros(3);
  Kq_aq = zeros(3);
%
  Kp_ap(1,1) = -1/Rad1_p;
  Kp_ap(2,2) = -1/Rad2_p;
%
  Kq_aq(1,1) = -1/Rad1_q;
  Kq_aq(2,2) = -1/Rad2_q;
%
  T_ap_e = [Darboux1_p;
            Darboux2_p;
            Darboux3_p];
%
  T_aq_e = [Darboux1_q;
            Darboux2_q;
            Darboux3_q];
%
  T_aq_ap = T_aq_e * T_ap_e';
%
  Kq_ap = T_aq_ap' * Kq_aq * T_aq_ap;
%
  Kp_Kq_ap = Kp_ap + Kq_ap;
%
% Kp_Kq_e = T_ap_e' * (Kp_Kq_ap) * T_ap_e;
%
  Kp_Kq_ap_inv = [inv(Kp_Kq_ap(1:2,1:2)) [0;0]; [0 0 0]];
  Kp_Kq_e_inv = T_ap_e' * Kp_Kq_ap_inv * T_ap_e;
%
  Kp_e = T_ap_e' * Kp_ap * T_ap_e;
  Kq_e = T_aq_e' * Kq_aq * T_aq_e;
%
  p = 1; q = 2;
  K_pq(pq,p,:,:) = Kp_e;
  K_pq(pq,q,:,:) = Kq_e;
  K_inv(pq,:,:) = Kp_Kq_e_inv;
