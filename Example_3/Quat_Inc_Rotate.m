function [Qp] = Quat_Inc_Rotate(Qp, dtheta)
%
% this function adjusts the particles' quaterions after the particles
% are rotated by the increments "dtheta"
%
%----------- INPUT ---------------------
% Qp = Nx4 array of particle quaternions
% dtheta = Nx3 increments of particle rotations
%
%----------- OUTPUT --------------------
% Qp = Nx4 array of new particle quaternions
%
% Compute the new quaternion of the particle frame:
%          Qp_new -> Qp + dQp = Qp + (1/2) dtheta*Qp, where "*" is
% the quaternion composition operator.  The updated
% quaternion Qpi_new will rotate the (Jager) xi-plane
% (whose normal is in the z-direction) into the (global)
% contact tangent plane
%
  wq21 = dtheta(:,1).*Qp(:,1);
  wq31 = dtheta(:,2).*Qp(:,1);
  wq41 = dtheta(:,3).*Qp(:,1);
  wq22 = dtheta(:,1).*Qp(:,2);
  wq32 = dtheta(:,2).*Qp(:,2);
  wq42 = dtheta(:,3).*Qp(:,2);
  wq23 = dtheta(:,1).*Qp(:,3);
  wq33 = dtheta(:,2).*Qp(:,3);
  wq43 = dtheta(:,3).*Qp(:,3);
  wq24 = dtheta(:,1).*Qp(:,4);
  wq34 = dtheta(:,2).*Qp(:,4);
  wq44 = dtheta(:,3).*Qp(:,4);
%
  Qp(:,1) = Qp(:,1) + 0.5*(-wq22 - wq33 - wq44);
  Qp(:,2) = Qp(:,2) + 0.5*( wq21 - wq43 + wq34);
  Qp(:,3) = Qp(:,3) + 0.5*( wq31 + wq42 - wq24);
  Qp(:,4) = Qp(:,4) + 0.5*( wq41 - wq32 + wq23);
%
% normalize the new quaternions
  Qp = ((1./sqrt(sum(Qp.^2,2)))*[1 1 1 1]) .* Qp;
