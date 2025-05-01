function [ istep, icont, ddef, ddefm, dstress, ...
           stress_target, Reset_target_stress, ...
           istop, increment, max_increments, nStress, iStress, ...
           lap_Time, start_Time] = ...
         Control_alt(istep, icontr, icont, defrat, igoal, finalv, idim, ...
                     I, J, Time);
%
  istep = istep + 1;
%
% whether to stop the simulation
  istop = 0;
%
% increment (time steps) within the new control period
  increment = 0;
%
% initialize output arguments
  max_increments = [];
% icont = [];
  ddef = [];
  ddefm = zeros(3,3);
  dstress = zeros(3,3);
  nStress = [];
  iStress = [];
  Reset_target_stress = [];
  stress_target = [];
  lap_Time = [];
  start_Time = [];
%
  if istep <= size(icontr,1)
%
%   continue the simulation
    istop = 0;
%
%   convert the icontr numeric value to a string
    Control_String = sprintf('%06i',icontr(istep));
%
%   whether the ith component is stress-controlled (=1) or 
%   strain-controlled (=0)
    icont = zeros(1,6);
    icont(1) = Control_String(1)=='1';
    icont(2) = Control_String(2)=='1';
    icont(3) = Control_String(3)=='1';
    icont(4) = Control_String(4)=='1';
    icont(5) = Control_String(5)=='1';
    icont(6) = Control_String(6)=='1';
%
    if idim==2
      if icont(3)~=0 || icont(5)~=0 || icont(6)~=0
        disp(cstrcat('WARNING: Your input values of icontr are ignored ', ...
                      'for your 2D simulation.'))
      end
      icont(3) = 0;
      icont(5) = 0;
      icont(6) = 0;
    end
%
%   number of controlled stresses
    nStress = sum(icont);
%
%   the stresses being controlled.  Integers in set [1:6]
    iStress = sort(find(icont==1));
%
%   the prescribed deformation increments
    ddef = zeros(1,6);
    ddef(find(icont==0)) = defrat(istep,find(icont==0));
%
%   the prescribed stress increments, as a 1x6 vector
    dstress2 = zeros(1,6);
    dstress2(find(icont==1)) = defrat(istep,find(icont==1));
    if idim==2
      dstress2([3, 5, 6]) = 0;
    end
%
%   stress increment as a 3x3 matrix
    dstress(1,1) = dstress2(1);
    dstress(2,2) = dstress2(2);
    dstress(3,3) = dstress2(3);
    dstress(1,2) = dstress2(4);
    dstress(1,3) = dstress2(5);
    dstress(2,3) = dstress2(6);
    dstress(2,1) = dstress(1,2);
    dstress(3,1) = dstress(1,3);
    dstress(3,2) = dstress(2,3);
%
%   we will need to reset the target stresses
    Reset_target_stress = 1;
%
%   change in the deformation gradient
    ddefm = zeros(3,3);
    ddefm(1,1) = ddef(1);
    ddefm(2,2) = ddef(2);
    ddefm(3,3) = ddef(3);
    ddefm(1,2) = ddef(4);
    ddefm(1,3) = ddef(5);
    ddefm(2,3) = ddef(6);
%
%   correct any invalid input for 2D simulations
    if idim==2
      icont([3, 5, 6]) = 0;
      ddef([3, 5, 6]) = 0;
      ddefm(3,3) = 0;
      ddefm(1,3) = 0;
      ddefm(2,3) = 0;
    end
%
%   establish termination criteria for this control step
    if igoal(istep) == 70
      max_increments = finalv(istep);
    elseif igoal(istep) == 80
      lap_Time = finalv(istep);
      start_Time = Time;
    end
% 
  else
%   stop the simulation
    istop = 1;
  end
