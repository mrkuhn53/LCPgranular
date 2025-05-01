function [algori, iinteg, ncownt, iout, ...
          istart, iend, idef, iupdtm, icirct, imodel, nplatn, nloop1, ...
          k, alpha, mu, frictw, rho, sep, pcrit, ...
          xseed, rmsvel, pdif, tmax, A_1, dt, ...
          ivers, iexact, isub, gravty, ifree, rfree, ...
          icontr, defrat, igoal, krotat, finalv, ipts, ...
          idump, iflexc, imicro, ibodyf, defdot, ipts2, iplot, ...
          Only_Fast_Method, Type_Constraint, Min_Coord_No, ...
          Max_Repetitions, M_lambda_max_enumerate, Sep, lambda_tolerance, ...
          Max_stress_deviation, Max_b_p, K_walls, lambda_threshold, ...
          if_Inst_Output, Choose_MaxEnt_Path] = ...
 Read_LCP_RunFile(RunFile_path) 
%
% |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
% This function reads an OVAL or DEMPLA RunFile and returns the various
% input variables and arrays in the file
%
% Input:
%   RunFile_path   the character string of the full RunFile name,
%                  including the full path to the RunFile
%
% Dependencies:  None
%
% Function call...
%  [algori, iinteg, ncownt, iout, ...
%   istart, iend, idef, iupdtm, icirct, imodel, nplatn, nloop1, ...
%   k, alpha, mu, frictw, rho, sep, pcrit, ...
%   xseed, rmsvel, pdif, tmax, A_1, dt, ...
%   ivers, iexact, isub, gravty, ifree, rfree, ...
%   icontr, defrat, igoal, krotat, finalv, ipts, ...
%   idump, iflexc, imicro, ibodyf, defdot, ipts2, iplot] = ...
%  Read_RunFile(RunFile_path);
%
  algori =   0;
  iinteg =   0;
  ncownt =   0;
  iout(1) =  0;
  iout(2) =  0;
  iout(3) =  0;
  istart =   0;
  iend =     0;
  idef =     0;
  iupdtm =   0;
  icirct =   0;
  imodel =   0;
  nplatn =   0;
  nloop1 =   0;
  k =        0;
  alpha =    0;
  mu =       0;
  frictw =   0;
  rho =      0;
  sep =      0;
  pcrit(1) = 0;
  pcrit(2) = 0;
  pcrit(3) = 0;
  xseed =    0;
  rmsvel =   0;
  pdif =     0;
  tmax =     0;
  A_1 =      0;
  dt =       0;
%
  iexact =   0;
  isub =     0;
  ifree =    zeros(8);
  gravty(1) = 0;
  gravty(2) = 0;
  gravty(3) = 0;
  ivers =     iinteg;
  rfree =     zeros(8);
%
  icontr = 0;
  defrat = zeros(1,6);
  igoal =  0;
  krotat = 0;
  finalv = 0;
  ipts =   0;
  idump =  0;
  iflexc = 0;
  imicro = 0;
  ibodyf = 0;
  defdot = 0;
  ipts2 =  0;
  iplot =  0;
%
  if ~exist(RunFile_path, 'file')
    disp('File does not exist for Read_RunFile.m:')
    disp(RunFile_path)
  end
%
  File_no = fopen(RunFile_path, 'rt');
  fgetl(File_no); % Skip the first line
%
  algori = fgetl(File_no); algori = algori(1:16); algori = str2num(algori);
  iinteg = fgetl(File_no); iinteg = iinteg(1:16); iinteg = str2num(iinteg);
  ivers = iinteg;
  ncownt = fgetl(File_no); ncownt = ncownt(1:16); ncownt = str2num(ncownt);
  iout2 =  fgetl(File_no); iout2 =  iout2(1:16);  iout2 =  str2num(iout2);
  iout(2) = iout2;
  iout3 =  fgetl(File_no); iout3 =  iout3(1:16);  iout3 =  str2num(iout3);
  iout(3) = iout3;
  istart = fgetl(File_no); istart = istart(1:16); istart = str2num(istart);
  iend =   fgetl(File_no); iend =   iend(1:16);   iend =   str2num(iend);
  idef =   fgetl(File_no); idef =   idef(1:16);   idef =   str2num(idef);
  iupdtm = fgetl(File_no); iupdtm = iupdtm(1:16); iupdtm = str2num(iupdtm);
  icirct = fgetl(File_no); icirct = icirct(1:16); icirct = str2num(icirct);
  imodel = fgetl(File_no); imodel = imodel(1:16); imodel = str2num(imodel);
  nplatn = fgetl(File_no); nplatn = nplatn(1:16); nplatn = str2num(nplatn);
  nloop1 = fgetl(File_no); nloop1 = nloop1(1:16); nloop1 = str2num(nloop1);
%
  if ivers==3 || ivers==4
    iexact = fgetl(File_no); iexact = iexact(1:16); iexact = str2num(iexact);
    isub =   fgetl(File_no); isub =   isub(1:16);   isub =   str2num(isub);
    ifree6 = fgetl(File_no); ifree6 = ifree6(1:16); ifree6 = str2num(ifree6);
    ifree7 = fgetl(File_no); ifree7 = ifree7(1:16); ifree7 = str2num(ifree7);
    ifree8 = fgetl(File_no); ifree8 = ifree8(1:16); ifree8 = str2num(ifree8);
  elseif ivers==0 || ivers==1
    iexact = 0;
    isub = 0;
    ifree6 = 0;
    ifree7 = 0;
    ifree8 = 0;
  end
%
  k = fgetl(File_no); k = k(1:16); k = str2num(k);
  alpha = fgetl(File_no); alpha = alpha(1:16); alpha = str2num(alpha);
  mu = fgetl(File_no); mu = mu(1:16); mu = str2num(mu);
  frictw = fgetl(File_no); frictw = frictw(1:16); frictw = str2num(frictw);
  rho = fgetl(File_no); rho = rho(1:16); rho = str2num(rho);
  sep = fgetl(File_no); sep = sep(1:16); sep = str2num(sep);
  pcrit1 = fgetl(File_no); pcrit1 = pcrit1(1:16); pcrit1 = str2num(pcrit1);
  pcrit2 = fgetl(File_no); pcrit2 = pcrit2(1:16); pcrit2 = str2num(pcrit2);
  pcrit3 = fgetl(File_no); pcrit3 = pcrit3(1:16); pcrit3 = str2num(pcrit3);
  pcrit(1) = pcrit1;
  pcrit(2) = pcrit2;
  pcrit(3) = pcrit3;
  xseed = fgetl(File_no); xseed = xseed(1:16); xseed = str2num(xseed);
  rmsvel = fgetl(File_no); rmsvel = rmsvel(1:16); rmsvel = str2num(rmsvel);
  pdif = fgetl(File_no); pdif = pdif(1:16); pdif = str2num(pdif);
  tmax = fgetl(File_no); tmax = tmax(1:16); tmax = str2num(tmax);
  A_1 = fgetl(File_no); A_1 = A_1(1:16); A_1 = str2num(A_1);
  dt = fgetl(File_no); dt = dt(1:16); dt = str2num(dt);
%
  if ivers==4
    gravty1 = fgetl(File_no); gravty1 = gravty1(1:16);gravty1=str2num(gravty1);
    gravty2 = fgetl(File_no); gravty2 = gravty2(1:16);gravty2=str2num(gravty2);
    gravty3 = fgetl(File_no); gravty3 = gravty3(1:16);gravty3=str2num(gravty3);
    gravty(1) = gravty1;
    gravty(2) = gravty2;
    gravty(3) = gravty3;
    rfree4 = fgetl(File_no); rfree4 = rfree4(1:16); rfree4 = str2num(rfree4);
    rfree5 = fgetl(File_no); rfree5 = rfree5(1:16); rfree5 = str2num(rfree5);
    rfree6 = fgetl(File_no); rfree6 = rfree6(1:16); rfree6 = str2num(rfree6);
    rfree7 = fgetl(File_no); rfree7 = rfree7(1:16); rfree7 = str2num(rfree7);
    rfree8 = fgetl(File_no); rfree8 = rfree8(1:16); rfree8 = str2num(rfree8);
    rfree(4) = rfree4;
    rfree(5) = rfree5;
    rfree(6) = rfree6;
    rfree(7) = rfree7;
    rfree(8) = rfree8;
  else
    gravty(1) = 0;
    gravty(2) = 0;
    gravty(3) = 0;
    rfree(4) = 0;
    rfree(5) = 0;
    rfree(6) = 0;
    rfree(7) = 0;
    rfree(8) = 0;
  end
%
% Skip 6 lines
  fgetl(File_no);
  fgetl(File_no);
  fgetl(File_no);
  fgetl(File_no);
  fgetl(File_no);
%
  nFields = 18;
  TEMPLATE = '%g';
  [Data, Count] = fscanf(File_no, TEMPLATE);
%
  icontr =   Data( 1:nFields:Count);
  defrat1 =  Data( 2:nFields:Count);
  defrat2 =  Data( 3:nFields:Count);
  defrat3 =  Data( 4:nFields:Count);
  defrat4 =  Data( 5:nFields:Count);
  defrat5 =  Data( 6:nFields:Count);
  defrat6 =  Data( 7:nFields:Count);
  igoal =    Data( 8:nFields:Count);
  krotat =   Data( 9:nFields:Count);
  finalv =   Data(10:nFields:Count);
  ipts =     Data(11:nFields:Count);
  idump =    Data(12:nFields:Count);
  iflexc =   Data(13:nFields:Count);
  imicro =   Data(14:nFields:Count);
  ibodyf =   Data(15:nFields:Count);
  defdot =   Data(16:nFields:Count);
  ipts2 =    Data(17:nFields:Count);
  iplot =    Data(18:nFields:Count);
%
  defrat = [defrat1, defrat2, defrat3, defrat4, defrat5, defrat6];
%
  fclose(File_no);
%
% alias the names of some input information
  Only_Fast_Method = ncownt;
  Type_Constraint = iout(2);
  Min_Coord_No = iout(3);
  Max_Repetitions = icirct;
  M_lambda_max_enumerate = nplatn;
  if_Inst_Output = nloop1;
  Choose_MaxEnt_Path = iexact;
  Sep = sep;
  lambda_tolerance = pcrit(1);
  Max_stress_deviation = pcrit(2);
  Max_b_p = pcrit(3);
  K_walls = xseed;
  lambda_threshold = rmsvel;
%
  if K_walls==0
    K_walls = zeros(6);
  end
%
  if Max_stress_deviation <=  0 || ~all(all(~K_walls))
    Max_stress_deviation = inf;
  end
%
  if Max_b_p <=  0
    Max_b_p = inf;
  end
%
  if lambda_threshold <= 0
    lambda_threshold = inf;
  end
