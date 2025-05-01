  function [Primary_version, Secondary_version, Tertiary_version, ...
            dt, dtmax, mass, Shape, np, e, sf, porosity, vcell, Model, ...
            knh, kn, kt, frict, G, Poisson, A1, ...
            LineNo, Line, Timer, istep, nupd, ipt, ...
            xloops, chi1, chi2, psi, sweep] = ...
    Read_Oval_screen_output(Oval_output_file_path)
%
% |+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|
%
% This function reads a file that contains the screen output of an
% OVAL or DEMPLA simulation
%
% Input:
%   Oval_output_file_path   the character string of the full RunFile name,
%                           including the full path to the RunFile
%
% Dependencies:  None
%
% Usage: when running OVAL or DEMPLA, redirect the screen output to a file,
% with these linux/unix commands:
%   dempla > screen_output_file
%
% where "dempla" is the path to the executable dempla file.  After running
% the dempla simulation, then call the function to retrieve the screen
% output:
%
% [Primary_version, Secondary_version, Tertiary_version, ...
%  dt, dtmax, mass, Shape, np, e, sf, porosity, vcell, Model, ...
%  knh, kn, kt, frict, G, Poisson, A1, ...
%  LineNo, Line, Timer, istep, nupd, ipt, ...
%  xloops, chi1, chi2, psi, sweep] = ...
% Read_Oval_screen_output(Oval_output_file_path);
%
% This function reads the screen output from an OVAL run.  This output
% can be stored to a file with redirection:
%   oval > output_file
%
  [Circle, Oval, Ellipse, Sphere, Ovoid, Nobby, Bumpy] = Shapes3();
%
  if ~exist(Oval_output_file_path, 'file')
    disp('File does not exist for Read_Oval_screen_output.m:')
    disp(Oval_output_file_path)
  end
%
  File_no = fopen(Oval_output_file_path, 'rt');
    fgetl(File_no); % Skip the first line
    Header =        fgetl(File_no); % Grab the second line
    oval_position = strfind(Header,'oval-');
    f_position =    strfind(Header,'.f');
    Version =       Header(oval_position+5:f_position-1);
    Dots =          strfind(Version, '.');
    Primary_version = str2num(Version(1:Dots(1)-1));
    Secondary_version = str2num(Version(Dots(1)+1:Dots(2)-1));
    Tertiary_version = str2num(Version(Dots(2)+1:end));
%
    if Primary_version==0 && Secondary_version==7 && Tertiary_version>=177
      for i = 3:19
        fgetl(File_no);  % Skip lines 3 - 19
      end
%
      dt_line = fgetl(File_no);
      dt = str2num(dt_line(strchr(dt_line,'0123456789')(1):end));
%
      dtmax_line = fgetl(File_no);
      dtmax = str2num(dtmax_line(strchr(dtmax_line,'0123456789')(1):end));
%
      mass_line = fgetl(File_no);
      mass = str2num(mass_line(strchr(mass_line,'0123456789')(1):end));
%
      fgetl(File_no);  % Skip a lines
%
      shape_line = fgetl(File_no);
%
      if ~isempty(strfind(shape_line, 'Circular'))
        Shape = Circle;
      elseif ~isempty(strfind(shape_line, 'Four'))
        Shape = Oval;
      elseif ~isempty(strfind(shape_line, 'Elliptical'))
        Shape = Ellipse;
      elseif ~isempty(strfind(shape_line, 'Spherical'))
        Shape = Sphere;
      elseif ~isempty(strfind(shape_line, 'Ovoid'))
        Shape = Ovoid;
      elseif ~isempty(strfind(shape_line, 'Nobby'))
        Shape = Nobby;
      elseif ~isempty(strfind(shape_line, 'Bumpy'))
        Shape = Bumpy;
      end
%
      fgetl(File_no);  % Skip a lines
%
      np_line = fgetl(File_no);
      np = str2num(np_line(strchr(np_line,'0123456789')(1):end));
%
      e_line = fgetl(File_no);
      e = str2num(e_line(strchr(e_line,'0123456789')(1):end));
%
      sf_line = fgetl(File_no);
      sf = str2num(sf_line(strchr(sf_line,'0123456789')(1):end));
%
      por_line = fgetl(File_no);
      porosity = str2num(por_line(strchr(por_line,'0123456789')(1):end));
%
      vcell_line = fgetl(File_no);
      vcell = str2num(vcell_line(strchr(vcell_line,'0123456789')(1):end));
%
      fgetl(File_no);  % Skip a lines
%
      model_line = fgetl(File_no);
%
      if ~isempty(strfind(model_line, 'Linear'))
        Model = 0;
      elseif ~isempty(strfind(model_line, 'Mindlin'))
        Model = 5;
      elseif ~isempty(strfind(model_line, 'conical'))
        Model = 7;
      elseif ~isempty(strfind(model_line, 'Kalker'))
        Model = 8;
      elseif ~isempty(strfind(model_line, 'general'))
        Model = 9;
      elseif ~isempty(strfind(model_line, 'Jager'))
        Model = 6;
      end
%
      knh = 0; kn = 0; kt = 0; frict = 0; G = 0; Poisson = 0; A1 = 0;
      LineNo = 0; Timer =  0; istep =  0; nupd =   0; ipt =    0;
      xloops = 0; chi1 =   0; chi2 =   0; psi =    0; sweep =  0;
%
      if Model==0
        knh_line = fgetl(File_no);
        knh = str2num(knh_line(strchr(knh_line,'0123456789')(1):end));
%
        kn_line = fgetl(File_no);
        kn = str2num(kn_line(strchr(kn_line,'0123456789')(1):end));
%     
        kt_line = fgetl(File_no);
        kt = str2num(kt_line(strchr(kt_line,'0123456789')(1):end));
        kt = kt*kn;
%     
        frict_line = fgetl(File_no);
        frict = str2num(frict_line(strchr(frict_line,'0123456789')(1):end));
      elseif Model>=5 && Model<=9
        G_line = fgetl(File_no);
        G = str2num(G_line(strchr(G_line,'0123456789')(1):end));
%     
        Poisson_line = fgetl(File_no);
        Poisson=str2num(Poisson_line(strchr(Poisson_line,'0123456789')(1):end));
%     
        frict_line = fgetl(File_no);
        frict = str2num(frict_line(strchr(frict_line,'0123456789')(1):end));
      end
%
      if Model==7 || Model==9
        A1_line = fgetl(File_no);
        A1 = str2num(A1_line(strchr(A1_line,'0123456789')(1):end));
      end
%
      if Model==9
        palpha_line = fgetl(File_no);
        palpha = str2num(A1_line(strchr(palpha_line,'0123456789')(1):end));
      end
%
      fgetl(File_no); % Skip a line
%
      M_line = fgetl(File_no);
      M = str2num(M_line(strchr(M_line,'0123456789')(1):end));
%
      ovrlap_line = fgetl(File_no);
      ovrlap = str2num(ovrlap_line(strchr(ovrlap_line,'0123456789')(1):end));
%
      fgetl(File_no); % Skip 3 lines
      fgetl(File_no);
      fgetl(File_no);
    end
%
    if 1
      if 1
        nFields = 9;
        TEMPLATE = '%*5c%g%g%g%g%g%g%g%g%g\n'; % ignore the first 5 columns
      elseif Model>=5 && Model<=9
        nFields = 10;
        TEMPLATE = '%g';
      else
%       Note: we place a maximum of 999 on istep
        nFields = 10;
        TEMPLATE = '%5i%11g%3i%6i%8i%6g%10g%10g%10g%6g\n';
      end
%
      [Data, Count] = fscanf(File_no, TEMPLATE);
%
      if 0
%       old style, doesn't work with ***** in first 5 columns
        LineNo = Data( 1:nFields:Count); % old style, doesn't work with *****
        Line = LineNo + 1;
        Timer =  Data( 2:nFields:Count);
        istep =  Data( 3:nFields:Count);
        nupd =   Data( 4:nFields:Count);
        ipt =    Data( 5:nFields:Count);
        xloops = Data( 6:nFields:Count);
        chi1 =   Data( 7:nFields:Count);
        chi2 =   Data( 8:nFields:Count);
        psi =    Data( 9:nFields:Count);
        sweep =  Data( 10:nFields:Count);
      elseif 1
        Line = [1:Count/nFields]';
        LineNo = Line - 1;
        Timer =  Data( 1:nFields:Count);
        istep =  Data( 2:nFields:Count);
        nupd =   Data( 3:nFields:Count);
        ipt =    Data( 4:nFields:Count);
        xloops = Data( 5:nFields:Count);
        chi1 =   Data( 6:nFields:Count);
        chi2 =   Data( 7:nFields:Count);
        psi =    Data( 8:nFields:Count);
        sweep =  Data( 9:nFields:Count);
      end
    elseif 0
%     A weird kluge so that we can read "istep" values between 1000 and 9999
      nFields = 24;
      TEMPLATE = '%g';
      TEMPLATE = '%5i%16c%6i%8i%6g%10g%10g%10g%7g';
      [Data, Count] = fscanf(File_no, TEMPLATE);
%
      LineNo = Data( 1:nFields:Count);
      Line = LineNo + 1;
%     Timer =  Data( 2:nFields:Count);
%     istep =  Data( 3:nFields:Count);
      nupd =   Data( 18:nFields:Count);
      ipt =    Data( 19:nFields:Count);
      xloops = Data( 20:nFields:Count);
      chi1 =   Data( 21:nFields:Count);
      chi2 =   Data( 22:nFields:Count);
      psi =    Data( 23:nFields:Count);
      sweep =  Data( 24:nFields:Count);
%
%     for Index = 1:length(sweep);
%       TimerString = Data((nFields*(Index-1)+2):(nFields*(Index-1)+13))';
%     end
%     char(TimerString)
%     length(sweep)
    end

  fclose(File_no);
