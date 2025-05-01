%function [Near_list_2, Min_Separation, Theta] = ...
function [ Near_list_2, ...
           u_update, xcell_update, nbins, bin_list, surround_bins, ...
           Threshold_separation, iupdat] = ...
         Near_Neighbors(N, Shape, Size, u, Qp, s_rad, b_rad, xlocal, ...
                        xcell, periodic_directions, rmax, ravg, njoin, ...
                        Sep, BinRatio)
% place the particles in bins
%
% basic particle shapes
  [Circle, Oval, Ellipse, Sphere, Ovoid, Nobby, Bumpy] = Shapes3();
%
% 2D or 3D assembly
  if Shape==Circle || Shape==Nobby
    idim = 2;
  elseif Shape==Sphere || Shape==Bumpy
    idim = 3;
  else
    disp('Only disks, spheres, nobbies, or bumpies allowed.')
    ERROR_Create_Assembly_D_1
  end
%
% reset values that will be used for determining when "Near_Neighbors" must
% be run again
  iupdat = 0;
  xcell_update = xcell;
  u_update = u;
%
% we will be finding the minimum separation (or maximum overlap) among
% the particles
  Min_Separation = inf;
%
% initialize
  Theta = [];
%
  if   (idim==2 && all(periodic_directions(1:2))) ...
    || (idim==3 && all(periodic_directions(1:3))) 
%   number of bins
    nbins = floor((1/2/rmax)*BinRatio*diag(xcell));
    if idim==2; nbins(3) = 0; end
%
    dims = 1:idim;
%
%   bin widths
%   bin_width = zeros(3,1);
%   bin_width(dims) = diag(xcell(dims,dims)) ./ nbins(dims);
%
%   the size of a bin, expressed as a 3x3 matrix, including shearing
    bin_cell = eye(3);
    for i = 1:idim
      bin_cell(:,i) = (1/nbins(i)) * xcell(:,i);
    end
%
%   assign particles to the bins in each direction
    bins = zeros(N,3);
    bins = floor(inv(bin_cell)*u')';
    for i = 1:idim
      bins(:,i) = mod(bins(:,i),nbins(i));
    end
%
%   assign bins to particles
%   bins = zeros(N,3);
%   bins(:,dims) = floor(u(:,dims) ./ (ones(N,1)*bin_width(dims)'));
  end
%
% initialize the cell array of near neighbors.  For example Near_list{10}
% is a list (vector) of particles that are near particle 10
  Near_list = cell(N,1);
%
  if idim==2
%   initialize the list of bins that are near to a given bin.
%   For example surround_bins{3,5} is a list of the 9 bins that are 
%   near to to bin 3,5:  the bins 3,4; 3,5; 3,6; 2,4; etc.
    surround_bins = cell(nbins(1),nbins(2));
%
    for bin1 = 1:nbins(1)
      for bin2 = 1:nbins(2)
%       create a list of the 9 surrounding bins, including the central bin
        surround_bins{bin1,bin2} = [bin1-1, bin2-1, 0;
                                    bin1-1, bin2  , 0;
                                    bin1-1, bin2+1, 0;
                                    bin1  , bin2-1, 0;
                                    bin1  , bin2  , 0;
                                    bin1  , bin2+1, 0;
                                    bin1+1, bin2-1, 0;
                                    bin1+1, bin2  , 0;
                                    bin1+1, bin2+1, 0];
%       in Octave/Matlab matrix indices must be 1 or higher.  To use modulo
%       remainders we must convert to indices that begin with 0
        surround_bins{bin1,bin2} = surround_bins{bin1,bin2} - 1;
%
        if   (idim==2 && all(periodic_directions(1:2))) ...
          || (idim==3 && all(periodic_directions(1:3))) 
%         use modulo remainders to account for periodic boundaries
          surround_bins{bin1,bin2}(:,1) = ...
            mod(surround_bins{bin1,bin2}(:,1),nbins(1));
          surround_bins{bin1,bin2}(:,2) = ...
            mod(surround_bins{bin1,bin2}(:,2),nbins(2));
          surround_bins{bin1,bin2}(:,3) = 0;
        end
      end
    end
%
%   lists of particles that are inside of each of the bins
    bin_list = cell(nbins(1),nbins(2));
%
%   lists of particles in the various bins
    Ones = ones(N,1);
    for i = 0:nbins(1)-1
      for j = 0:nbins(2)-1
%       find all particles that are in bin i,j,0
        bin_list{i+1,j+1} = find(all(bins==Ones*[i,j,0],2));
      end
    end
%
%   for each i particle, search through the 9 nearby bins, and add the
%   particles in those bins to the list of particles near to particle i
    for p = 1:N
      for j = 1:9
        Near_list{p} = [Near_list{p}, ...
                      bin_list{surround_bins{bins(p,1)+1,bins(p,2)+1}(j,1)+1, ...
                               surround_bins{bins(p,1)+1,bins(p,2)+1}(j,2)+1}'];
      end
    end
%
  elseif idim==3
%   initialize the list of bins that are near to a given bin.
%   For example surround_bins{3,5} is a list of the 27 bins that are 
%   near to to bin 3,5:  the bins 3,4; 3,5; 3,6; 2,4; etc.
    surround_bins = cell(nbins(1),nbins(2),nbins(3));
%
    for bin1 = 1:nbins(1)
      for bin2 = 1:nbins(2)
        for bin3 = 1:nbins(3)
%         create a list of the 27 surrounding bins, including the central bin
          surround_bins{bin1,bin2,bin3} = [bin1-1, bin2-1, bin3-1;
                                           bin1-1, bin2  , bin3-1;
                                           bin1-1, bin2+1, bin3-1;
                                           bin1  , bin2-1, bin3-1;
                                           bin1  , bin2  , bin3-1;
                                           bin1  , bin2+1, bin3-1;
                                           bin1+1, bin2-1, bin3-1;
                                           bin1+1, bin2  , bin3-1;
                                           bin1+1, bin2+1, bin3-1;
                                           bin1-1, bin2-1, bin3  ;
                                           bin1-1, bin2  , bin3  ;
                                           bin1-1, bin2+1, bin3  ;
                                           bin1  , bin2-1, bin3  ;
                                           bin1  , bin2  , bin3  ;
                                           bin1  , bin2+1, bin3  ;
                                           bin1+1, bin2-1, bin3  ;
                                           bin1+1, bin2  , bin3  ;
                                           bin1+1, bin2+1, bin3  ;
                                           bin1-1, bin2-1, bin3+1;
                                           bin1-1, bin2  , bin3+1;
                                           bin1-1, bin2+1, bin3+1;
                                           bin1  , bin2-1, bin3+1;
                                           bin1  , bin2  , bin3+1;
                                           bin1  , bin2+1, bin3+1;
                                           bin1+1, bin2-1, bin3+1;
                                           bin1+1, bin2  , bin3+1;
                                           bin1+1, bin2+1, bin3+1];
%         in Octave/Matlab matrix indices must be 1 or higher.  To use modulo
%         remainders we must convert to indices that begin with 0
          surround_bins{bin1,bin2,bin3} = surround_bins{bin1,bin2,bin3} - 1;
%
          if   (idim==2 && all(periodic_directions(1:2))) ...
            || (idim==3 && all(periodic_directions(1:3))) 
%           use modulo remainders to account for periodic boundaries
            surround_bins{bin1,bin2,bin3}(:,1) = ...
              mod(surround_bins{bin1,bin2,bin3}(:,1),nbins(1));
            surround_bins{bin1,bin2,bin3}(:,2) = ...
              mod(surround_bins{bin1,bin2,bin3}(:,2),nbins(2));
            surround_bins{bin1,bin2,bin3}(:,3) = ...
              mod(surround_bins{bin1,bin2,bin3}(:,3),nbins(3));
          end
        end
      end
    end
%
%   lists of particles that are inside of each of the bins
    bin_list = cell(nbins(1),nbins(2),nbins(3));
%
%   lists of particles in the various bins
    Ones = ones(N,1);
    for i = 0:nbins(1)-1
      for j = 0:nbins(2)-1
        for k = 0:nbins(3)-1
%         find all particles that are in bin i,j,k
          bin_list{i+1,j+1,k+1} = find(all(bins==Ones*[i,j,k],2));
        end
      end
    end
%
%   for each i particle, search through the 27 nearby bins, and add the
%   particles in those bins to the list of particles near to particle i
    for p = 1:N
      sub_list = cell(27,1);
      for j = 1:27
        sub_list{j} = ...
          bin_list{ ...
             surround_bins{bins(p,1)+1,bins(p,2)+1,bins(p,3)+1}(j,1)+1, ...
             surround_bins{bins(p,1)+1,bins(p,2)+1,bins(p,3)+1}(j,2)+1, ...
             surround_bins{bins(p,1)+1,bins(p,2)+1,bins(p,3)+1}(j,3)+1}';
      end
      Near_list{p} = horzcat(sub_list{:});
    end
  end
%
% create an array M_near x 2 of particles p and q near-neighbor pairs
%
% find the lengths of the Near_list for each particle
  Near_List_Lengths = cellfun(@length, Near_list);
%
% the bottoms and tops of the list of all near-neighors
  Bots = [1; cumsum(Near_List_Lengths)(1:end-1)+1];
  Tops = cumsum(Near_List_Lengths);
%
% create the array M_near x 2 of particles p and q of near-neighbors
  Near_list_2 = zeros(sum(Near_List_Lengths),2);
  for p = 1:N
    Near_list_2(Bots(p):Tops(p),1) = p;
    Near_list_2(Bots(p):Tops(p),2) = Near_list{p};;
  end
%
% remove the qp variants, leaving only a pq list, with p < q
  pq_rows = find(Near_list_2(:,1) < Near_list_2(:,2));
  Near_list_2 = Near_list_2(pq_rows,:);
%
% now, we have a list of candidates for near-neighbors.  We will cull
% the list, so that the list only includes particles that are within a
% threshold distance of separation
%
% branch vectors between the centers of the two particles
  rxcntr = u(Near_list_2(:,2),:) - u(Near_list_2(:,1),:);
%
  if   (idim==2 && all(periodic_directions(1:2))) ...
    || (idim==3 && all(periodic_directions(1:3))) 
%   the number of periodic cell boundaries between the particle pair
%   (for example, when the pair straddles a periodic boundary)
    rcells = (inv(xcell) * rxcntr')';
%
%   the integer number of periodic boundaries crossed by the particle pair.
    nintx = round(rcells);
%
%   adjust the distance between the two particles, adjusting for the 
%   number of periodic boundaries crossed by the particle pair. The
%   branch vector between the centers of the two particles
    rxcntr = rxcntr - (xcell * nintx')';
  end
%
% threshold separation
  Threshold_separation = Sep * ravg;
%
  if Shape==Circle || Shape==Sphere
%   separation between particles in the near-neighbor list
    Separation = sqrt(sum(rxcntr.^2,2)) ...
                 - (Size(Near_list_2(:,1)) + Size(Near_list_2(:,2)));
%
    Min_Separation = min(Separation);
%
    Near_list_2 = Near_list_2(find(Separation < Threshold_separation),:);
%
%   we append 1's to represent the 1st components of the circle/sphere shapes
    Near_list_2 = [Near_list_2, ones(size(Near_list_2))];
%
  elseif Shape==Nobby
%
%   create a njoin^2 x 2 array of combinations of the nobs of p and q
    Combinations = zeros(njoin^2,2);
    for i = 1:njoin
      Combinations(((i-1)*njoin+1):i*njoin,1) = i;
      Combinations(((i-1)*njoin+1):i*njoin,2) = 1:njoin;
    end
    nCombinations = njoin^2;
%
%   convert quaternions to orientation angles
    Theta = 2*atan2(Qp(:,4), Qp(:,1));
    cosTheta = cos(Theta);
    sinTheta = sin(Theta);
%
%   create a shadow near-neighbor list to collect the culled neighbors
    Near_list_3 = zeros(size(Near_list_2,1)*njoin^2,4);
%
%   total number of near-neighbor combinations
    tCombination = 0;
%
    for pq = 1:size(Near_list_2,1)
%
%     the p and q particles
      p = Near_list_2(pq,1); q = Near_list_2(pq,2);
%
%     the vectors for the center of the p-nob to the center of the q-nob
      rxcntr_2 = zeros(njoin^2,3);
%
%     equal to rxcntr vector from the center of p to the center of q
%     plus vector from center of q to the q-nob minus vector from center
%     of p to p-nob
      rxcntr_2 = ones(njoin^2,1)*rxcntr(pq,:) ...
                 + ...
                 (Size(q)*b_rad(Combinations(:,2)))*[1 1 1] ...
                 .* ...
                 [  cosTheta(q)*xlocal(1,Combinations(:,2))' ...
                  - sinTheta(q)*xlocal(2,Combinations(:,2))', ...
                    sinTheta(q)*xlocal(1,Combinations(:,2))' ...
                  + cosTheta(q)*xlocal(2,Combinations(:,2))', ...
                  zeros(njoin^2,1)] ...
                 - ...
                 (Size(p)*b_rad(Combinations(:,1)))*[1 1 1] ...
                 .* ...
                 [  cosTheta(p)*xlocal(1,Combinations(:,1))' ...
                  - sinTheta(p)*xlocal(2,Combinations(:,1))', ...
                    sinTheta(p)*xlocal(1,Combinations(:,1))' ...
                  + cosTheta(p)*xlocal(2,Combinations(:,1))', ...
                  zeros(njoin^2,1)];
%
%     separation between particles in the near-neighbor l1st
      Separation = sqrt(sum(rxcntr_2.^2,2)) ...
                   - (  Size(p)*s_rad(Combinations(:,1)) ...
                      + Size(q)*s_rad(Combinations(:,2)));
%
%     if the pair is a near-neighbor, append to the temporary list of
%     near-neighbors
      List = Separation < Threshold_separation;
%
      Min_Separation = min(Min_Separation, min(Separation));
%
%     cull the list
      nList = sum(List);
      List = find(List);
%
%     append the near neighbors
      Near_list_3(tCombination+1:tCombination+nList,:) = ...
        [p*ones(size(List)), q*ones(size(List)), ...
         Combinations(List,1), Combinations(List,2)];
      tCombination = tCombination + nList;
%
    end
%
    Near_list_2 = Near_list_3(1:tCombination,:);
%
  elseif Shape==Bumpy
%
%   rotation matrices of N particles in a 3x3xN array
    [QRot] = RotMatrix_from_Quat(Qp);
%
%   create a njoin^2 x 2 array of combinations of the bumps of p and q
    Combinations = zeros(njoin^2,2);
    for i = 1:njoin
      Combinations(((i-1)*njoin+1):i*njoin,1) = i;
      Combinations(((i-1)*njoin+1):i*njoin,2) = 1:njoin;
    end
    nCombinations = njoin^2;
%
%   create a shadow near-neighbor list to collect the culled neighbors
    Near_list_3 = zeros(size(Near_list_2,1)*nCombinations,4);
%
%   total number of near-neighbor combinations
    tCombination = 0;
%
    for pq = 1:size(Near_list_2,1)
%
%     the p and q particles
      p = Near_list_2(pq,1); q = Near_list_2(pq,2);
%
%     the vectors for the center of the p-nob to the center of the q-nob
      rxcntr_2 = zeros(njoin^2,3);
%
%     equal to rxcntr vector from the center of p to the center of q
%     plus vector from center of q to the q-nob minus vector from center
%     of p to p-nob
      rxcntr_2 = ones(njoin^2,1)*rxcntr(pq,:) ...
                 + ...
                 (Size(q)*b_rad(Combinations(:,2)))*[1 1 1] ...
                 .* ...
                 (QRot(:,:,q)' * xlocal(:,Combinations(:,2)))' ...
                 - ...
                 (Size(p)*b_rad(Combinations(:,1)))*[1 1 1] ...
                 .* ...
                 (QRot(:,:,p)' * xlocal(:,Combinations(:,1)))';
%
%     separation between particles in the near-neighbor l1st
      Separation = sqrt(sum(rxcntr_2.^2,2)) ...
                   - (  Size(p)*s_rad(Combinations(:,1)) ...
                      + Size(q)*s_rad(Combinations(:,2)));
%
%     if the pair is a near-neighbor, append to the temporary list of
%     near-neighbors
      List = Separation < Threshold_separation;
%
      Min_Separation = min(Min_Separation, min(Separation));
%
%     cull the list
      nList = sum(List);
      List = find(List);
%
%     append the near neighbors
      Near_list_3(tCombination+1:tCombination+nList,:) = ...
        [p*ones(size(List)), q*ones(size(List)), ...
         Combinations(List,1), Combinations(List,2)];
      tCombination = tCombination + nList;
%
    end
%
    Near_list_2 = Near_list_3(1:tCombination,:);
  end
