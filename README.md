The package LCPgranular was used by the author to solve the three 
example problems in the paper:
  Matthew R. Kuhn, "Quasi-static loading of granular media as a linear 
  complementarity problem"

The author ran the examples in a Linux environment.  The open source 
GNU Octave software was used to run the three examples.  After downloading 
the entire package, the examples can be run with Octave, within their 
respective directories, using the following scripts:

Example_1 : LCP2_Analyze_Three_Disks
Example_2 : Analyze_BCC
Example_3 : Driver_alt

Some examples require installation of the "optim" package from Octave Forge.
Data produced in Example_3 can be unbundled with the Octave script 
Plot_Results.m.

This package also includes the RunFile and StartFile for using the DEM
simulation software OVAL/DEMPLA.  This DEM software is in Fortra and is
available at 
  https://github.com/mrkuhn53/dempla
"dempla" must be compiled (e.g., gfortran) before running the simulation,
as explained in the DEMPLA manual.  The RunFile and StartFile simulate 
the 49-disk problem in Example 3 of the paper.  These files are contained
in the direcctor "Oval":
  RunFile =   Biaxial_Comp_k16_oval
  StartFile = Dcircls_49_a-57
The simulation is run with dempla in Mode 1.
%
The Figure 7(8?) of the paper was produced with the plotting software
gnuplot, and the figure uses LaTeX symbols.  The gnuplot figure is
created with the Linuxcommand: 
  gnu2pdf Disks_Results
