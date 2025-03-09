import pytest
import numpy as np
import os
from io import StringIO
import tempfile
import ase
from VisualizePhonon.vibrational_analysis import VibrationalMode, VibrationAnalysis
from VisualizePhonon.vibrational_analysis_io import read_file, save_xsf


class TestVibrationalMode:
    """VibrationalModeクラスのテスト"""

    def test_init_real_frequency(self):
        """実数周波数の初期化テスト"""
        atoms = ase.Atoms(["H", "H", "O"])
        freq = 2.5
        eigenvector = np.ones((3, 3))
        mode = VibrationalMode(atoms, freq, eigenvector)

        assert mode.frequency == freq
        assert np.array_equal(mode.eigenvector, eigenvector)
        assert mode.frequency == 2.5

    def test_init_imaginary_frequency(self):
        """負周波数の初期化テスト"""
        atoms = ase.Atoms(["H", "H", "O"])
        freq = -3.6
        eigenvector = np.ones((3, 3))
        mode = VibrationalMode(atoms, freq, eigenvector)

        assert mode.frequency == freq
        assert np.array_equal(mode.eigenvector, eigenvector)
        assert mode.frequency == -3.6

    def test_repr(self):
        """__repr__メソッドのテスト"""
        # 実数周波数
        atoms = ase.Atoms(["H"])
        mode_real = VibrationalMode(atoms, 2.5, np.ones((1, 3)))
        assert "VibrationalMode(frequency=2.500000 THz, real)" in repr(
            mode_real)

        # 虚数周波数
        atoms = ase.Atoms(["H"])
        mode_imag = VibrationalMode(atoms, -3.6, np.ones((1, 3)))
        assert "VibrationalMode(frequency=-3.600000 THz, imaginary)" in repr(mode_imag)


class TestVibrationAnalysis:
    """VibrationAnalysisクラスのテスト"""

    @pytest.fixture
    def mock_outcar(self):
        """モックOUTCARファイルを作成する"""
        f = open('OUTCAR', 'r')
        data = f.read()
        f.close()

        outcar_content = """
vasp.5.4.4.18Apr17-6-g9f103f2a35 (build May 29 2022 01:12:50) complex          
  
 executed on             LinuxIFC date 2025.03.09  00:51:20
 running on  128 total cores
 distrk:  each k-point on  128 cores,    1 groups
 distr:  one band on NCORES_PER_BAND=   1 cores,  128 groups


--------------------------------------------------------------------------------------------------------


 INCAR:
 POTCAR:   PAW_RPBE O 08Apr2002                   
 POTCAR:   PAW_RPBE H 15Jun2001                   

 ----------------------------------------------------------------------------- 
|                                                                             |
|           W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!           |
|           W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!           |
|           W    W  A    A  R    R  N N  N  II  N N  N  G       !!!           |
|           W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !            |
|           WW  WW  A    A  R   R   N   NN  II  N   NN  G    G                |
|           W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!           |
|                                                                             |
|      For optimal performance we recommend to set                            |
|        NCORE= 4 - approx SQRT( number of cores)                             |
|      NCORE specifies how many cores store one orbital (NPAR=cpu/NCORE).     |
|      This setting can  greatly improve the performance of VASP for DFT.     |
|      The default,   NCORE=1            might be grossly inefficient         |
|      on modern multi-core architectures or massively parallel machines.     |
|      Do your own testing !!!!                                               |
|      Unfortunately you need to use the default for GW and RPA calculations. |
|      (for HF NCORE is supported but not extensively tested yet)             |
|                                                                             |
 ----------------------------------------------------------------------------- 

 POTCAR:   PAW_RPBE O 08Apr2002                   
   VRHFIN =O: s2p4                                                              
   LEXCH  = PE                                                                  
   EATOM  =   432.3788 eV,   31.7789 Ry                                         
                                                                                
   TITEL  = PAW_RPBE O 08Apr2002                                                
   LULTRA =        F    use ultrasoft PP ?                                      
   IUNSCR =        0    unscreen: 0-lin 1-nonlin 2-no                           
   RPACOR =     .000    partial core radius                                     
   POMASS =   16.000; ZVAL   =    6.000    mass and valenz                      
   RCORE  =    1.520    outmost cutoff radius                                   
   RWIGS  =    1.550; RWIGS  =     .820    wigner-seitz radius (au A)           
   ENMAX  =  400.000; ENMIN  =  300.000 eV                                      
   ICORE  =        2    local potential                                         
   LCOR   =        T    correct aug charges                                     
   LPAW   =        T    paw PP                                                  
   EAUG   =  605.392                                                            
   DEXC   =     .000                                                            
   RMAX   =    2.264    core radius for proj-oper                               
   RAUG   =    1.300    factor for augmentation sphere                          
   RDEP   =    1.550    radius for radial grids                                 
   QCUT   =   -5.520; QGAM   =   11.041    optimization parameters              
                                                                                
   Description                                                                  
     l     E      TYP  RCUT    TYP  RCUT                                        
     0   .000     23  1.200                                                     
     0  -.700     23  1.200                                                     
     1   .000     23  1.520                                                     
     1   .600     23  1.520                                                     
     2   .000      7  1.500                                                     
  local pseudopotential read in
  atomic valenz-charges read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
    PAW grid and wavefunctions read in
 
   number of l-projection  operators is LMAX  =           4
   number of lm-projection operators is LMMAX =           8
 
 POTCAR:   PAW_RPBE H 15Jun2001                   
   VRHFIN =H: ultrasoft test                                                    
   LEXCH  = PE                                                                  
   EATOM  =    12.4884 eV,     .9179 Ry                                         
                                                                                
   TITEL  = PAW_RPBE H 15Jun2001                                                
   LULTRA =        F    use ultrasoft PP ?                                      
   IUNSCR =        0    unscreen: 0-lin 1-nonlin 2-no                           
   RPACOR =     .000    partial core radius                                     
   POMASS =    1.000; ZVAL   =    1.000    mass and valenz                      
   RCORE  =    1.100    outmost cutoff radius                                   
   RWIGS  =     .700; RWIGS  =     .370    wigner-seitz radius (au A)           
   ENMAX  =  250.000; ENMIN  =  200.000 eV                                      
   RCLOC  =     .701    cutoff for local pot                                    
   LCOR   =        T    correct aug charges                                     
   LPAW   =        T    paw PP                                                  
   EAUG   =  400.000                                                            
   RMAX   =    2.174    core radius for proj-oper                               
   RAUG   =    1.200    factor for augmentation sphere                          
   RDEP   =    1.112    radius for radial grids                                 
   QCUT   =   -5.749; QGAM   =   11.498    optimization parameters              
                                                                                
   Description                                                                  
     l     E      TYP  RCUT    TYP  RCUT                                        
     0   .000     23  1.100                                                     
     0   .500     23  1.100                                                     
     1  -.300     23  1.100                                                     
  local pseudopotential read in
  atomic valenz-charges read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
    PAW grid and wavefunctions read in
 
   number of l-projection  operators is LMAX  =           3
   number of lm-projection operators is LMMAX =           5
 
 PAW_RPBE O 08Apr2002                   :
 energy of atom  1       EATOM= -432.3788
 kinetic energy error for atom=    0.1156 (will be added to EATOM!!)
 PAW_RPBE H 15Jun2001                   :
 energy of atom  2       EATOM=  -12.4884
 kinetic energy error for atom=    0.0098 (will be added to EATOM!!)
 
 
 POSCAR: H2O _2                                  
  positions in cartesian coordinates
  No initial velocities read in
 exchange correlation table for  LEXCH =        8
   RHO(1)=    0.500       N(1)  =     2000
   RHO(2)=  100.500       N(2)  =     4000
 


--------------------------------------------------------------------------------------------------------


 ion  position               nearest neighbor table
   1  0.000  0.000  0.000-   2 0.97   3 0.97
   2  0.075  0.904  0.000-   1 0.97
   3  0.075  0.096  0.000-   1 0.97
 
  LATTYP: Found a simple cubic cell.
 ALAT       =     8.0000000000
  
  Lattice vectors:
  
 A1 = (   8.0000000000,   0.0000000000,   0.0000000000)
 A2 = (   0.0000000000,   8.0000000000,   0.0000000000)
 A3 = (   0.0000000000,   0.0000000000,   8.0000000000)


Analysis of symmetry for initial positions (statically):
=====================================================================
 Subroutine PRICEL returns:
 Original cell was already a primitive cell.
 

 Routine SETGRP: Setting up the symmetry group for a 
 simple cubic supercell.


 Subroutine GETGRP returns: Found  4 space group operations
 (whereof  4 operations were pure point group operations)
 out of a pool of 48 trial point group operations.


The static configuration has the point symmetry C_2v.


Analysis of symmetry for dynamics (positions and initial velocities):
=====================================================================
 Subroutine PRICEL returns:
 Original cell was already a primitive cell.
 

 Routine SETGRP: Setting up the symmetry group for a 
 simple cubic supercell.


 Subroutine GETGRP returns: Found  4 space group operations
 (whereof  4 operations were pure point group operations)
 out of a pool of 48 trial point group operations.


The dynamic configuration has the point symmetry C_2v.


 Subroutine INISYM returns: Found  4 space group operations
 (whereof  4 operations are pure point group operations),
 and found     1 'primitive' translations

 
 
 KPOINTS: Gamma-point only                        

Automatic generation of k-mesh.
Space group operators:
 irot       det(A)        alpha          n_x          n_y          n_z        tau_x        tau_y        tau_z
    1     1.000000     0.000000     1.000000     0.000000     0.000000     0.000000     0.000000     0.000000
    2    -1.000000   180.000000     0.000000     0.000000     1.000000     0.000000     0.000000     0.000000
    3     1.000000   180.000000     1.000000     0.000000     0.000000     0.000000     0.000000     0.000000
    4    -1.000000   180.000000     0.000000     1.000000     0.000000     0.000000     0.000000     0.000000
 
 Subroutine IBZKPT returns following result:
 ===========================================
 
 Found      1 irreducible k-points:
 
 Following reciprocal coordinates:
            Coordinates               Weight
  0.000000  0.000000  0.000000      1.000000
 
 Following cartesian coordinates:
            Coordinates               Weight
  0.000000  0.000000  0.000000      1.000000
 
 
 Subroutine IBZKPT_HF returns following result:
 ==============================================
 
 Found      1 k-points in 1st BZ
 the following      1 k-points will be used (e.g. in the exchange kernel)
 Following reciprocal coordinates:   # in IRBZ
  0.000000  0.000000  0.000000    1.00000000   1 t-inv F


--------------------------------------------------------------------------------------------------------




 Dimension of arrays:
   k-points           NKPTS =      1   k-points in BZ     NKDIM =      1   number of bands    NBANDS=    128
   number of dos      NEDOS =    301   number of ions     NIONS =      3
   non local maximal  LDIM  =      4   non local SUM 2l+1 LMDIM =      8
   total plane-waves  NPLWV = 157464
   max r-space proj   IRMAX =      1   max aug-charges    IRDMAX=   4189
   dimension x,y,z NGX =    54 NGY =   54 NGZ =   54
   dimension x,y,z NGXF=   108 NGYF=  108 NGZF=  108
   support grid    NGXF=   108 NGYF=  108 NGZF=  108
   ions per type =               1   2
   NGX,Y,Z   is equivalent  to a cutoff of  11.22, 11.22, 11.22 a.u.
   NGXF,Y,Z  is equivalent  to a cutoff of  22.44, 22.44, 22.44 a.u.

 SYSTEM =  H2O vibration                           
 POSCAR =  H2O _2                                  

 Startparameter for this run:
   NWRITE =      2    write-flag & timer
   PREC   = a         normal or accurate (medium, high low for compatibility)
   ISTART =      0    job   : 0-new  1-cont  2-samecut
   ICHARG =      2    charge: 1-file 2-atom 10-const
   ISPIN  =      1    spin polarized calculation?
   LNONCOLLINEAR =      F non collinear calculations
   LSORBIT =      F    spin-orbit coupling
   INIWAV =      1    electr: 0-lowe 1-rand  2-diag
   LASPH  =      F    aspherical Exc in radial PAW
   METAGGA=      F    non-selfconsistent MetaGGA calc.

 Electronic Relaxation 1
   ENCUT  =  400.0 eV  29.40 Ry    5.42 a.u.  13.05 13.05 13.05*2*pi/ulx,y,z
   ENINI  =  400.0     initial cutoff
   ENAUG  =  605.4 eV  augmentation charge cutoff
   NELM   =     60;   NELMIN=  2; NELMDL= -5     # of ELM steps 
   EDIFF  = 0.1E-07   stopping-criterion for ELM
   LREAL  =      F    real-space projection
   NLSPLINE    = F    spline interpolate recip. space projectors
   LCOMPAT=      F    compatible to vasp.4.4
   GGA_COMPAT  = T    GGA compatible to vasp.4.4-vasp.4.6
   LMAXPAW     = -100 max onsite density
   LMAXMIX     =    2 max onsite mixed and CHGCAR
   VOSKOWN=      0    Vosko Wilk Nusair interpolation
   ROPT   =    0.00000   0.00000
 Ionic relaxation
   EDIFFG = 0.1E-06   stopping-criterion for IOM
   NSW    =     73    number of steps for IOM
   NBLOCK =      1;   KBLOCK =     73    inner block; outer block 
   IBRION =      6    ionic relax: 0-MD 1-quasi-New 2-CG
   NFREE  =      2    steps in history (QN), initial steepest desc. (CG)
   ISIF   =      2    stress and relaxation
   IWAVPR =     11    prediction:  0-non 1-charg 2-wave 3-comb
   ISYM   =      2    0-nonsym 1-usesym 2-fastsym
   LCORR  =      T    Harris-Foulkes like correction to forces

   POTIM  = 0.0150    time-step for ionic-motion
   TEIN   =    0.0    initial temperature
   TEBEG  =    0.0;   TEEND  =   0.0 temperature during run
   SMASS  =  -3.00    Nose mass-parameter (am)
   estimated Nose-frequenzy (Omega)   =  0.10E-29 period in steps =****** mass=  -0.146E-26a.u.
   SCALEE = 1.0000    scale energy and forces
   NPACO  =    256;   APACO  = 16.0  distance and # of slots for P.C.
   PSTRESS=    0.0 pullay stress

  Mass of Ions in am
   POMASS =  16.00  1.00
  Ionic Valenz
   ZVAL   =   6.00  1.00
  Atomic Wigner-Seitz radii
   RWIGS  =  -1.00 -1.00
  virtual crystal weights 
   VCA    =   1.00  1.00
   NELECT =       8.0000    total number of electrons
   NUPDOWN=      -1.0000    fix difference up-down

 DOS related values:
   EMIN   =  10.00;   EMAX   =-10.00  energy-range for DOS
   EFERMI =   0.00
   ISMEAR =     0;   SIGMA  =   0.20  broadening in eV -4-tet -1-fermi 0-gaus

 Electronic relaxation 2 (details)
   IALGO  =     38    algorithm
   LDIAG  =      T    sub-space diagonalisation (order eigenvalues)
   LSUBROT=      F    optimize rotation matrix (better conditioning)
   TURBO    =      0    0=normal 1=particle mesh
   IRESTART =      0    0=no restart 2=restart with 2 vectors
   NREBOOT  =      0    no. of reboots
   NMIN     =      0    reboot dimension
   EREF     =   0.00    reference energy to select bands
   IMIX   =      4    mixing-type and parameters
     AMIX     =   0.40;   BMIX     =  1.00
     AMIX_MAG =   1.60;   BMIX_MAG =  1.00
     AMIN     =   0.10
     WC   =   100.;   INIMIX=   1;  MIXPRE=   1;  MAXMIX= -45

 Intra band minimization:
   WEIMIN = 0.0010     energy-eigenvalue tresh-hold
   EBREAK =  0.20E-10  absolut break condition
   DEPER  =   0.30     relativ break condition  

   TIME   =   0.40     timestep for ELM

  volume/ion in A,a.u.               =     170.67      1151.72
  Fermi-wavevector in a.u.,A,eV,Ry     =   0.409275  0.773417  2.279054  0.167506
  Thomas-Fermi vector in A             =   1.364147
 
 Write flags
   LWAVE        =      T    write WAVECAR
   LDOWNSAMPLE  =      F    k-point downsampling of WAVECAR
   LCHARG       =      T    write CHGCAR
   LVTOT        =      F    write LOCPOT, total local potential
   LVHAR        =      F    write LOCPOT, Hartree potential only
   LELF         =      F    write electronic localiz. function (ELF)
   LORBIT       =      0    0 simple, 1 ext, 2 COOP (PROOUT), +10 PAW based schemes


 Dipole corrections
   LMONO  =      F    monopole corrections only (constant potential shift)
   LDIPOL =      F    correct potential (dipole corrections)
   IDIPOL =      0    1-x, 2-y, 3-z, 4-all directions 
   EPSILON=  1.0000000 bulk dielectric constant

 Exchange correlation treatment:
   GGA     =    --    GGA type
   LEXCH   =     8    internal setting for exchange type
   VOSKOWN=      0    Vosko Wilk Nusair interpolation
   LHFCALC =     F    Hartree Fock is set to
   LHFONE  =     F    Hartree Fock one center treatment
   AEXX    =    0.0000 exact exchange contribution

 Linear response parameters
   LEPSILON=     F    determine dielectric tensor
   LRPA    =     F    only Hartree local field effects (RPA)
   LNABLA  =     F    use nabla operator in PAW spheres
   LVEL    =     F    velocity operator in full k-point grid
   LINTERFAST=   F  fast interpolation
   KINTER  =     0    interpolate to denser k-point grid
   CSHIFT  =0.1000    complex shift for real part using Kramers Kronig
   OMEGAMAX=  -1.0    maximum frequency
   DEG_THRESHOLD= 0.2000000E-02 threshold for treating states as degnerate
   RTIME   =   -0.100 relaxation time in fs
  (WPLASMAI=    0.000 imaginary part of plasma frequency in eV, 0.658/RTIME)
   DFIELD  = 0.0000000 0.0000000 0.0000000 field for delta impulse in time
 
 Orbital magnetization related:
   ORBITALMAG=     F  switch on orbital magnetization
   LCHIMAG   =     F  perturbation theory with respect to B field
   DQ        =  0.001000  dq finite difference perturbation B field
   LLRAUG    =     F  two centre corrections for induced B field



--------------------------------------------------------------------------------------------------------


 finite differences with symmetry
 charge density and potential will be updated during run
 non-spin polarized calculation
 Variant of blocked Davidson
 Davidson routine will perform the subspace rotation
 perform sub-space diagonalisation
    after iterative eigenvector-optimisation
 modified Broyden-mixing scheme, WC =      100.0
 initial mixing is a Kerker type mixing with AMIX =  0.4000 and BMIX =      1.0000
 Hartree-type preconditioning will be used
 using additional bands          124
 reciprocal scheme for non local part
 calculate Harris-corrections to forces 
   (improved forces if not selfconsistent)
 use gradient corrections 
 use of overlap-Matrix (Vanderbilt PP)
 Gauss-broadening in eV      SIGMA  =   0.20


--------------------------------------------------------------------------------------------------------


  energy-cutoff  :      400.00
  volume of cell :      512.00
      direct lattice vectors                 reciprocal lattice vectors
     8.000000000  0.000000000  0.000000000     0.125000000  0.000000000  0.000000000
     0.000000000  8.000000000  0.000000000     0.000000000  0.125000000  0.000000000
     0.000000000  0.000000000  8.000000000     0.000000000  0.000000000  0.125000000

  length of vectors
     8.000000000  8.000000000  8.000000000     0.125000000  0.125000000  0.125000000


 
 k-points in units of 2pi/SCALE and weight: Gamma-point only                        
   0.00000000  0.00000000  0.00000000       1.000
 
 k-points in reciprocal lattice and weights: Gamma-point only                        
   0.00000000  0.00000000  0.00000000       1.000
 
 position of ions in fractional coordinates (direct lattice) 
   0.00000000  0.00000000  0.00000000
   0.07451015  0.90403665  0.00000000
   0.07451015  0.09596335  0.00000000
 
 position of ions in cartesian coordinates  (Angst):
   0.00000000  0.00000000  0.00000000
   0.59608120  7.23229320  0.00000000
   0.59608120  0.76770680  0.00000000
 


--------------------------------------------------------------------------------------------------------


 k-point  1 :   0.0000 0.0000 0.0000  plane waves:    9315

 maximum and minimum number of plane-waves per node :      9315     9315

 maximum number of plane-waves:      9315
 maximum index in each direction: 
   IXMAX=   13   IYMAX=   13   IZMAX=   13
   IXMIN=  -13   IYMIN=  -13   IZMIN=  -13


 serial   3D FFT for wavefunctions
 parallel 3D FFT for charge:
    minimum data exchange during FFTs selected (reduces bandwidth)


 total amount of memory used by VASP MPI-rank0    39891. kBytes
=======================================================================

   base      :      30000. kBytes
   nonl-proj :       1644. kBytes
   fftplans  :       1080. kBytes
   grid      :       7007. kBytes
   one-center:          9. kBytes
   wavefun   :        151. kBytes
 
     INWAV:  cpu time    0.0000: real time    0.0000
 Broyden mixing: mesh for mixing (old mesh)
   NGX = 27   NGY = 27   NGZ = 27
  (NGX  =108   NGY  =108   NGZ  =108)
  gives a total of  19683 points

 initial charge density was supplied:
 charge density of overlapping atoms calculated
 number of electron       8.0000000 magnetization 
 keeping initial charge density in first step


--------------------------------------------------------------------------------------------------------


 Maximum index for augmentation-charges          293 (set IRDMAX)


--------------------------------------------------------------------------------------------------------


 First call to EWALD:  gamma=   0.222
 Maximum number of real-space cells 3x 3x 3
 Maximum number of reciprocal cells 3x 3x 3

    FEWALD:  cpu time    0.0046: real time    0.0046


--------------------------------------- Iteration      1(   1)  ---------------------------------------


    POTLOK:  cpu time    4.2251: real time    4.4764
    SETDIJ:  cpu time    0.0037: real time    0.0038
     EDDAV:  cpu time   14.0357: real time   14.8238
       DOS:  cpu time    0.0037: real time    0.0040
    --------------------------------------------
      LOOP:  cpu time   18.2682: real time   19.3081

 eigenvalue-minimisations  :   256
 total energy-change (2. order) :-0.1903270E+01  (-0.2726502E+03)
 number of electron       8.0000000 magnetization 
 augmentation part        8.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.43793512
  Ewald energy   TEWEN  =        24.54851968
  -Hartree energ DENC   =      -385.76897518
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        32.23254366
  PAW double counting   =       247.56105404     -249.42875511
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -128.70605406
  atomic energy  EATOM  =       457.22046137
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =        -1.90327049 eV

  energy without entropy =       -1.90327049  energy(sigma->0) =       -1.90327049


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   2)  ---------------------------------------


     EDDAV:  cpu time   17.6461: real time   18.3860
       DOS:  cpu time    0.0007: real time    0.0007
    --------------------------------------------
      LOOP:  cpu time   17.6467: real time   18.3866

 eigenvalue-minimisations  :   384
 total energy-change (2. order) :-0.1434109E+02  (-0.1434109E+02)
 number of electron       8.0000000 magnetization 
 augmentation part        8.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.43793512
  Ewald energy   TEWEN  =        24.54851968
  -Hartree energ DENC   =      -385.76897518
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        32.23254366
  PAW double counting   =       247.56105404     -249.42875511
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -143.04714084
  atomic energy  EATOM  =       457.22046137
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -16.24435727 eV

  energy without entropy =      -16.24435727  energy(sigma->0) =      -16.24435727


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   3)  ---------------------------------------


     EDDAV:  cpu time    1.6789: real time    1.8616
       DOS:  cpu time    0.0006: real time    0.0006
    --------------------------------------------
      LOOP:  cpu time    1.6795: real time    1.8622

 eigenvalue-minimisations  :   256
 total energy-change (2. order) :-0.5414561E-02  (-0.5414561E-02)
 number of electron       8.0000000 magnetization 
 augmentation part        8.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.43793512
  Ewald energy   TEWEN  =        24.54851968
  -Hartree energ DENC   =      -385.76897518
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        32.23254366
  PAW double counting   =       247.56105404     -249.42875511
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -143.05255540
  atomic energy  EATOM  =       457.22046137
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -16.24977183 eV

  energy without entropy =      -16.24977183  energy(sigma->0) =      -16.24977183


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   4)  ---------------------------------------


     EDDAV:  cpu time    4.8905: real time    4.9724
       DOS:  cpu time    0.0004: real time    0.0004
    --------------------------------------------
      LOOP:  cpu time    4.8908: real time    4.9728

 eigenvalue-minimisations  :   512
 total energy-change (2. order) :-0.5554011E-05  (-0.5554010E-05)
 number of electron       8.0000000 magnetization 
 augmentation part        8.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.43793512
  Ewald energy   TEWEN  =        24.54851968
  -Hartree energ DENC   =      -385.76897518
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        32.23254366
  PAW double counting   =       247.56105404     -249.42875511
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -143.05256096
  atomic energy  EATOM  =       457.22046137
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -16.24977739 eV

  energy without entropy =      -16.24977739  energy(sigma->0) =      -16.24977739


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   5)  ---------------------------------------


     EDDAV:  cpu time    1.0164: real time    1.0295
       DOS:  cpu time    0.0004: real time    0.0004
    CHARGE:  cpu time    0.0512: real time    0.0885
    MIXING:  cpu time    0.0088: real time    0.0117
    --------------------------------------------
      LOOP:  cpu time    1.0768: real time    1.1301

 eigenvalue-minimisations  :   256
 total energy-change (2. order) :-0.7389644E-11  (-0.5698553E-11)
 number of electron       7.9999967 magnetization 
 augmentation part        0.8297670 magnetization 


--------------------------------------------------------------------------------------------------------


--------------------------------------------------------------------------------------------------------


 Finite differences progress:
  Degree of freedom:   6/  6
  Displacement:        2/  2
  Total:              12/ 12


--------------------------------------------------------------------------------------------------------


 INTERNAL STRAIN TENSORS FROM DISPLACED ATOMS
 ============================================

 INTERNAL STRAIN TENSOR FOR ION    1 for displacements in x,y,z  (eV/Angst):
          X           Y           Z          XY          YZ          ZX
  --------------------------------------------------------------------------------
  x    24.16255    28.30049     0.03786     0.00000     0.00000     0.00000
  y     0.00000     0.00000     0.00000    37.03527     0.00000     0.00000
  z     0.00000     0.00000     0.00000     0.00000     0.00000    -0.22355

 INTERNAL STRAIN TENSOR FOR ION    2 for displacements in x,y,z  (eV/Angst):
          X           Y           Z          XY          YZ          ZX
  --------------------------------------------------------------------------------
  x   -12.08592   -14.15263    -0.01906    14.31295     0.00000     0.00000
  y    10.95335    28.22858    -0.00604   -18.52124     0.00000     0.00000
  z     0.00000     0.00000     0.00000     0.00000    -0.14689     0.11204

 INTERNAL STRAIN TENSOR FOR ION    3 for displacements in x,y,z  (eV/Angst):
          X           Y           Z          XY          YZ          ZX
  --------------------------------------------------------------------------------
  x   -12.08592   -14.15263    -0.01906   -14.31295     0.00000     0.00000
  y   -10.95335   -28.22858     0.00604   -18.52124     0.00000     0.00000
  z     0.00000     0.00000     0.00000     0.00000     0.14689     0.11204


--------------------------------------------------------------------------------------------------------


 
 Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------
 
 
   1 f  =  115.029900 THz   722.754181 2PiTHz 3836.984350 cm-1   475.725600 meV
             X         Y         Z           dx          dy          dz
      0.000000  0.000000  0.000000    -0.000000   -0.268872    0.000000  
      0.596081  7.232293  0.000000    -0.417632    0.537994   -0.000000  
      0.596081  0.767707  0.000000     0.417632    0.537994   -0.000000  
 
   2 f  =  111.640489 THz   701.457881 2PiTHz 3723.925763 cm-1   461.708116 meV
             X         Y         Z           dx          dy          dz
      0.000000  0.000000  0.000000     0.192070   -0.000000   -0.000000  
      0.596081  7.232293  0.000000    -0.384335    0.577790    0.000000  
      0.596081  0.767707  0.000000    -0.384335   -0.577790   -0.000000  
 
   3 f  =   47.399735 THz   297.821318 2PiTHz 1581.084921 cm-1   196.029617 meV
             X         Y         Z           dx          dy          dz
      0.000000  0.000000  0.000000     0.271617   -0.000000    0.000000  
      0.596081  7.232293  0.000000    -0.544935   -0.407625    0.000000  
      0.596081  0.767707  0.000000    -0.544935    0.407625   -0.000000  
 
   4 f  =    0.571381 THz     3.590091 2PiTHz   19.059211 cm-1     2.363042 meV
             X         Y         Z           dx          dy          dz
      0.000000  0.000000  0.000000     0.000000    0.000000    0.885951  
      0.596081  7.232293  0.000000     0.000000    0.000000    0.327942  
      0.596081  0.767707  0.000000    -0.000000    0.000000    0.327942  
 
   5 f  =    0.261288 THz     1.641723 2PiTHz    8.715638 cm-1     1.080602 meV
             X         Y         Z           dx          dy          dz
      0.000000  0.000000  0.000000     0.000000   -0.930287    0.000000  
      0.596081  7.232293  0.000000    -0.032115   -0.257393    0.000000  
      0.596081  0.767707  0.000000     0.032115   -0.257393    0.000000  
 
   6 f/i=    0.039361 THz     0.247310 2PiTHz    1.312929 cm-1     0.162783 meV
             X         Y         Z           dx          dy          dz
      0.000000  0.000000  0.000000     0.943045    0.000000   -0.000000  
      0.596081  7.232293  0.000000     0.235230   -0.000274   -0.000000  
      0.596081  0.767707  0.000000     0.235230    0.000274   -0.000000  
 
   7 f/i=    3.984794 THz    25.037197 2PiTHz  132.918407 cm-1    16.479788 meV
             X         Y         Z           dx          dy          dz
      0.000000  0.000000  0.000000    -0.000000    0.000000    0.463779  
      0.596081  7.232293  0.000000    -0.000000   -0.000000   -0.626462  
      0.596081  0.767707  0.000000     0.000000   -0.000000   -0.626462  
 
   8 f/i=    4.355057 THz    27.363629 2PiTHz  145.269056 cm-1    18.011074 meV
             X         Y         Z           dx          dy          dz
      0.000000  0.000000  0.000000     0.000000   -0.000000   -0.000000  
      0.596081  7.232293  0.000000     0.000000    0.000000    0.707107  
      0.596081  0.767707  0.000000    -0.000000    0.000000   -0.707107  
 
   9 f/i=    4.648552 THz    29.207716 2PiTHz  155.059014 cm-1    19.224874 meV
             X         Y         Z           dx          dy          dz
      0.000000  0.000000  0.000000    -0.000000   -0.249546    0.000000  
      0.596081  7.232293  0.000000     0.569695    0.379883   -0.000000  
      0.596081  0.767707  0.000000    -0.569695    0.379883    0.000000  
 

 ELASTIC MODULI CONTR FROM IONIC RELAXATION (kBar)
 Direction    XX          YY          ZZ          XY          YZ          ZX
 --------------------------------------------------------------------------------
 XX         -44.5872    -52.2500     -0.0700      0.0000     -0.0000      0.0000
 YY         -52.2500   -134.4072      0.0284     -0.0000     -0.0000     -0.0000
 ZZ          -0.0700      0.0284     -0.0003     -0.0000      0.0000     -0.0000
 XY           0.0000      0.0000      0.0000    -68.4799      0.0000      0.0000
 YZ          -0.0000     -0.0000      0.0000      0.0000     -0.0000     -0.0000
 ZX           0.0000     -0.0000     -0.0000      0.0000     -0.0000     -0.8904
 --------------------------------------------------------------------------------



--------------------------------------------------------------------------------------------------------


  LATTYP: Found a simple cubic cell.
 ALAT       =     8.0000000000
  
  Lattice vectors:
  
 A1 = (   8.0000000000,   0.0000000000,   0.0000000000)
 A2 = (   0.0000000000,   8.0000000000,   0.0000000000)
 A3 = (   0.0000000000,   0.0000000000,   8.0000000000)


Analysis of symmetry for initial positions (statically):
=====================================================================
 Subroutine PRICEL returns:
 Original cell was already a primitive cell.
 

 Routine SETGRP: Setting up the symmetry group for a 
 simple cubic supercell.


 Subroutine GETGRP returns: Found  4 space group operations
 (whereof  4 operations were pure point group operations)
 out of a pool of 48 trial point group operations.


The static configuration has the point symmetry C_2v.


Analysis of symmetry for dynamics (positions and initial velocities):
=====================================================================
 Subroutine PRICEL returns:
 Original cell was already a primitive cell.
 

 Routine SETGRP: Setting up the symmetry group for a 
 simple cubic supercell.


 Subroutine GETGRP returns: Found  4 space group operations
 (whereof  4 operations were pure point group operations)
 out of a pool of 48 trial point group operations.


The dynamic configuration has the point symmetry C_2v.


 Subroutine INISYM returns: Found  4 space group operations
 (whereof  4 operations are pure point group operations),
 and found     1 'primitive' translations

 
 KPOINTS: Gamma-point only                        

Automatic generation of k-mesh.
Space group operators:
 irot       det(A)        alpha          n_x          n_y          n_z        tau_x        tau_y        tau_z
    1     1.000000     0.000000     1.000000     0.000000     0.000000     0.000000     0.000000     0.000000
    2    -1.000000   180.000000     0.000000     0.000000     1.000000     0.000000     0.000000     0.000000
    3     1.000000   180.000000     1.000000     0.000000     0.000000     0.000000     0.000000     0.000000
    4    -1.000000   180.000000     0.000000     1.000000     0.000000     0.000000     0.000000     0.000000
 
 Subroutine IBZKPT returns following result:
 ===========================================
 
 Found      1 irreducible k-points:
 
 Following reciprocal coordinates:
            Coordinates               Weight
  0.000000  0.000000  0.000000      1.000000
 
 Following cartesian coordinates:
            Coordinates               Weight
  0.000000  0.000000  0.000000      1.000000
 
 writing wavefunctions
     LOOP+:  cpu time   53.9186: real time   55.2177
    4ORBIT:  cpu time    0.0000: real time    0.0000

 total amount of memory used by VASP MPI-rank0    39891. kBytes
=======================================================================

   base      :      30000. kBytes
   nonl-proj :       1644. kBytes
   fftplans  :       1080. kBytes
   grid      :       7007. kBytes
   one-center:          9. kBytes
   wavefun   :        151. kBytes
 
  
  
 General timing and accounting informations for this job:
 ========================================================
  
                  Total CPU time used (sec):      584.746
                            User time (sec):      564.399
                          System time (sec):       20.347
                         Elapsed time (sec):      607.447
  
                   Maximum memory used (kb):      380764.
                   Average memory used (kb):           0.
  
                          Minor page faults:       206805
                          Major page faults:            0
                 Voluntary context switches:        15485
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(outcar_content)
            outcar_path = f.name

        yield outcar_path

        # テスト後に一時ファイルを削除
        os.unlink(outcar_path)

    def test_load_from_outcar(self):
        """OUTCARからの読み込みテスト"""
        analysis = read_file(os.path.split(__file__)[0]+"/OUTCAR")  # H2O
        modes = analysis.modes

        assert len(modes) == 9

        # 実数周波数のテスト
        assert modes[0].frequency == 3836.98435
        assert modes[1].frequency == 3723.925763
        assert modes[2].frequency == 1581.084921

    def test_exclude_imag(self):
        """虚数周波数除外オプションのテスト"""
        analysis = read_file(os.path.split(__file__)[
                             0]+"/OUTCAR", exclude_imag=True)
        modes = analysis.modes

        assert len(modes) == 5
        assert all(mode.frequency >= 0 for mode in modes)

    def test_get_mode(self):
        """test get_mode"""
        analysis = read_file(os.path.split(__file__)[0]+"/OUTCAR")
        mode = analysis.get_mode(0)
        assert mode.frequency == 3836.98435

        # 範囲外のインデックスでの例外テスト
        with pytest.raises(IndexError):
            analysis.get_mode(10)

    def test_get_frequencies(self):
        """get_frequenciesメソッドのテスト"""
        analysis = read_file(os.path.split(__file__)[0]+"/OUTCAR")

        freqs = analysis.get_frequencies()
        assert len(freqs) == 9
        assert freqs[0] == 3836.98435
        assert freqs[-1] == -155.059014

    def test_get_real_modes(self):
        """実数モードの取得テスト"""
        analysis = read_file(os.path.split(__file__)[0]+"/OUTCAR")
        modes = analysis.modes

        real_modes = [mode for mode in modes if mode.frequency >= 0]
        assert len(real_modes) == 5
        assert real_modes[0].frequency == 3836.98435
        assert real_modes[1].frequency == 3723.925763

    def test_get_imaginary_modes(self):
        """虚数モードの取得テスト"""
        analysis = read_file(os.path.split(__file__)[0]+"/OUTCAR")
        modes = analysis.modes

        imag_modes = [mode for mode in modes if mode.frequency < 0]
        assert len(imag_modes) == 4
        assert imag_modes[0].frequency == -1.312929

    def test_get_eigenvectors(self):
        """固有ベクトルの取得テスト"""
        analysis = read_file(os.path.split(__file__)[0]+"/OUTCAR")

        eigenvectors = analysis.get_eigenvectors()
        assert eigenvectors.shape == (9, 3, 3)  # 3モード, 3イオン, 3次元

        # 最初のモードの最初のイオンの変位ベクトルをチェック
        assert np.array_equal(eigenvectors[0, 0], [-0., -0.268872,  0.])

# 修正が必要なバグを見つけるテスト


def test_class_fix_methods():
    """クラス実装のバグを見つけるためのテスト"""
    # get_real_modesとget_imaginary_modesメソッドのバグを検出
    atoms = ase.Atoms(["H"])
    analysis = VibrationAnalysis()
    analysis.modes = [
        VibrationalMode(atoms, 1.0, np.ones((1, 3))),
        VibrationalMode(atoms, -2.0, np.ones((1, 3)))
    ]

    # バグ：selfの参照が正しくない
    # 以下のテストは現在のコードでは失敗するが、修正後は成功する
    try:
        real_modes = analysis.get_real_modes()
        assert len(real_modes) == 1
        assert real_modes[0].frequency == 1.0

        imag_modes = analysis.get_imaginary_modes()
        assert len(imag_modes) == 1
        assert imag_modes[0].frequency == -2.0
    except NameError:
        # modesが定義されていないエラーが出るはず
        print("バグ検出: get_real_modes()とget_imaginary_modes()メソッドでselfを参照していません")
        print("修正案: 'modes'を'self.modes'に変更してください")


def test_save_xsf():
    """save_xsf関数のテスト"""
    # テスト用のASE原子オブジェクトを作成
    from ase import Atoms

    atoms = Atoms('H2O',
                  positions=[[0.0, 0.0, 0.0],
                             [0.0, 0.77, 0.59],
                             [0.0, 0.77, -0.59]],
                  cell=[[3.0, 0.0, 0.0],
                        [0.0, 3.0, 0.0],
                        [0.0, 0.0, 3.0]])

    eigenvector = np.array([[0.1, 0.0, 0.0],
                            [0.1, 0.0, 0.0],
                            [0.1, 0.0, 0.0]])

    frequency = np.array([5, 6, 7])

    vibmode = VibrationalMode(atoms, frequency, eigenvector)
    vibanalysis = VibrationAnalysis()
    vibanalysis.set_params([vibmode])

    # create temporary files
    with tempfile.NamedTemporaryFile(delete=False, suffix='.xsf') as temp_file:
        temp_filename = temp_file.name

    try:
        result = vibanalysis.save_xsf(temp_filename, 0, scale=2.0)
        assert os.path.exists(temp_filename)  # check if file exists
        # read & check the content
        with open(temp_filename, 'r') as f:
            file_content = f.read()
        assert result == file_content
        assert "CRYSTAL" in result
        assert "PRIMVEC" in result
        assert "PRIMCOORD" in result
        assert "3 1" in result
        assert "0.2000000000000000" in result

    finally:
        # remove test files
        if os.path.exists(temp_filename):
            os.remove(temp_filename)


def test_save_xsf2():
    """save_xsf関数のテスト"""
    # テスト用のASE原子オブジェクトを作成
    from ase import Atoms

    atoms = Atoms('H2O',
                  positions=[[0.0, 0.0, 0.0],
                             [0.0, 0.77, 0.59],
                             [0.0, 0.77, -0.59]],
                  cell=[[3.0, 0.0, 0.0],
                        [0.0, 3.0, 0.0],
                        [0.0, 0.0, 3.0]])

    eigenvector = np.array([[0.1, 0.0, 0.0],
                            [0.1, 0.0, 0.0],
                            [0.1, 0.0, 0.0]])

    frequency = np.array([5, 6, 7])

    vibmode = VibrationalMode(atoms, frequency, eigenvector)

    # create temporary files
    with tempfile.NamedTemporaryFile(delete=False, suffix='.xsf') as temp_file:
        temp_filename = temp_file.name

    try:
        result = save_xsf(temp_filename, vibmode, scale=2.0)
        assert os.path.exists(temp_filename)  # check if file exists
        # read & check the content
        with open(temp_filename, 'r') as f:
            file_content = f.read()
        assert result == file_content
        assert "CRYSTAL" in result
        assert "PRIMVEC" in result
        assert "PRIMCOORD" in result
        assert "3 1" in result
        assert "0.2000000000000000" in result

    finally:
        # remove test files
        if os.path.exists(temp_filename):
            os.remove(temp_filename)


if __name__ == "__main__":
    pytest.main(["-v"])
