!     HELIOS simulations.  PA Butler  April 2014.  Modified for active gas volume mode June 2022.
!     Fortran version created February 2025

      PROGRAM HELIOS
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Q_store_pass(100000), Q_store_loops(1000000)
      
      COMMON/stopping/energy4(1000), dedx4(1000), pi, ii_max
      
      REAL*8 M1, M2, M3, M4, mass
      REAL*8 M1_stop, M4_stop, inhomo_Bz, inhomo_Bz_percent
      REAL*8 :: multi_scatter, terpolate, gauss

      CHARACTER*200 dummy, M4_stop_file, M4_stop_filename, folder_name, reaction_label, file_in, file_out
      CHARACTER*1 creturn
      INTEGER beam_stop, eject_stop 
      INTEGER sintheta_weighting, active_flag, ii_eject_scat
      INTEGER cent, plot
      INTEGER error_traj_low, error_det_overrun, error_slow, error_large_radius, error_dedx, error_E4low 
      INTEGER error_miss_detector, error_no_zvertex, error_no_converge, error_poor_Q, loops, good_loops
      INTEGER pass_through, good_pass, check_pad, radius_i, pitch_x, pitch_z, pitch_radius_pad
      INTEGER xi, zi, kz
      INTEGER, parameter :: seed = 86456
	        
      creturn = ACHAR(13) 
  	  
      sintheta_weighting = 0.
      plot = 0
      ii_max = 0
      pi = 4.*ATAN(1.)

      CALL SRAND(seed)

      OPEN(4,FILE = 'reaction.txt')
      READ(4,100)folder_name
      READ(4,100)reaction_label

100   FORMAT(A80)

      file_in = TRIM(folder_name) // '/input_' // TRIM(reaction_label) // '.txt'
      file_out = TRIM(folder_name) // '/output_' // TRIM(reaction_label) // '.txt'


      OPEN(1, FILE = file_in)
      OPEN(2, FILE = file_out)


      Ez = 2500. / 0.15
! electric field if active target (V/m)

      DO 

       READ(1, 100, END = 302)dummy
       READ(1,*)M1, M2, M4, EoA, Qvalue

       READ(1, 100)dummy
       READ(1,*)Bz, charge_state, Rmagnetcm, inhomo_Bz_percent

       READ(1, 100)dummy
       READ(1,*)Rdetcm, z1detcm, z2detcm

       READ(1, 100)dummy
       READ(1,*)pitch_zmm, pitch_xmm, fwhm_E4keV, E4_thresh_MeV

       READ(1, 100)dummy
       READ(1,*)target, M1_stop, M4_stop_file, beam_stop, eject_stop

       READ(1, 100)dummy
       READ(1,*)ii_scat_max, Z_beam, ii_eject_scat, fwhm_beam_straggle_percent, fwhm_eject_straggle_percent
    
       READ(1, 100)dummy
       READ(1,*)fwhm_beam_spot_mm, fwhm_E0_percent, fwhm_diverg_mrad

       READ(1, 100)dummy
       READ(1,*)thetacm_step

       IF (thetacm_step.EQ.0.) THEN
 
        thetacm = 35.
        thetacm_step = 180.
 
       ENDIF
    
       IF (thetacm_step.LT.0.) THEN
 
        thetacm = -thetacm_step
        thetacm_step = 180
 
       ENDIF
    

    
    
       M4_stop_filename = 'srim_files/' // TRIM(M4_stop_file)
    
       OPEN(3, FILE = M4_stop_filename)

       Bx = 0
       By = 0
       Z_target = 0.
       A_target = 0.
    
       IF ((M4_stop_file.EQ."proton_srim.txt").OR.(M4_stop_file.EQ."deuterium_srim.txt").OR.(M4_stop_file.EQ."carbon_srim.txt") &
            .OR.(M4_stop_file.EQ."proton_on_isobutane_srim.txt")) THEN
        Z_target = 6.
        A_target = 12.
       ENDIF
    
       IF ((M4_stop_file.EQ."proton_on_CD2_srim.txt").OR.(M4_stop_file.EQ."deuteron_on_CD2_srim.txt").OR.(M4_stop_file.EQ. &
          "carbon_on_carbon_srim.txt")) THEN
        Z_target = 6.
        A_target = 12.
       ENDIF
    
       IF (M4_stop_file.EQ."proton_on_deuterium_srim.txt") THEN
        Z_target = 1.
        A_target = 2.
       ENDIF
    
       IF ((M4_stop_file.EQ."deuteron_on_deuterium_srim.txt").OR.(M4_stop_file.EQ."He-3_on_deuterium_srim.txt")) THEN
        Z_target = 1.
        A_target = 2.
       ENDIF
    
       IF ((M4_stop_file.EQ."alpha_on_helium_srim.txt").OR.(M4_stop_file.EQ."alpha_on_He-CO2_srim.txt")) THEN
        Z_target = 2.
        A_target = 4.
       ENDIF
    
       IF ((M4_stop_file.EQ."proton_on_helium_srim.txt").OR.(M4_stop_file.EQ."proton_on_He-CO2_srim.txt")) THEN
        Z_target = 2.
        A_target = 4.
       ENDIF
    
       IF (Z_target.EQ.0.) THEN
       WRITE(2, 101)
101     FORMAT("target not specified")
        GOTO 1001
       ENDIF
    
       scatter_factor_active = 1.8
    
       active_flag = 0
       IF (target .GT. 1) THEN 
        active_flag = 1
       ENDIF

! if value of target thickness is greater than 1, then this value is assumed to be the gas pressure (in torr) of an active target
        
       density_760 = 0.
!  this statement has been inserted February 2025
	   
       IF (active_flag.EQ.1) THEN
        s_travel_increment = 0.01
        
        
        
        IF (ii_eject_scat.LT.0) THEN
            scatter_factor_active = 1.0
            ii_eject_scat = 1
        ENDIF
!       if scatter flag < 0 do not modify multiple scattering for gas target
         
          
    
        IF (Z_target.EQ.1. .AND. A_target.EQ.2.) THEN
            density_760 = 0.1702 / 1000
        ENDIF
!       deuterium
    
        IF (Z_target.EQ.2. .AND. A_target.EQ.4.) THEN
            density_760 = 0.1662 / 1000
        ENDIF
!       helium
    
        IF (Z_target.EQ.6. .AND. A_target .EQ. 12.) THEN
            density_760 = 2.53 / 1000
            s_travel_increment = 0.001
        ENDIF
!       isobutane
    
        IF (density_760.EQ.0.) THEN
         WRITE(2,102)
102      FORMAT("gas density not specified")
         GOTO 1001
        END IF
    
       ENDIF

!    gas density in g/cm^3 at 15C, atmospheric pressure

    

       density = density_760 * target / 760.

!   actual gas density in g/cm^3
   


       
       ii = 0
       DO 
        ii = ii + 1
        READ(3,*, END = 301)energy4(ii), dedx4(ii)
       END DO
	   
301    ii_max = ii

       CLOSE(3)

       sigma_beam_straggle = fwhm_beam_straggle_percent / 235.
       sigma_beam_spot = fwhm_beam_spot_mm / 2350.
       sigma_eject_straggle = fwhm_eject_straggle_percent / 235.




       Rmagnet = Rmagnetcm / 100.

! radius of magnet
    
       scatter_factor = 1.
	   
       z_interact_startcm = 0.
!    this line inserted February 2025
    
       IF (active_flag.EQ.1) THEN
    
        scatter_factor = scatter_factor_active
    
! this emperically adjusts the value of sigma for multiple scattering so that the Kantele subroutine
! agrees with the Fano prescription, which reproduces the data of Kuhn et al. NIM B4 (1984) 332
    
         z_interact_startcm = z2detcm
     
         IF (z1detcm.LT.0.) THEN
            z2detcm = z1detcm - 50.
            window_positioncm = z1detcm - 5.
         ELSE
            z2detcm = z1detcm + 50.
            window_positioncm = -5.
         ENDIF

! window_positioncm is position of gas window, in cm
 
       ENDIF
    
! for active target the Si detector is assumed to be 50cm long
! in this case the interaction region is between 0 and z_interact_startcm
    
       window_thickness = 0.1

! for active target this is window foil thickness in mg/cm**2
    

   
    
       Rpadcm = 5.5
       cent = 0

       pitch_rpadmm = 0.
!   this line inserted February 2025	   

       IF (cent.EQ.0 .AND. active_flag.EQ.1) THEN
        Rpadcm = Rdetcm
        Rdetcm = 3.
       ENDIF
    
       Rdet = Rdetcm / 100.

! radius of Si detector (m)
    
       Rpad = Rpadcm / 100.

! radius of inner pads (m)
        

       IF (active_flag.EQ.1) THEN
        FWHM_z_padmm = pitch_zmm
        pitch_rpadmm = pitch_xmm
        pitch_zmm = 0.95
        pitch_xmm = 2.0
       ENDIF
    
       IF (FWHM_z_padmm.LT.0.02 .AND. active_flag.EQ.1) THEN

        pitch_zmm = 0.01
        pitch_xmm = 0.01
        pitch_rpadmm = 0.01

       ENDIF
    
! for active target the Si detector radius, pitch in z and x direction are fixed
! in this case the uncertainty in the z-value measured by the pad detectors is FWHM_z_padmm
! if FWHM_z_padmm is very small (0.01mm) then it is assumed that pitch_zmm and pitch_xmm are also very small

       inhomo_Bz = inhomo_Bz_percent / 100. * Bz

       z1det = z1detcm / 100.
       z2det = z2detcm / 100.

! front position and rear position of Si detector wrt target


       pitch_z = INT(pitch_zmm * 100.)

! pitch of z strips of Si detector


       pitch_x = INT(pitch_xmm * 100.)

! pitch of x strips of Si detector

       pitch_radius_pad = INT(pitch_rpadmm * 100.)
    
! pitch of pad sensors
    

       sigma_E4 = fwhm_E4keV / 2350.

! uncertainty in ejectile energy measurement

      
!       pi = 3.1415927
       z_factor = SQRT(2. * 1.60218E-13 / (1.66054E-27)) * 1.66054E-27 * 2. * pi / 1.60218E-19

       M3 = M1 + M2 - M4
       E0beam = EoA * M1
    
       ET0 = E0beam + Qvalue
       factor1 = (1. + M1 / M2 * Qvalue / ET0)
       factor2 = (M1 + M2) * (M3 + M4)
       A0 = M1 * M4 / factor2 * E0beam / ET0
       B0 = M1 * M3 / factor2 * E0beam / ET0
       C0 = M2 * M3 / factor2 * factor1
       D0 = M2 * M4 / factor2 * factor1
    
    
    
    
       VV = SQRT(E0beam**2 + 2 * E0beam * M1 * 931.36814) / (E0beam + M1 * 931.36814) * 2.9979246E+10
       XX = VV / (3.6E+08 * Z_beam**0.45)
       YY = 1 - .00119 * (Z_target - 6) * SQRT(XX) + .00001 * (Z_target - 6)**2 * XX
       charge_state_beam = Z_beam * (1 - EXP(-1.25 * XX + .32 * XX**2 - .11 * XX**3)) * YY
    

       sigma_E0beam = fwhm_E0_percent / 235. * E0beam
       sigma_angle_deg = fwhm_diverg_mrad / 2.35 * (180. / (1000. * pi))


       Vcm = M1 / (M1 + M2) * 4.633981637 * SQRT(E0beam / M1) * 2997924.6
       factor4 = M4 * 1.6603145E-27 * (Vcm * Vcm) / 2. * 6.241509074E12
       rho_max_factor = 2. * SQRT(2.) * SQRT(1.602176634E-13 * 1.6603145E-27) / 1.602176634E-19
       Tcyc_factor = 2 * pi * 1.6603145E-27 / 1.602176634E-19 * 1.00000E9
       Tcyc = Tcyc_factor * M4 * 1.00000E-9 / (charge_state * Bz)
       Helios_factor = Tcyc / (1.6603145E-27 * M4 * Vcm * 6.241509074E12)
       factor5 = 1. / Helios_factor / 100.
       AA = (M3 + M4) / M3
       BB = factor4 * AA - M2 * E0beam / (M1 + M2)
       CC = factor5 * AA
    
!         WRITE "z_factor"; z_factor; "   rho_max_factor"; rho_max_factor


    Tcyc_ns = Tcyc * 1E9

       WRITE(2,103)M1, M2, M3, M4, EoA, Qvalue
103    FORMAT("M1 = ",F4.0,"  M2 =",F4.0,"  M3 =",F4.0,"  M4 =",F4.0,"  Beam energy =",F5.2," MeV/A", "  Q value =",F5.1," MeV")

       WRITE(2,104)Bz, charge_state, Tcyc_ns
104    FORMAT("Bz =",F5.1," T", "  q =",F5.0," T cyc =",F6.3," ns")

       WRITE(2,105)AA, BB, CC
105    FORMAT("Qvalue (MeV) = ",F8.4, " *E4 (MeV) +",F8.4," - ",F8.4," *z (cm)")


       mass = M4 * 1.66054E-27

       q = charge_state * 1.60218E-19
      
   
       deltat = 1.0E-11

!  time interval for path integral (s)






 
    
       IF (active_flag.EQ.0) THEN

        WRITE(2,106)Rdetcm, z1detcm, z2detcm, Rmagnetcm, inhomo_Bz_percent  
106     FORMAT("radius Si detector =",F5.1,"cm", "  detector between ",F5.1," and ",F5.1,"cm", "    radius magnet =", F5.1,"cm",  &
             "    inhomogeneity in Bz =",F5.1,"%")
        WRITE(2,107)pitch_zmm, pitch_xmm, fwhm_E4keV, E4_thresh_MeV
107     FORMAT("Si pitch z =",F5.2,"mm", "    pitch x =",F5.2,"mm", "     Si FWHM E4 =",F5.2," keV    E4 threshold =",F5.2," MeV")
        WRITE(2,151) target
151     FORMAT("target thickness =", F5.2, " mg/cm**2",$)
        IF(ii_scat_max.GT.1)THEN
            WRITE(2,108)ii_scat_max
108         FORMAT("  no. of targets =",I5)
        ELSE
            WRITE(2,109)
109         FORMAT("  no. of targets = 1")
        ENDIF

       ELSE

        WRITE(2,112)Rpadcm, Rdetcm, z1detcm, z2detcm, Rmagnetcm, inhomo_Bz_percent
112     FORMAT("inner radius pad detector =",F6.1,"cm", "   radius Si detector =",F6.1,"cm", "  detector between ",F6.1,  &
              " and ",F6.1,"cm", "    radius magnet =",F6.1,"cm", "    inhomogeneity in Bz =",F6.1,"%")
        WRITE(2,113)pitch_zmm, pitch_xmm, FWHM_z_padmm, pitch_rpadmm, fwhm_E4keV, E4_thresh_MeV
113     FORMAT("Si pitch z =",F5.2,"mm", "    pitch x =",F5.2,"mm","     pad FWHM z =",F6.1,"mm","    pitch pad =",F6.1, &
                 "mm", "     Si FWHM E4 =",F6.1,"keV    E4 threshold =",F6.1,"MeV")
        WRITE(2,114)target, z_interact_startcm
114     FORMAT("gas pressure =", F6.1,"torr", "  between 0 and ",F6.1,"cm")

       ENDIF

       WRITE(2,110) M1_stop, M4_stop_file
110    FORMAT("beam stopping =",F6.1," MeV/mg/cm**2", "     ejectile stopping power file: ",A30,$)
       WRITE(2,111) beam_stop, eject_stop
111    FORMAT("  beam stopping flag =", I2,"  ejectile stopping flag =",I2)
       WRITE(2,115)Z_beam
115    FORMAT("Z beam =",F4.0,$)
       WRITE(2,116)charge_state_beam;
116    FORMAT("  beam charge state = ", F4.1,$)
       WRITE(2,117)Z_target, A_target
117    FORMAT("  Z target =",F4.0,"  A target =",F4.0,$)
       WRITE(2,118)ii_eject_scat
118    FORMAT("   ejectile multiple scattering flag =",I2,$)
       IF(ii_scat_max.GT.0)THEN
        WRITE(2,119) 
119     FORMAT("  beam multiple scattering flag = 1")
       ELSE
        WRITE(2,120) 
120     FORMAT("  beam multiple scattering flag = 0")
       ENDIF
       WRITE(2,121)scatter_factor
121    FORMAT("multiple scattering factor =", F6.1)  
       WRITE(2,122)fwhm_beam_straggle_percent, fwhm_eject_straggle_percent
122    FORMAT("FWHM beam straggling =",F5.0,"% of energy loss", "  FWHM ejectile straggling = ", F5.0, "% of energy loss")
       WRITE(2,123)fwhm_beam_spot_mm, fwhm_E0_percent, fwhm_diverg_mrad
123    FORMAT("FWHM beam spot =", F5.1, "mm", "  FWHM beam energy spread =", F5.1, "%", &
               "  FWHM beam divergence =", F5.1,"mrad")
       IF (active_flag .EQ. 1) THEN
        WRITE(2,124)window_thickness, window_positioncm
124     FORMAT("window foil thickness =", F5.2,"mg/cm^2","   window position = ", F6.1,"cm")
       ENDIF
    
       v_sum = 0.
       v_diff_sq_sum = 0.



    
       n_pass = 1
       IF (thetacm_step.GT.1. .AND. plot.EQ.0) THEN
        n_pass = INT(thetacm_step)
        thetacm_step = 0.2
       ENDIF
    
       IF (plot.EQ.0) THEN 
        thetacm = thetacm_step
       ENDIF
    
    
    
!  if thetacm_step is between 0 and 1 then step over theta
!  if thetacm_step > 1 then make multiple passes (e.g. to determine efficiency)
!  if thetacm_step = 0 then plot trajectory for theta = 30 degrees
!  if thetacm_step < 0 then plot trajectories at theta = |thetacm_step|
    
       WRITE(2,125)
125    FORMAT(" ")
       WRITE(2,126)n_pass
126    FORMAT("number of passes = ", I5)
    
       WRITE(2,127) 
127    FORMAT("theta cm     E4    zeta     t      z     Q value   Dz (pad) (Si)  vertex  corr. t       z     Q value ",  &
               "Delta Q rho_max rad_max distance E4 loss theta  efficiency")
       WRITE(2,128)
128    FORMAT(" (deg.)     (MeV) (deg.)   (ns)   (cm)    (MeV)     (cm)    (cm)    (cm)       (ns)    (cm)    (MeV)", &
               "    (keV)    (cm)    (cm)    (cm)   (MeV)   (deg)   (%)")
    
       error_traj_low = 0
       error_det_overrun = 0
       error_slow = 0
       error_large_radius = 0
       error_dedx = 0
       error_E4low = 0
       error_miss_detector = 0
       error_no_zvertex = 0
       error_no_converge = 0
       error_poor_Q = 0
       loops = 0
       good_loops = 0
    

       DO WHILE (thetacm < 175.)
    
        WRITE(*,352, ADVANCE='NO' ) creturn, thetacm
352     FORMAT(A, F4.1)

        pass_through = 0
        good_pass = 0
        
        E4_sum = 0.
        zeta_sum = 0.
        phi_deg_sum = 0.
        tns_sum = 0.
        z_vertexcm_sum = 0.
        zcm_sum = 0.
        Qvalue_uncorr_sum = 0.
        z_corrcm_sum = 0.
        z_changecm_sum = 0.
        z_change_padcm_sum = 0.
        Qvalue_corr_sum = 0.
        tns_corr_sum = 0.
        rho_max_cm_sum = 0.
        s_travel_tot_sum = 0.
        radius_max = 0.
        E4_loss_sum = 0.
        thetacm_meas_sum = 0.

!    Qvalue_corr = Qvalue
        
    
        DO WHILE (pass_through < n_pass)
    
            pass_through = pass_through + 1
            loops = loops + 1
                  
            E0 = gauss(E0beam, sigma_E0beam)

!    randomise beam energy because of energy spread
        

            IF(active_flag.EQ.0)THEN
                z_vertexcm = 0
            
                tt = RAND() * target
                IF (beam_stop .EQ. 1) THEN 
                 E0 = E0 - M1_stop * tt
                ENDIF

!     beam loses energy randomly in the target


                IF (target .GT. 0) THEN
                    ii_target = INT(RAND() * FLOAT(ii_scat_max)) + 1
                    IF (ii_target.GT.ii_scat_max) THEN 
                     ii_target = ii_scat_max
                    ENDIF
                    D_thick = (ii_target - 1) * target + tt
                    sigma_beam_straggle_abs = sigma_beam_straggle * M1_stop * D_thick
                    E0 = gauss(E0, sigma_beam_straggle_abs)
                ENDIF

!    select one of targets randomly and randomise beam energy because of straggling
        
        
            ELSE
        
                z_vertexcm = RAND() * z_interact_startcm
                IF (plot.EQ.1 .AND. thetacm_step.EQ.180.) THEN 
                 z_vertexcm = 0.5 * z_interact_startcm
                ENDIF  
            
!    for active target reaction occurs beween 0 and z_vertex
            
                gas_lengthcm = z_vertexcm - window_positioncm
 
!    beam traverses this distance to interaction point
                       

 
                D_thick = gas_lengthcm * density * 1000. + window_thickness
 
!     effective gas thickness for beam stopping, in mg/cm^2
!   add nominal amount to take into account gas window foil
            
            
                IF (beam_stop.EQ.1) THEN 
                 E0 = E0 - M1_stop * D_thick
                ENDIF

!    beam loses energy in the target
!    for active target should not correct for beam energy loss as this can be determined.


                sigma_beam_straggle_abs = sigma_beam_straggle * M1_stop * D_thick
                
!     WRITE #2, M1_stop, gas_lengthcm, D_thick, sigma_beam_straggle_abs

                E0 = gauss(E0, sigma_beam_straggle_abs)

!    for active target estimate beam energy straggling.
        

            ENDIF
        

            ET = E0 + Qvalue
            factor1 = (1. + M1 / M2 * Qvalue / ET)
            factor2 = (M1 + M2) * (M3 + M4)
            A = M1 * M4 / factor2 * E0 / ET
            B = M1 * M3 / factor2 * E0 / ET
            C = M2 * M3 / factor2 * factor1
            D = M2 * M4 / factor2 * factor1
    

            E4 = (A + C + 2 * SQRT(A * C) * COS((180. - thetacm) * pi / 180.)) * ET
            E3 = (B + D + 2 * SQRT(A * C) * COS(thetacm * pi / 180.)) * ET
            sinpsi = SIN(thetacm * pi / 180.) / SQRT(E3 / (ET * D))
            psi = ASIN(sinpsi) * 180. / pi
            factor3 = E3 * (1. - M3 / M1) + E4 * (1. - M4 / M1) - Qvalue
            cos_zeta_plus_psi = M1 / (2 * SQRT(E3 * E4 * M3 * M4)) * factor3
            zeta0 = ACOS(cos_zeta_plus_psi) * 180. / pi - psi
            
            IF (z1detcm.LT.0 .AND. zeta0.LT.91) THEN
                loops = loops - 1
                WRITE(2,129)
129             FORMAT("zeta approaching 90 degrees, terminate calculation")
                GOTO 11
            ENDIF
            

!    2-body kinematics
            
            E4_start = E4
            
!   initial energy of ejectile


            zeta1 = gauss(zeta0, sigma_angle_deg)

!   randomise ejectile angle because of beam divergence


            zeta = zeta1
            
            
            rho_max_0 = rho_max_factor * SQRT(E4 * M4) * SIN(zeta * pi / 180.) / (charge_state * Bz)
            IF (rho_max_0 .GT. Rmagnet + 0.01) THEN
                loops = loops - 1
                WRITE(2,130)
130             FORMAT("rho_max too large, terminate calculation")
                GOTO 11
            ENDIF
            
!   if trajectory radius of ejectile exceeds detector radius by 1cm discontinue theta loop
            
            sigma_scat_beam_deg = 0.
            sigma_scat_beam_foil_deg = 0.

!  bug corrected 21/7/24 - these variables not initialised if scattering flags set to zero

            IF (ii_scat_max .GT. 0 .AND. target .GT. 0.) THEN
                sigma_scat_beam_deg = multi_scatter(Z_target, A_target, charge_state_beam, E0, D_thick) /    &
                      scatter_factor * (180. / pi)
                
  
  
                IF (active_flag .EQ. 1) THEN
                    sigma_scat_beam_foil_deg = multi_scatter(6.D0, 12.D0, charge_state_beam, E0, 0.05D0) /    &
                         scatter_factor * (180. / pi)
                ENDIF


!   for active target take into account beam scattering in window foil

                sigma_scat_beam_deg_tot = SQRT(sigma_scat_beam_deg ** 2 + sigma_scat_beam_foil_deg ** 2)
                zeta = gauss(zeta1, sigma_scat_beam_deg_tot)
            ENDIF
            
                                    

!    randomise ejectile angle because of multiple scattering of beam

            IF (E4 .LT. 0.) THEN 
             E4 = 0.
            ENDIF
            IF (E4 .LT. E4_thresh_MeV) THEN
                error_E4low = error_E4low + 1
                GOTO 10
            ENDIF
            
!     check that light particle energy is above threshold


            M4_stop = terpolate(E4)
            
            IF (M4_stop .EQ. 0.) THEN
                error_dedx = error_dedx + 1
                GOTO 10
            ENDIF
 
!     stopping power of ejectile


            IF (active_flag .EQ. 0) THEN


                IF (target .GT. 0.) THEN
                    IF (zeta .GT. 90.) THEN
                        D_thick = -tt / COS(zeta * pi / 180.)
                    ELSE
                        D_thick = (target - tt) / COS(zeta * pi / 180.)
                    ENDIF
                    IF (ii_eject_scat .EQ. 1) THEN
                        sigma_scat_deg = multi_scatter(Z_target, A_target, charge_state, E4, D_thick) / scatter_factor * (180. / pi)
                        zeta = gauss(zeta, sigma_scat_deg)
                    ENDIF

!         multiple scattering of ejectile in the target



                    sigma_eject_straggle_abs = sigma_eject_straggle * M4_stop * D_thick
                    E4 = gauss(E4, sigma_eject_straggle_abs)

!     straggling of ejectile in target


                ENDIF


                IF (eject_stop .EQ. 1 .AND. target .GT. 0) THEN
                    IF (zeta .GT. 90.) THEN
                        E4 = E4 + M4_stop * tt / COS(zeta * pi / 180.)
                    ELSE
                        E4 = E4 - M4_stop * (target - tt) / COS(zeta * pi / 180.)
                    ENDIF
                ENDIF

!     energy loss of ejectile in the target


            ENDIF
            IF (E4 .LT. 0) THEN 
             E4 = 0.
            ENDIF
            IF (E4 .LT. E4_thresh_MeV) THEN
                error_E4low = error_E4low + 1
                GOTO 10
            ENDIF

!    check that light particle energy is above threshold after losing energy in target
        
        
            CONTINUE
            sigma_beam_spot_tot = sigma_beam_spot
 
            IF (active_flag .EQ. 1) THEN
                scatter_displace = sigma_scat_beam_deg * pi / 180. * gas_lengthcm / 200.
                scatter_displace_foil = sigma_scat_beam_foil_deg * pi / 180. * gas_lengthcm / 100.
                sigma_beam_spot_tot = SQRT(sigma_beam_spot ** 2 + scatter_displace ** 2 + scatter_displace_foil ** 2)
            END IF
            
            
!     for active target take into account displacement in x,y due to beam scattering.
!     For multiple scattering of beam in gas, reduce gas length by factor of 2 as this incrementally increases
            
            x = gauss(0.D0, sigma_beam_spot_tot)
            y = gauss(0.D0, sigma_beam_spot_tot)
            radius = SQRT(x * x + y * y)

!    IF radius > 2. * sigma_beam_spot_tot THEN GOTO 5
 
!    randomise starting x,y because of finite beam spot
            
            x_start = x * 1000.
            y_start = y * 1000.

            IF (active_flag .EQ. 0) THEN
                z = 0. 
            ELSE       
                z_sigma_pad = FWHM_z_padmm / 2350.
                z = z_vertexcm / 100           
            ENDIF
        


            v = SQRT(2. * E4 * 1.60218E-13 / mass)

!   velocity of ejectile+


            IF (plot .EQ. 0) THEN
                phi = 2 * pi * RAND()
            ELSE
!     phi = pi#
                phi = 5. * pi / 4.
!   (to compare with 46Ar(p,p') in Bradt et al.)
            ENDIF

!    random azimuthal angle


            vz = v * COS(zeta * pi / 180.)
            vx = v * SIN(zeta * pi / 180.) * COS(phi)
            vy = v * SIN(zeta * pi / 180.) * SIN(phi)
        

        
            check_pad = 1
            n = 0
            nprev = 0
            t = 0
            s_travel = 0
            s_travel_tot = 0
        
            x0 = x
            y0 = y
            z0 = z
            sum_d_theta = 0.
            sum_d_phi = 0.


            gas_thick_tot = 0.
            Delta_x_tot = 0.
            Delta_y_tot = 0
            Delta_z_tot = 0
            Delta_s_tot = 0

                
            DO WHILE (t .LT. Tcyc * 0.5 .OR. radius .GT. Rdet)

!    loop until distance from beam axis < radius of detector. The first condition ensures that particle is returning to the beam axis


                IF (((z1det .LT. 0. .AND. z .LT. z1det) .OR. (z1det .GT. 0 .AND. z .GT. z1det)) .AND. radius .LT. Rdet) THEN
                    error_traj_low = error_traj_low + 1
                    GOTO 10
                ENDIF

!    check that trajectory is not too low for finite size detector


                IF ((z2det .GT. 0 .AND. z .GT. z2det) .OR. (z2det .LT. 0 .AND. z .LT. z2det)) THEN
                    error_det_overrun = error_det_overrun + 1
                    GOTO 10
                ENDIF

!    check that have not overran detector


                IF (t .GT. Tcyc * 2.) THEN
                    error_slow = error_slow + 1
                    GOTO 10
                ENDIF

!    check that trajectory time does not excessively exceed cyclotron time


                IF (radius .GT. Rmagnet) THEN
                    error_large_radius = error_large_radius + 1
                    GOTO 10
                ENDIF

!     check that radius of trajectory is within magnet radius
                
                
                z_interact_start = z_interact_startcm / 100.
                IF (radius .GT. Rpad .AND. check_pad .EQ. 1 .AND. active_flag .EQ. 1 .AND. ((z1detcm .LT. 0 .AND.   &
                        z .GT. (z1detcm + 5.)) .OR. (z1detcm .GT. 0. .AND. z .LT. (z1detcm - 5.)))) THEN

                
                    check_pad = 0
                    z_meas_pad = gauss(z, z_sigma_pad)
                    
!     randomise z measurement at pad because of uncertainty of its determination in active target
                
                    radius_i = INT(radius * 1.0E5)
                    kradius = INT(radius_i / pitch_radius_pad)
                    R_meas_pad = (kradius * pitch_radius_pad) / 1.E5 + 0.5 * pitch_rpadmm / 1000.
                    
!     measured radius of pad is binned because of sensor pitch
                                                                  
                ENDIF

!     as ejectile trajectory crosses the inner radius of the pad detector, record values of z and R
                                                                                                            
        
                s_iteration = v * deltat

!    distance travelled each iteration in m
            

                s_travel = s_travel + s_iteration
                s_travel_tot = s_travel_tot + s_iteration
                
!     total distance travelled by ejectile
            
                r_from_target = SQRT(x * x + y * y + z * z)
                Bzz = Bz - (((r_from_target / 0.5) ** 2) * inhomo_Bz)
                Bx = x * z * inhomo_Bz / (0.5 ** 2)
                By = y * z * inhomo_Bz / (0.5 ** 2)

!    include effect of field inhomogeneity
                

            
                      
                IF (active_flag .EQ. 1 .AND. s_travel .GT. s_travel_increment) THEN
            
                    gas_thick = s_travel * 100. * density * 1000.

!    thickness of gas travelled through in ~ 1cm increments (~ 1mm for isobutane), in mg/cm**2
                

                    gas_thick_tot = gas_thick_tot + gas_thick
            
                    E4_0 = E4
                    
                    M4_stop = terpolate(E4)
                    
                    IF (M4_stop .EQ. 0) THEN
                        error_dedx = error_dedx + 1
                        GOTO 10
                    ENDIF
            
                    IF (eject_stop .EQ. 1) THEN
                        E4 = E4 - M4_stop * gas_thick
                    ENDIF
                
!     ejectile loses energy in the gas
!     if active target should not correct for ejectile energy loss as this can be determined.
                    
                    IF (E4 .LT. 0) THEN 
                     E4 = 0.
                    ENDIF
                    IF (E4 .LT. E4_thresh_MeV) THEN
                        error_E4low = error_E4low + 1
                        GOTO 10
                    ENDIF
                    
!     check that light particle energy is above threshold after losing energy in gas
                    
                    
                    sigma_eject_straggle_abs = sigma_eject_straggle * M4_stop * gas_thick
                    E4 = gauss(E4, sigma_eject_straggle_abs)
!    straggling of ejectile in gas

                    IF (E4 .LT. 0) THEN 
                     E4 = 0.
                    ENDIF
                    IF (E4 .LT. E4_thresh_MeV) THEN
                        error_E4low = error_E4low + 1
                        GOTO 10
                    ENDIF

!     check that light particle energy is above threshold after straggling

            

                
            
            
                    factor = SQRT(E4 / E4_0)
                    vx = factor * vx
                    vy = factor * vy
                    vz = factor * vz
                

                    IF (ii_eject_scat .EQ. 1) THEN
            
                        r_travel = SQRT((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)
 
                        cos_theta_z = (z - z0) / r_travel
                        sin_theta_z = SQRT(1 - cos_theta_z * cos_theta_z)
                        sin_phi_z = (y - y0) / (r_travel * sin_theta_z)
                        cos_phi_z = (x - x0) / (r_travel * sin_theta_z)
                    
                        sigma_scat_rad = multi_scatter(Z_target, A_target, charge_state, E4, gas_thick) / scatter_factor
                    
                    
                        d_theta = gauss(0.D0, sigma_scat_rad)
                                                                                          
                        sum_d_theta = sum_d_theta + d_theta
                        
                        Delta_x = (z - z0) * cos_phi_z * sum_d_theta
                        Delta_y = (z - z0) * sin_phi_z * sum_d_theta
                        Delta_z = -r_travel * sin_theta_z * sum_d_theta
                        
                    
                        x = x + Delta_x
                        y = y + Delta_y
                        z = z + Delta_z
                    
                        
                   
!   multiple scattering of ejectile in the gas. The angular change has to be a running sum over the whole path
                    
                    
                        
                        Delta_x_tot = Delta_x_tot + Delta_x * 1000.
                        Delta_y_tot = Delta_y_tot + Delta_y * 1000.
                        Delta_z_tot = Delta_z_tot + Delta_z * 1000.
                        
                    
!     tests
 
                    ENDIF
                    
                    s_travel = 0.
            
                    x0 = x
                    y0 = y
                    z0 = z
                        
                ENDIF

                CONTINUE
                       
                vx = -q * vy * Bzz * deltat / mass + q * vz * By * deltat / mass + vx
                vy = -q * vz * Bx * deltat / mass + q * vx * Bzz * deltat / mass + vy
                vz = -q * vx * By * deltat / mass + q * vy * Bx * deltat / mass + vz
                
                IF (active_flag .EQ. 1) THEN
                    vz = vz - q * Ez * deltat / mass
                ENDIF

                IF (n .EQ. nprev .AND. plot .EQ. 1) THEN

!    plotting variables

                    x1 = x * xscale + 150.
                    y1 = y * yscale + 250.
                    z1 = z * zscale + z_plot_offset
                ENDIF

                x = vx * deltat + x
                y = vy * deltat + y
                z = vz * deltat + z
                radius = SQRT(x * x + y * y)
                
                IF (radius .GT. radius_max) THEN 
                 radius_max = radius
                ENDIF
!    distance from beam axis
  

                n = n + 1
                t = t + deltat

!           total time


                IF (n .EQ. nprev + 10 .AND. plot .EQ. 1) THEN
 
!   plot x-y and z-y every 100 increments
  
                    x2 = x * xscale + 150.
                    y2 = y * yscale + 250.
                    z2 = z * zscale + z_plot_offset
!                    LINE (x1, y1)-(x2, y2), 4
!                    plots in x, y plane, colour red on white background
!                    LINE (z1, y1)-(z2, y2), 4
!                    plots in z, y plane, colour red on white background
!                     plot trajectory
!                     can remove the x,y or z,y plot for clarity
                    x_plot = x * 1000.
                    y_plot = y * 1000.
                    z_plot = z * 1000.
                    E4_out = E4
!                    WRITE #5, E4_out; "   "; x_plot; "   "; y_plot; "   "; z_plot
                    nprev = n
                ENDIF

            END DO
 
 
            IF ((z1det .LT. 0. .AND. z .GT. z1det) .OR. (z1det .LT. 0 .AND. z .LT. z2det) .OR. (z1det .GT. 0. .AND.  &
                           z .LT. z1det) .OR. (z1det .GT. 0. .AND. z .GT. z2det)) THEN
                error_miss_detector = error_miss_detector + 1
                GOTO 10
            ENDIF
 
 
!   check that intersection is at location of detector
            
            
            IF (active_flag .EQ. 1 .AND. check_pad .EQ. 1) THEN
                error_no_zvertex = error_no_zvertex + 1
                GOTO 10
            ENDIF
            
!    check that there has been a valid estimate made of the z position of the reaction vertex
 

            vend = SQRT(vx * vx + vy * vy + vz * vz)

!    check that velocity has not changed

            E4_meas = gauss(E4, sigma_E4)

!    randomise measured energy because of detector resolution
            
            IF (E4 .LT. 0) THEN 
             E4 = 0.
            ENDIF
            IF (E4_meas .LT. E4_thresh_MeV) THEN
                error_E4low = error_E4low + 1
                GOTO 10
            ENDIF

!   check that light particle energy is above threshold after randomising energy
            

            
            zi = INT(z * 1.0E5)
            kz = INT(zi / pitch_z)
            z_meas = (kz * pitch_z) / 1.E5 + 0.5 * pitch_zmm / 1000.

            xi = INT(x * 1.0E5)
            kx = INT(xi / pitch_x)
            x_meas = (kx * pitch_x) / 1.E5 + 0.5 * pitch_xmm / 1000.

!              bin measured z and x using Si array because of detector pitch
            
            
            Q_corr = 0
            z_corr = 0
                        
            iter_max = 10
            IF (active_flag .EQ. 0) THEN
                iter_max = 2
                z_vertex = 0.
                Delta_z_pad = 0.
                z_meas_pad = 0.
            ENDIF
                                    
            DO 200 iter = 1, iter_max
            
                z_corr_prev = z_corr
                
                IF (active_flag .EQ. 1) THEN
                    
                    IF (iter .LT. 6) THEN 
                     z_corr = ((AA * E4_meas + BB - Q_corr) / CC) / 100.
                    ENDIF
                    
!              estimate of total distance between vertex and intersection of ejectile trajectory with beam axis.  Starting Q value is nominal

                    zeta_meas = ACOS(z_corr * charge_state * Bz / (z_factor * SQRT(E4_meas * M4)))
                    rho_max_meas = rho_max_factor * SQRT(E4_meas * M4) * SIN(zeta_meas) / (charge_state * Bz)
                    alpha = 2. * ASIN(R_meas_pad / rho_max_meas)
                    Delta_z_pad = -z_corr * alpha / (2. * pi)
                    z_vertex = z_meas_pad + Delta_z_pad
                    
!                     WRITE #2, "iter "; iter%; "z_vertex "; z_vertex#
                    
!                    this extrapolates the value of z at the pad inner radius to the vertex position, using PAB algorithm.
                                
                ENDIF
                                                                                                                    
                IF (active_flag .EQ. 0) THEN 
                 z_corr = z_meas
                ENDIF

                R_meas = SQRT(x_meas * x_meas + y * y)


                DO 201 jj = 1, 5
                    zeta_meas = ACOS(z_corr * charge_state * Bz / (z_factor * SQRT(E4_meas * M4)))
                    rho_max_meas = rho_max_factor * SQRT(E4_meas * M4) * SIN(zeta_meas) / (charge_state * Bz)
                    alpha = 2. * ASIN(R_meas / rho_max_meas)
                    Delta_z_Si = z_corr * alpha / (2. * pi)
                    z_corr = z_meas + Delta_z_Si - z_vertex
                    
!                    WRITE #2, "iter "; iter%; "z_corr "; z_corr#
201             CONTINUE
!              this corrects the measured z at the Si array to the true z, using PAB algorithm
                
            
                Q_corr = AA * E4_meas + BB - CC * z_corr * 100.
 
                IF (ABS(z_corr - z_corr_prev) .LT. 0.0001) THEN 
                  GOTO 155
                ENDIF
!           check for convergence
            
200         CONTINUE
            
            
            error_no_converge = error_no_converge + 1
            GOTO 10
!       reject if does not converge after iter_max% iterations
            
155         CONTINUE
            
            z_change_padcm = Delta_z_pad * 100.
            z_changecm = Delta_z_Si * 100.
        
            rho_max_cm = rho_max_meas * 100.
        
            Qvalue_uncorr = AA * E4_meas + BB - CC * (z_meas - z_meas_pad) * 100.
            Qvalue_corr = Q_corr
            
            
            IF (ABS(Qvalue_uncorr - Qvalue_corr) .GT. 2.) THEN

!         original criteria was that error flag = 1 if the absolute value of the difference > 10 MeV (August 2022)
!         reduced to > 1 MeV later, but seems to cut off forward angles for (alpha,p) reaction.  Suggest (4th Jan 2023) 2 MeV as a compromise.

                error_poor_Q = error_poor_Q + 1
                GOTO 10
            ENDIF
            
            
!      reject if calculate Q value that is clearly in error
            
            z_vertexcm = z_vertex * 100.
            
            phi_deg = phi * 180. / pi
                       
            tns = t * 1E9
            tns_corr = tns * 2. * pi / (2. * pi - alpha)
        

            IF (ABS(Tcyc_ns - tns_corr) .GT. 1) THEN
                error_slow = error_slow + 1
                GOTO 10
            ENDIF

!               reject slow events at edge of detector
        

            ET_meas = E0beam + Qvalue_corr
            cosangle = (E4_meas / ET_meas - A0 - C0) / (2 * SQRT(A0 * C0))
            thetacm_meas = 180. - ACOS(cosangle) * 180. / pi
            
            

            zcm = (z_meas - z_meas_pad) * 100.
            z_corrcm = z_corr * 100.
            
            good_loops = good_loops + 1
            good_pass = good_pass + 1

            Q_store_pass(good_pass) = Qvalue_corr
            Q_store_loops(good_loops) = Qvalue_corr
            v_diff_sq_sum = (vend - v) * (vend - v) + v_diff_sq_sum
            v_sum = v_sum + v
            
            
            E4_sum = E4_sum + E4
            zeta_sum = zeta_sum + zeta
            phi_deg_sum = phi_deg_sum + phi_deg
            tns_sum = tns_sum + tns
            z_vertexcm_sum = z_vertexcm_sum + z_vertexcm
            zcm_sum = zcm_sum + zcm
            Qvalue_uncorr_sum = Qvalue_uncorr_sum + Qvalue_uncorr
            z_corrcm_sum = z_corrcm_sum + z_corrcm
            z_change_padcm_sum = z_change_padcm_sum + z_change_padcm
            z_changecm_sum = z_changecm_sum + z_changecm
            Qvalue_corr_sum = Qvalue_corr_sum + Qvalue_corr
            tns_corr_sum = tns_corr_sum + tns_corr
            rho_max_cm_sum = rho_max_cm_sum + rho_max_cm
            s_travel_tot_sum = s_travel_tot_sum + s_travel_tot
            E4_loss_sum = E4_loss_sum + (E4_start - E4)
            thetacm_meas_sum = thetacm_meas_sum + thetacm_meas
            
!            
!             IF (thetacm .GT. 10.99 AND thetacm < 11.01) THEN
!                 WRITE #2, Q_store_pass(good_pass&); "    "; E4#; "    "; E4_meas#; "    "; z_corr#
!              END IF
!               For small theta_cm the large z correction is strongly correlated to the measured energy in the Si array, resulting in a smaller overall spread in Q values.
!           

10          CONTINUE

        END DO
        
        IF (good_pass .GT. 0) THEN
 
            efficiency = FLOAT(good_pass) / FLOAT(n_pass) * 100.
            E4_av = E4_sum / FLOAT(good_pass)
            zeta_av = zeta_sum / FLOAT(good_pass)
            phi_deg_av = phi_deg_sum / FLOAT(good_pass)
            tns_av = tns_sum / FLOAT(good_pass)
            z_vertexcm_av = z_vertexcm_sum / FLOAT(good_pass)
            zcm_av = zcm_sum / FLOAT(good_pass)
            Qvalue_uncorr_av = Qvalue_uncorr_sum / FLOAT(good_pass)
            z_corrcm_av = z_corrcm_sum / FLOAT(good_pass)
            z_change_padcm_av = z_change_padcm_sum / FLOAT(good_pass)
            z_changecm_av = z_changecm_sum / FLOAT(good_pass)
            Qvalue_corr_av = Qvalue_corr_sum / FLOAT(good_pass)
            tns_corr_av = tns_corr_sum / FLOAT(good_pass)
            rho_max_cm_av = rho_max_cm_sum / FLOAT(good_pass)
            s_travelcm_tot_av = s_travel_tot_sum / FLOAT(good_pass) * 100.
            E4_loss_av = E4_loss_sum / FLOAT(good_pass)
            thetacm_meas_av = thetacm_meas_sum / FLOAT(good_pass)
        
            Q_pass_sum = 0.
 
    
            DO 202 n = 1, good_pass
                Q_pass_sum = Q_pass_sum + Q_store_pass(n)
202         CONTINUE


            Q_pass_average = Q_pass_sum / good_pass
        
            sigma_pass_sq = 0.
            DO 203 n = 1, good_pass
                sigma_pass_sq = sigma_pass_sq + (Q_store_pass(n) - Q_pass_average) * (Q_store_pass(n) - Q_pass_average)
203         CONTINUE

            sigma_pass = SQRT(sigma_pass_sq / FLOAT(good_pass)) * 1000.
            fwhm_pass = 2.35 * sigma_pass
            
            radius_maxcm = radius_max * 100.

!    tmp$ = "###,.##  ###,.## ###,.## ###,.#  ##,.##  ##,.##   ##,.##  ##,.##  ##,.##     ##,.##  ##,.##  ##,.###   ###,.# ###,.## ###,.## ###,.## ##,.### ###,.##   ###"           
        
            WRITE(2,131) thetacm, E4_av,  zeta_av, tns_av, zcm_av,  Qvalue_uncorr_av, z_change_padcm, z_changecm_av,  &
                   z_vertexcm_av, tns_corr_av, z_corrcm_av, Qvalue_corr_av, fwhm_pass, rho_max_cm_av, radius_maxcm,        &
                   s_travelcm_tot_av, E4_loss_av, thetacm_meas_av, efficiency
131         FORMAT(      F5.2,6X, F5.2,2X,F7.2,4X, F4.1,1X F7.2,3X, F5.2,3X,          F5.2,3X,        F6.2,4X,       F6.2,2X,      &  
                     F6.2,2X,     F6.2,2X,     F6.2,3X,        F6.1,3X,   F5.2,3X,       F5.2,3X,        &   
                     F6.2,2X,           F6.3,3X,    F5.2,2X,         F5.0)

        ELSE
        
            IF (Rdetcm .LT. 1. .OR. Rpadcm .LT. 1.) THEN
                radius_maxcm = radius_max * 100.
                WRITE(2,132) thetacm, radius_maxcm
132             FORMAT("theta c. of m. =", F3.1,"degrees    maximum radius =",F5.2,"cm")
            ENDIF
            
        ENDIF
        
        thetacm = thetacm + thetacm_step
       END DO
11     CONTINUE

       IF (good_loops .EQ. 0) THEN
        WRITE(2,141)
141     FORMAT("no good events")
        GOTO 601
       ENDIF

       v_diff_rms_percent = 100. * SQRT(v_diff_sq_sum) / v_sum * SQRT(FLOAT(good_loops))

       Q_sum = 0.
    
       
       DO 204 n = 1, good_loops
        Q_sum = Q_sum + Q_store_loops(n)
204    CONTINUE


       Q_average = Q_sum / good_loops
    

       WRITE(2,133)
133    FORMAT(" ")
       sigma_sq = 0.
       DO 205 n = 1, good_loops
        sigma_sq = sigma_sq + (Q_store_loops(n) - Q_average) * (Q_store_loops(n) - Q_average)
205    CONTINUE
    
    

       WRITE(2,133)
    
       sigma = SQRT(sigma_sq / FLOAT(good_loops)) * 1000.
       fwhm = 2.35 * sigma
    

    


       WRITE(2,134)v_diff_rms_percent
134    FORMAT("average rms velocity error = ",F6.4," %")

       WRITE(2,135)Q_average, sigma, fwhm
135    FORMAT("average measured Q value = ",F6.3, " MeV    standard deviation = ",F4.1," keV  FWHM = ",F5.1," keV")

       WRITE(2,133)


601    CONTINUE
    

       WRITE(2,136) loops, good_loops
136    FORMAT("total number of trajectories", I6, "  successful trajectories", I6)
       WRITE(2,137)error_traj_low, error_det_overrun, error_slow
137    FORMAT("low trajectory",I6, "  detector overrun",I6,"  too slow",I6,$)
       WRITE(2,138)error_large_radius, error_dedx
138    FORMAT("  too large radius",I6, "  de/dx lookup error", I6,$)
       WRITE(2,139)error_E4low, error_miss_detector
139    FORMAT("  E4 low", I6, "  miss detector", I6,$)
       WRITE(2,140)error_no_zvertex, error_no_converge, error_poor_Q
140    FORMAT("  no z vertex", I6, "  extrapolation fail", I6, "  poor Q value", I6)
    
       WRITE(2,133)
       WRITE(2,133)
       WRITE(2,133)
    

      END DO
	  
302   CONTINUE
      IF (plot .EQ. 0) THEN 
        WRITE(*,353) 
353     FORMAT(" ")
       GOTO 1001
      ENDIF

1001  END PROGRAM


      REAL*8 FUNCTION gauss (AM, S)
      IMPLICIT REAL*8(A-H,O-Z)
      A = 0.0
      DO 1 I = 1, 12
      A = A + RAND()
1     CONTINUE
      gauss = (A - 6.0) * S + AM
      RETURN
      END

      REAL*8 FUNCTION terpolate (x)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/stopping/energy4(1000), dedx4(1000), pi, ii_max

      IST = 0

      DO 1 I = 1, ii_max - 1
      IF (energy4(I) .LE. x .AND. energy4(I + 1) .GE. x) THEN
        IST = I
        GO TO 2
      END IF
1     CONTINUE

2     CONTINUE

      IF (IST .EQ. 0) THEN
       terpolate = 0.
       GO TO 120
      ENDIF

      ZM = LOG(dedx4(IST + 1) / dedx4(IST)) / LOG(energy4(IST + 1) / energy4(IST))
      C = LOG(dedx4(IST)) - ZM * LOG(energy4(IST))
      Y = ZM * LOG(x) + C
      terpolate = EXP(Y)
120   CONTINUE
      RETURN
      END

      REAL*8 FUNCTION multi_scatter (Z, A, Q, T, D)
      IMPLICIT REAL*8(A-H,O-Z)
      XK = .0393 / 1000. * Z * (Z + 1.) * Q ** 2 * D / A / T ** 2
      B = 2.119 * Z ** (-5. / 3.) * (Z + 1.) * D / A
      XLB = LOG(B) + 11.51
      IF (XLB .GT. 0) THEN
        F = SQRT(XK * XLB)
      ELSE
        F = 0.
      ENDIF
!      multi_scatter = 1.36 * F
!  original statement
      multi_scatter = F
!   RMS scattering angle in radians, taken from Kantele's handbook
      RETURN
      END



