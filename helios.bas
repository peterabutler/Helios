REM HELIOS simulations.  PA Butler  April 2014.  Modified for active gas volume mode June 2022.

DIM Q_store_pass(1 TO 100000), Q_store_loops(1 TO 1000000)
COMMON SHARED pi#, ii_max%
DIM SHARED energy4(1 TO 1000), dedx4(1 TO 1000)

SCREEN 9

xframe% = 640
yframe% = 350
xscale = xframe%
yscale = xscale
zscale = xscale
plot% = 0
file_plot_open% = 0
REM   plot trajectory

sintheta_weighting% = 0

RANDOMIZE TIMER

OPEN "reaction.txt" FOR INPUT AS #4
LINE INPUT #4, folder_name$
LINE INPUT #4, reaction_label$

file_in$ = folder_name$ + "\input_" + reaction_label$ + ".txt"
OPEN file_in$ FOR INPUT AS #1

file_out$ = folder_name$ + "\output_" + reaction_label$ + ".txt"
OPEN file_out$ FOR OUTPUT AS #2

Ez = 2500 / 0.15
REM electric field if active target (V/m)

DO WHILE NOT EOF(1)

    INPUT #1, dummy$
    INPUT #1, M1, M2, M4, EoA, Qvalue

    INPUT #1, dummy$
    INPUT #1, Bz, charge_state, Rmagnetcm, inhomo_Bz_percent

    INPUT #1, dummy$
    INPUT #1, Rdetcm, z1detcm, z2detcm

    INPUT #1, dummy$
    INPUT #1, pitch_zmm, pitch_xmm, fwhm_E4keV, E4_thresh_MeV

    INPUT #1, dummy$
    INPUT #1, target, M1_stop, M4_stop_file$, beam_stop%, eject_stop%

    INPUT #1, dummy$
    INPUT #1, ii_scat_max%, Z_beam, ii_eject_scat%, fwhm_beam_straggle_percent, fwhm_eject_straggle_percent
    
    INPUT #1, dummy$
    INPUT #1, fwhm_beam_spot_mm, fwhm_E0_percent, fwhm_diverg_mrad

    INPUT #1, dummy$
    INPUT #1, thetacm_step

    IF thetacm_step = 0 THEN
        plot% = 1
        thetacm = 35
        thetacm_step = 180
        COLOR , 63
    END IF
    
    IF thetacm_step < 0 THEN
        plot% = 1
        thetacm = -thetacm_step
        thetacm_step = 180
        COLOR , 63
    END IF
    
    IF plot% = 1 AND file_plot_open% = 0 THEN
        file_plot$ = folder_name$ + "\plot_" + reaction_label$ + ".txt"
        OPEN file_plot$ FOR OUTPUT AS #5
        file_plot_open% = 1
    END IF
    IF plot% = 1 THEN PRINT #5, "E4            x             y             z"
    
    IF z1detcm > 0 THEN z_plot_offset = 300
    IF z1detcm < 0 THEN z_plot_offset = 600
    
    REM   plot trajectory if plot% = 1
    
    

    M4_stop_filename$ = "srim_files\" + M4_stop_file$

    OPEN M4_stop_filename$ FOR INPUT AS #3

    Bx = 0
    By = 0
    Z_target = 0
    A_target = 0
    
    IF M4_stop_file$ = "proton_srim.txt" OR M4_stop_file$ = "deuterium_srim.txt" OR M4_stop_file$ = "carbon_srim.txt" OR M4_stop_file$ = "proton_on_isobutane_srim.txt" THEN
        Z_target = 6.
        A_target = 12.
    END IF
    
    IF M4_stop_file$ = "proton_on_CD2_srim.txt" OR M4_stop_file$ = "deuteron_on_CD2_srim.txt" OR M4_stop_file$ = "carbon_on_carbon_srim.txt" THEN
        Z_target = 6.
        A_target = 12.
    END IF
    
    IF M4_stop_file$ = "proton_on_deuterium_srim.txt" THEN
        Z_target = 1.
        A_target = 2.
    END IF
    
    IF M4_stop_file$ = "deuteron_on_deuterium_srim.txt" OR M4_stop_file$ = "He-3_on_deuterium_srim.txt" THEN
        Z_target = 1.
        A_target = 2.
    END IF
    
    IF M4_stop_file$ = "alpha_on_helium_srim.txt" OR M4_stop_file$ = "alpha_on_He-CO2_srim.txt" THEN
        Z_target = 2.
        A_target = 4.
    END IF
    
    IF M4_stop_file$ = "proton_on_helium_srim.txt" OR M4_stop_file$ = "proton_on_He-CO2_srim.txt" THEN
        Z_target = 2.
        A_target = 4.
    END IF
    
    IF Z_target = 0 THEN
        PRINT #2, "target not specified"
        STOP
    END IF
    
    scatter_factor_active = 1.8
    
    active_flag% = 0
    IF target > 1 THEN active_flag% = 1

    REM if value of target thickness is greater than 1, then this value is assumed to be the gas pressure (in torr) of an active target
        
    density_760 = 0.
    REM this statement has been inserted February 2025

    IF active_flag% = 1 THEN
        s_travel_increment = 0.01
        density_760 = 0.
        
        
        IF ii_eject_scat% < 0 THEN
            scatter_factor_active = 1.0
            ii_eject_scat% = 1
        END IF
        REM if scatter flag < 0 do not modify multiple scattering for gas target
         
          
    
        IF Z_target = 1. AND A_target = 2. THEN
            density_760 = 0.1702 / 1000
        END IF
        REM deuterium
    
        IF Z_target = 2. AND A_target = 4. THEN
            density_760 = 0.1662 / 1000
        END IF
        REM helium
    
        IF Z_target = 6. AND A_target = 12. THEN
            density_760 = 2.53 / 1000
            s_travel_increment = 0.001
        END IF
        REM isobutane
    
        IF density_760 = 0. THEN
            PRINT #2, "gas density not specified"
            STOP
        END IF
    
    END IF

    REM gas density in g/cm^3 at 15C, atmospheric pressure

    

    density = density_760 * target / 760

    REM actual gas density in g/cm^3
   
 
       
    ii% = 0
    DO WHILE NOT EOF(3)
        ii% = ii% + 1
        INPUT #3, energy4(ii%), dedx4(ii%)
    LOOP
    ii_max% = ii%

    CLOSE #3

    sigma_beam_straggle = fwhm_beam_straggle_percent / 235.
    sigma_beam_spot = fwhm_beam_spot_mm / 2350.
    sigma_eject_straggle = fwhm_eject_straggle_percent / 235.




    Rmagnet = Rmagnetcm / 100

    REM radius of magnet
    
    scatter_factor = 1.
	
    z_interact_startcm = 0.
    REM    this line inserted February 2025 

    IF active_flag% = 1 THEN
    
        scatter_factor = scatter_factor_active
    
        REM this emperically adjusts the value of sigma for multiple scattering so that the Kantele subroutine
        REM agrees with the Fano prescription, which reproduces the data of Kuhn et al. NIM B4 (1984) 332
    
        z_interact_startcm = z2detcm
     
        IF (z1detcm < 0) THEN
            z2detcm = z1detcm - 50.
            window_positioncm = z1detcm - 5.
        ELSE
            z2detcm = z1detcm + 50.
            window_positioncm = -5.
        END IF
        REM window_positioncm is position of gas window, in cm
    END IF
    
    REM for active target the Si detector is assumed to be 50cm long
    REM in this case the interaction region is between 0 and z_interact_startcm
    
    window_thickness = 0.1
    REM for active target this is window foil thickness in mg/cm**2
    

   
    
    Rpadcm = 5.5
    cent% = 0
	
	pitch_rpadmm = 0.
    REM this line inserted February 2025

    IF cent% = 0 AND active_flag% = 1 THEN
        Rpadcm = Rdetcm
        Rdetcm = 3.
    END IF
    
    Rdet = Rdetcm / 100

    REM radius of Si detector (m)
    
    Rpad = Rpadcm / 100

    REM radius of inner pads (m)
        

    IF active_flag% = 1 THEN
        FWHM_z_padmm = pitch_zmm
        pitch_rpadmm = pitch_xmm
        pitch_zmm = 0.95
        pitch_xmm = 2.0
    END IF
    
    IF FWHM_z_padmm < 0.02 AND active_flag% = 1 THEN

        pitch_zmm = 0.01
        pitch_xmm = 0.01
        pitch_rpadmm = 0.01

    END IF
    
    REM for active target the Si detector radius, pitch in z and x direction are fixed
    REM in this case the uncertainty in the z-value measured by the pad detectors is FWHM_z_padmm
    REM if FWHM_z_padmm is very small (0.01mm) then it is assumed that pitch_zmm and pitch_xmm are also very small

    inhomo_Bz = inhomo_Bz_percent / 100. * Bz

    z1det = z1detcm / 100
    z2det = z2detcm / 100

    REM front position and rear position of Si detector wrt target


    pitch_z% = pitch_zmm * 100.

    REM pitch of z strips of Si detector


    pitch_x% = pitch_xmm * 100.

    REM pitch of x strips of Si detector

    pitch_radius_pad% = pitch_rpadmm * 100.
    
    REM pitch of pad sensors
    

    sigma_E4 = fwhm_E4keV / 2350.

    REM uncertainty in ejectile energy measurement


    pi# = 3.1415927#
    z_factor# = SQR(2. * 1.60218E-13 / (1.66054E-27)) * 1.66054E-27 * 2. * pi# / 1.60218E-19

    M3 = M1 + M2 - M4
    E0beam = EoA * M1
    
    ET0 = E0beam + Qvalue
    factor1# = (1. + M1 / M2 * Qvalue / ET0)
    factor2# = (M1 + M2) * (M3 + M4)
    A0# = M1 * M4 / factor2# * E0beam / ET0
    B0# = M1 * M3 / factor2# * E0beam / ET0
    C0# = M2 * M3 / factor2# * factor1#
    D0# = M2 * M4 / factor2# * factor1#
    
    
    
    
    VV = SQR(E0beam ^ 2 + 2 * E0beam * M1 * 931.36814) / (E0beam + M1 * 931.36814) * 2.9979246E+10
    XX = VV / (3.6E+08 * Z_beam ^ .45)
    YY = 1 - .00119 * (Z_target - 6) * SQR(XX) + .00001 * (Z_target - 6) ^ 2 * XX
    charge_state_beam = Z_beam * (1 - EXP(-1.25 * XX + .32 * XX ^ 2 - .11 * XX ^ 3)) * YY
    

    sigma_E0beam = fwhm_E0_percent / 235. * E0beam
    sigma_angle_deg = fwhm_diverg_mrad / 2.35 * (180. / (1000. * pi#))


    Vcm# = M1 / (M1 + M2) * 4.633981637 * SQR(E0beam / M1) * 2997924.6#
    factor4# = M4 * 1.6603145E-27 * (Vcm# * Vcm#) / 2. * 6.241509074E12
    rho_max_factor# = 2. * SQR(2) * SQR(1.602176634E-13 * 1.6603145E-27) / 1.602176634E-19
    Tcyc_factor# = 2 * pi# * 1.6603145E-27 / 1.602176634E-19 * 1.00000E9
    Tcyc# = Tcyc_factor# * M4 * 1.00000E-9 / (charge_state * Bz)
    Helios_factor# = Tcyc# / (1.6603145E-27 * M4 * Vcm# * 6.241509074E12)
    factor5# = 1 / Helios_factor# / 100
    AA# = (M3 + M4) / M3
    BB# = factor4# * AA# - M2 * E0beam / (M1 + M2)
    CC# = factor5# * AA#
    
    REM PRINT #2, "z_factor"; z_factor#; "   rho_max_factor"; rho_max_factor#


    Tcyc_ns = Tcyc# * 1E9
    PRINT #2, "M1 ="; M1; " M2 ="; M2; "M3 ="; M3; " M4 ="; M4; "Beam energy ="; EoA; "MeV/A"; "  Q value ="; Qvalue; "MeV"
    PRINT #2, "Bz ="; Bz; "T"; "  q ="; charge_state; " T cyc ="; Tcyc_ns; "ns"

    tmp$ = "Qvalue (MeV) = ##,.#### *E4 (MeV) + ##,.#### - ##,.#### *z (cm)"
    PRINT #2, USING tmp$; AA#; BB#; CC#

    mass = M4 * 1.66054E-27

    q = charge_state * 1.60218E-19
    
    



    deltat = 1.0E-11

    REM time interval for path integral (s)


    tmp$ = "###,.##  ###,.## ###,.## ###,.#  ##,.##  ##,.##   ##,.##  ##,.##  ##,.##     ##,.##  ##,.##  ##,.###   ###,.# ###,.## ###,.## ###,.## ##,.### ###,.##   ###"



 
    
    IF active_flag% = 0 THEN
        PRINT #2, "radius Si detector ="; Rdetcm; "cm"; "  detector between "; z1detcm; "and "; z2detcm; "cm"; "    radius magnet ="; Rmagnetcm; "cm"; "    inhomogeneity in Bz ="; inhomo_Bz_percent; "%"
        PRINT #2, "Si pitch z ="; pitch_zmm; "mm"; "    pitch x ="; pitch_xmm; "mm"; "     Si FWHM E4 ="; fwhm_E4keV; "keV    E4 threshold ="; E4_thresh_MeV; "MeV"
        PRINT #2, "target thickness ="; target; "mg/cm**2";
        IF ii_scat_max% > 1 THEN
            PRINT #2, "  no. of targets ="; ii_scat_max%
        ELSE
            PRINT #2, "  no. of targets = 1"
        END IF
        PRINT #2, "beam stopping ="; M1_stop; "MeV/mg/cm**2"; "     ejectile stopping power file: "; M4_stop_file$;
        PRINT #2, "  beam stopping flag ="; beam_stop%; "  ejectile stopping flag ="; eject_stop%

    ELSE
        PRINT #2, "inner radius pad detector ="; Rpadcm; "cm"; "   radius Si detector ="; Rdetcm; "cm"; "  detector between "; z1detcm; "and "; z2detcm; "cm"; "    radius magnet ="; Rmagnetcm; "cm"; "    inhomogeneity in Bz ="; inhomo_Bz_percent; "%"
        PRINT #2, "Si pitch z ="; pitch_zmm; "mm"; "    pitch x ="; pitch_xmm; "mm"; "     pad FWHM z ="; FWHM_z_padmm; "mm"; "    pitch pad ="; pitch_rpadmm; "mm"; "     Si FWHM E4 ="; fwhm_E4keV; "keV    E4 threshold ="; E4_thresh_MeV; "MeV"
        PRINT #2, "gas pressure ="; target; "torr"; "  between 0 and "; z_interact_startcm; "cm"
        PRINT #2, "beam stopping ="; M1_stop; "MeV/mg/cm**2"; "     ejectile stopping power file: "; M4_stop_file$;
        PRINT #2, "  beam stopping flag ="; beam_stop%; "  ejectile stopping flag ="; eject_stop%
    END IF
    PRINT #2, "Z beam ="; Z_beam;
    PRINT #2, USING "  beam charge state = ###,.#"; charge_state_beam;
    PRINT #2, "  Z target ="; Z_target; "  A target ="; A_target;
    PRINT #2, "   ejectile multiple scattering flag ="; ii_eject_scat%;
    IF ii_scat_max% > 0 THEN
        PRINT #2, "beam multiple scattering flag = 1"
    ELSE
        PRINT #2, "beam multiple scattering flag = 0"
    END IF
    PRINT #2, "multiple scattering factor ="; scatter_factor
    
    PRINT #2, "FWHM beam straggling ="; fwhm_beam_straggle_percent; "% of energy loss"; "  FWHM ejectile straggling ="; fwhm_eject_straggle_percent; "% of energy loss"
    PRINT #2, "FWHM beam spot ="; fwhm_beam_spot_mm; "mm"; "  FWHM beam energy spread ="; fwhm_E0_percent; "%"; "  FWHM beam divergence ="; fwhm_diverg_mrad; "mrad"
    IF active_flag% = 1 THEN
        PRINT #2, "window foil thickness ="; window_thickness; "mg/cm^2"; "   window position = "; window_positioncm; "cm"
    END IF
    
    v_sum = 0.
    v_diff_sq_sum = 0.



    
    n_pass& = 1
    IF thetacm_step > 1. AND plot% = 0 THEN
        n_pass& = INT(thetacm_step)
        thetacm_step = 0.2
    END IF
    
    IF plot% = 0 THEN thetacm = thetacm_step
    
    
    
    REM if thetacm_step is between 0 and 1 then step over theta
    REM if thetacm_step > 1 then make multiple passes (e.g. to determine efficiency)
    REM if thetacm_step = 0 then plot trajectory for theta = 30 degrees
    REM if thetacm_step < 0 then plot trajectories at theta = |thetacm_step|
    
    PRINT #2, " "
    
    PRINT #2, "number of passes = "; n_pass&
    
    PRINT #2, "theta cm     E4    zeta     t      z     Q value   Dz (pad) (Si)  vertex  corr. t       z     Q value  Delta Q rho_max rad_max distance E4 loss theta  efficiency"
    PRINT #2, " (deg.)     (MeV) (deg.)   (ns)   (cm)    (MeV)     (cm)    (cm)    (cm)       (ns)    (cm)    (MeV)    (keV)    (cm)    (cm)    (cm)   (MeV)   (deg)   (%)"
    
    error_traj_low& = 0
    error_det_overrun& = 0
    error_slow& = 0
    error_large_radius& = 0
    error_dedx& = 0
    error_E4low& = 0
    error_miss_detector& = 0
    error_no_zvertex& = 0
    error_no_converge& = 0
    error_poor_Q& = 0
    loops& = 0
    good_loops& = 0
    
    

        
    DO WHILE thetacm < 175.
    
        tmp2$ = "###,.#"
        IF plot% = 1 THEN
            COLOR 2
        END IF
        LOCATE 7, 5: PRINT USING tmp2$; thetacm
        pass& = 0
        good_pass& = 0
        
        E4_sum = 0
        zeta_sum = 0
        phi_deg_sum = 0
        tns_sum = 0
        z_vertexcm_sum = 0
        zcm_sum = 0
        Qvalue_uncorr_sum = 0
        z_corrcm_sum = 0
        z_changecm_sum = 0
        z_change_padcm_sum = 0
        Qvalue_corr_sum = 0
        tns_corr_sum = 0
        rho_max_cm_sum = 0
        s_travel_tot_sum# = 0
        radius_max = 0
        E4_loss_sum = 0
        thetacm_meas_sum = 0
        REM Qvalue_corr = Qvalue
        


    
        DO WHILE pass& < n_pass&
    
            pass& = pass& + 1
            loops& = loops& + 1
        
            

            

            E0 = gauss(E0beam, sigma_E0beam)

            REM randomise beam energy because of energy spread
        

            IF active_flag% = 0 THEN
                z_vertexcm = 0
            
                tt = RND * target
                IF beam_stop% = 1 THEN E0 = E0 - M1_stop * tt

                REM beam loses energy randomly in the target


                IF target > 0 THEN
                    ii_target% = RND * ii_scat_max% + 1
                    IF ii_target% > ii_scat_max% THEN ii_target% = ii_scat_max%
                    D_thick = (ii_target% - 1) * target + tt
                    sigma_beam_straggle_abs = sigma_beam_straggle * M1_stop * D_thick
                    E0 = gauss(E0, sigma_beam_straggle_abs)
                END IF

                REM select one of targets randomly and randomise beam energy because of straggling
        
        
            ELSE
        
                z_vertexcm = RND * z_interact_startcm
                IF plot% = 1 AND thetacm_step = 180 THEN z_vertexcm = 0.5 * z_interact_startcm
            
                REM for active target reaction occurs beween 0 and z_vertex
            
                gas_lengthcm = z_vertexcm - window_positioncm
 
                REM beam traverses this distance to interaction point
                       

 
                D_thick = gas_lengthcm * density * 1000 + window_thickness
 
                REM effective gas thickness for beam stopping, in mg/cm^2
                REM add nominal amount to take into account gas window foil
            
            
                IF beam_stop% = 1 THEN E0 = E0 - M1_stop * D_thick

                REM beam loses energy in the target
                REM for active target should not correct for beam energy loss as this can be determined.


                sigma_beam_straggle_abs = sigma_beam_straggle * M1_stop * D_thick
                
                REM print #2, M1_stop, gas_lengthcm, D_thick, sigma_beam_straggle_abs
                E0 = gauss(E0, sigma_beam_straggle_abs)

                REM for active target estimate beam energy straggling.
        

            END IF
        

            ET = E0 + Qvalue
            factor1# = (1. + M1 / M2 * Qvalue / ET)
            factor2# = (M1 + M2) * (M3 + M4)
            A# = M1 * M4 / factor2# * E0 / ET
            B# = M1 * M3 / factor2# * E0 / ET
            C# = M2 * M3 / factor2# * factor1#
            D# = M2 * M4 / factor2# * factor1#
    

            E4# = (A# + C# + 2 * SQR(A# * C#) * COS((180. - thetacm) * pi# / 180.)) * ET
            E3# = (B# + D# + 2 * SQR(A# * C#) * COS(thetacm * pi# / 180.)) * ET
            sinpsi# = SIN(thetacm * pi# / 180.) / SQR(E3# / (ET * D#))
            psi# = arcsin#(sinpsi#) * 180. / pi#
            factor3# = E3# * (1. - M3 / M1) + E4# * (1. - M4 / M1) - Qvalue
            cos_zeta_plus_psi# = M1 / (2 * SQR(E3# * E4# * M3 * M4)) * factor3#
            zeta0# = arccos#(cos_zeta_plus_psi#) * 180. / pi# - psi#
            
            IF (z1detcm < 0 AND zeta0# < 91) THEN
                loops& = loops& - 1
                PRINT #2, "zeta approaching 90 degrees, terminate calculation"
                GOTO 11
            END IF
            

            REM 2-body kinematics
            
            E4_start = E4#
            
            REM initial energy of ejectile


            zeta1# = gauss(zeta0#, sigma_angle_deg)

            REM  randomise ejectile angle because of beam divergence


            zeta# = zeta1#
            
            
            rho_max_0 = rho_max_factor# * SQR(E4# * M4) * SIN(zeta# * pi# / 180.) / (charge_state * Bz)
            IF rho_max_0 > Rmagnet + 0.01 THEN
                loops& = loops& - 1
                PRINT #2, "rho_max too large, terminate calculation"
                GOTO 11
            END IF
            
            REM if trajectory radius of ejectile exceeds detector radius by 1cm discontinue theta loop
            
			sigma_scat_beam_deg = 0.
			sigma_scat_beam_foil_deg = 0.

			REM bug corrected 21/7/24 - these variables not initialised if scattering flags set to zero

            IF ii_scat_max% > 0 AND target > 0 THEN
                sigma_scat_beam_deg = multi_scatter(Z_target, A_target, charge_state_beam, E0, D_thick) / scatter_factor * (180. / pi#)
                
                
                IF active_flag% = 1 THEN
                    sigma_scat_beam_foil_deg = multi_scatter(6, 12, charge_state_beam, E0, 0.05) / scatter_factor * (180. / pi#)
                END IF

                REM for active target take into account beam scattering in window foil

                sigma_scat_beam_deg_tot = SQR(sigma_scat_beam_deg ^ 2 + sigma_scat_beam_foil_deg ^ 2)
                zeta# = gauss(zeta1#, sigma_scat_beam_deg_tot)
            END IF
            
                                    

            REM randomise ejectile angle because of multiple scattering of beam

            IF E4# < 0 THEN E4# = 0
            IF E4# < E4_thresh_MeV THEN
                error_E4low& = error_E4low& + 1
                GOTO 10
            END IF
            
            REM check that light particle energy is above threshold


            M4_stop = terpolate(E4#)
            
            IF M4_stop = 0 THEN
                error_dedx& = error_dedx& + 1
                GOTO 10
            END IF
 
            REM stopping power of ejectile


            IF active_flag% = 0 THEN


                IF target > 0 THEN
                    IF zeta# > 90. THEN
                        D_thick = -tt / COS(zeta# * pi# / 180.)
                    ELSE
                        D_thick = (target - tt) / COS(zeta# * pi# / 180.)
                    END IF
                    IF ii_eject_scat% = 1 THEN
                        sigma_scat_deg = multi_scatter(Z_target, A_target, charge_state, E4#, D_thick) / scatter_factor * (180. / pi#)
                        zeta# = gauss(zeta#, sigma_scat_deg)
                    END IF

                    REM multiple scattering of ejectile in the target


                    sigma_eject_straggle_abs = sigma_eject_straggle * M4_stop * D_thick
                    E4# = gauss(E4#, sigma_eject_straggle_abs)

                    REM straggling of ejectile in target


                END IF


                IF eject_stop% = 1 AND target > 0 THEN
                    IF zeta# > 90. THEN
                        E4# = E4# + M4_stop * tt / COS(zeta# * pi# / 180.)
                    ELSE
                        E4# = E4# - M4_stop * (target - tt) / COS(zeta# * pi# / 180.)
                    END IF
                END IF

                REM energy loss of ejectile in the target


            END IF
            IF E4# < 0 THEN E4# = 0
            IF E4# < E4_thresh_MeV THEN
                error_E4low& = error_E4low& + 1
                GOTO 10
            END IF

            REM check that light particle energy is above threshold after losing energy in target
        
        
            5
            sigma_beam_spot_tot = sigma_beam_spot
 
            IF active_flag% = 1 THEN
                scatter_displace = sigma_scat_beam_deg * pi# / 180. * gas_lengthcm / 200
                scatter_displace_foil = sigma_scat_beam_foil_deg * pi# / 180. * gas_lengthcm / 100
                sigma_beam_spot_tot = SQR(sigma_beam_spot ^ 2 + scatter_displace ^ 2 + scatter_displace_foil ^ 2)
            END IF
            
            
            REM for active target take into account displacement in x,y due to beam scattering.
            REM For multiple scattering of beam in gas, reduce gas length by factor of 2 as this incrementally increases
            
            x# = gauss(0., sigma_beam_spot_tot)
            y# = gauss(0., sigma_beam_spot_tot)
            radius = SQR(x# * x# + y# * y#)
            REM IF radius > 2. * sigma_beam_spot_tot THEN GOTO 5
 
            REM randomise starting x,y because of finite beam spot
            
            x_start = x# * 1000
            y_start = y# * 1000

            IF active_flag% = 0 THEN
                z# = 0.
            ELSE
        
                z_sigma_pad = FWHM_z_padmm / 2350.
                z# = z_vertexcm / 100

            
            END IF
        


            v# = SQR(2 * E4# * 1.60218E-13 / mass)

            REM velocity of ejectile+


            IF plot% = 0 THEN
                phi = 2 * pi# * RND
            ELSE
                REM phi = pi#
                phi = 5 * pi# / 4
                REM    (to compare with 46Ar(p,p') in Bradt et al.)
            END IF

            REM random azimuthal angle


            vz# = v# * COS(zeta# * pi# / 180.)
            vx# = v# * SIN(zeta# * pi# / 180.) * COS(phi)
            vy# = v# * SIN(zeta# * pi# / 180.) * SIN(phi)
        

        
            check_pad% = 1
            n% = 0
            nprev% = 0
            t = 0
            s_travel# = 0
            s_travel_tot# = 0
        
            x0# = x#
            y0# = y#
            z0# = z#
            sum_d_theta# = 0.
            sum_d_phi# = 0.


            gas_thick_tot = 0.
            Delta_x_tot# = 0.
            Delta_y_tot# = 0
            Delta_z_tot# = 0
            Delta_s_tot# = 0
            
                
            WHILE (t < Tcyc# * 0.5 OR radius > Rdet)

                REM loop until distance from beam axis < radius of detector. The first condition ensures that particle is returning to the beam axis


                IF ((z1det < 0 AND z# < z1det) OR (z1det > 0 AND z# > z1det)) AND radius < Rdet THEN
                    error_traj_low& = error_traj_low& + 1
                    GOTO 10
                END IF

                REM check that trajectory is not too low for finite size detector


                IF ((z2det > 0 AND z# > z2det) OR (z2det < 0 AND z# < z2det)) THEN
                    error_det_overrun& = error_det_overrun& + 1
                    GOTO 10
                END IF

                REM check that have not overran detector


                IF t > Tcyc# * 2. THEN
                    error_slow& = error_slow& + 1
                    GOTO 10
                END IF

                REM check that trajectory time does not excessively exceed cyclotron time


                IF radius > Rmagnet THEN
                    error_large_radius& = error_large_radius& + 1
                    GOTO 10
                END IF

                REM check that radius of trajectory is within magnet radius
                
                
                z_interact_start = z_interact_startcm / 100
                IF radius > Rpad AND check_pad% = 1 AND active_flag% = 1 AND ((z1detcm < 0 AND z# > (z1detcm + 5)) OR (z1detcm > 0 AND z# < (z1detcm - 5))) THEN
                
                    check_pad% = 0
                    z_meas_pad# = gauss(z#, z_sigma_pad)
                    
                    REM randomise z measurement at pad because of uncertainty of its determination in active target
                
                    radius_i& = radius * 1E5
                    kradius& = INT(radius_i& / pitch_radius_pad%)
                    R_meas_pad# = (kradius& * pitch_radius_pad%) / 1.E5 + 0.5 * pitch_rpadmm / 1000.
                    
                    REM measured radius of pad is binned because of sensor pitch
                                                                  
                END IF

                REM as ejectile trajectory crosses the inner radius of the pad detector, record values of z and R
                                                                                                            
        
                s_iteration# = v# * deltat

                REM distance travelled each iteration in m
            

                s_travel# = s_travel# + s_iteration#
                s_travel_tot# = s_travel_tot# + s_iteration#
                
                REM total distance travelled by ejectile
            
                r_from_target = SQR(x# * x# + y# * y# + z# * z#)
                Bzz = Bz - (r_from_target / 0.5) ^ 2 * inhomo_Bz
                Bx = x * z * inhomo_Bz / 0.5 ^ 2
                By = y * z * inhomo_Bz / 0.5 ^ 2

                REM include effect of field inhomogeneity
                

            
                      
                IF active_flag% = 1 AND s_travel# > s_travel_increment THEN
            
                    gas_thick = s_travel# * 100. * density * 1000

                    REM thickness of gas travelled through in ~ 1cm increments (~ 1mm for isobutane), in mg/cm**2
                

                    gas_thick_tot = gas_thick_tot + gas_thick
            
                    E4_0# = E4#
                    
                    M4_stop = terpolate(E4#)
                    
                    IF M4_stop = 0 THEN
                        error_dedx& = error_dedx& + 1
                        GOTO 10
                    END IF
            
                    IF eject_stop% = 1 THEN
                        E4# = E4# - M4_stop * gas_thick
                    END IF
                
                    REM ejectile loses energy in the gas
                    REM if active target should not correct for ejectile energy loss as this can be determined.
                    
                    IF E4# < 0 THEN E4# = 0
                    IF E4# < E4_thresh_MeV THEN
                        error_E4low& = error_E4low& + 1
                        GOTO 10
                    END IF
                    
                    REM check that light particle energy is above threshold after losing energy in gas
                    
                    
                    sigma_eject_straggle_abs = sigma_eject_straggle * M4_stop * gas_thick
                    E4# = gauss(E4#, sigma_eject_straggle_abs)
                    REM straggling of ejectile in gas

                    IF E4# < 0 THEN E4# = 0
                    IF E4# < E4_thresh_MeV THEN
                        error_E4low& = error_E4low& + 1
                        GOTO 10
                    END IF

                    REM check that light particle energy is above threshold after straggling

            

                
            
            
                    factor = SQR(E4# / E4_0#)
                    vx# = factor * vx#
                    vy# = factor * vy#
                    vz# = factor * vz#
                

                    IF ii_eject_scat% = 1 THEN
            
                        r_travel# = SQR((x# - x0#) ^ 2 + (y# - y0#) ^ 2 + (z# - z0#) ^ 2)
 
                        cos_theta_z = (z# - z0#) / r_travel#
                        sin_theta_z = SQR(1 - cos_theta_z * cos_theta_z)
                        sin_phi_z = (y# - y0#) / (r_travel# * sin_theta_z)
                        cos_phi_z = (x# - x0#) / (r_travel# * sin_theta_z)
                    
                        sigma_scat_rad# = multi_scatter(Z_target, A_target, charge_state, E4#, gas_thick) / scatter_factor
                    
                    
                        d_theta# = gauss(0, sigma_scat_rad#)
                                                                                          
                        sum_d_theta# = sum_d_theta# + d_theta#
                        
                        Delta_x# = (z# - z0#) * cos_phi_z * sum_d_theta#
                        Delta_y# = (z# - z0#) * sin_phi_z * sum_d_theta#
                        Delta_z# = -r_travel# * sin_theta_z * sum_d_theta#
                        
                    
                        x# = x# + Delta_x#
                        y# = y# + Delta_y#
                        z# = z# + Delta_z#
                    
                        
                   
                        REM multiple scattering of ejectile in the gas. The angular change has to be a running sum over the whole path
                    
                    
                        
                        Delta_x_tot# = Delta_x_tot# + Delta_x# * 1000
                        Delta_y_tot# = Delta_y_tot# + Delta_y# * 1000
                        Delta_z_tot# = Delta_z_tot# + Delta_z# * 1000
                        
                    
                        REM tests
 
                    END IF
                    
                    s_travel# = 0
            
                    x0# = x#
                    y0# = y#
                    z0# = z#
                        
                END IF
                99
                       
                vx# = -q * vy# * Bzz * deltat / mass + q * vz# * By * deltat / mass + vx#
                vy# = -q * vz# * Bx * deltat / mass + q * vx# * Bzz * deltat / mass + vy#
                vz# = -q * vx# * By * deltat / mass + q * vy# * Bx * deltat / mass + vz#
                
                IF active_flag% = 1 THEN
                    vz# = vz# - q * Ez * deltat / mass
                END IF

                IF n% = nprev% AND plot% = 1 THEN

                    REM plotting variables

                    x1 = x# * xscale + 150
                    y1 = y# * yscale + 250
                    z1 = z# * zscale + z_plot_offset
                END IF

                x# = vx# * deltat + x#
                y# = vy# * deltat + y#
                z# = vz# * deltat + z#
                radius = SQR(x# * x# + y# * y#)
                
                IF radius > radius_max THEN radius_max = radius

                REM distance from beam axis
  

                n% = n% + 1
                t = t + deltat

                REM total time


                IF n% = nprev% + 10 AND plot% = 1 THEN
 
                    REM plot x-y and z-y every 100 increments
  
                    x2 = x# * xscale + 150
                    y2 = y# * yscale + 250
                    z2 = z# * zscale + z_plot_offset
                    LINE (x1, y1)-(x2, y2), 4
                    REM plots in x, y plane, colour red on white background
                    LINE (z1, y1)-(z2, y2), 4
                    REM plots in z, y plane, colour red on white background
                    REM   plot trajectory
                    REM can remove the x,y or z,y plot for clarity
                    x_plot = x# * 1000
                    y_plot = y# * 1000
                    z_plot = z# * 1000
                    E4_out = E4#
                    PRINT #5, E4_out; "   "; x_plot; "   "; y_plot; "   "; z_plot
                    nprev% = n%
                END IF

            WEND
          
 
            IF (z1det < 0 AND z# > z1det) OR (z1det < 0 AND z# < z2det) OR (z1det > 0 AND z# < z1det) OR (z1det > 0 AND z# > z2det) THEN
                error_miss_detector& = error_miss_detector& + 1
                GOTO 10
            END IF
        
            REM check that intersection is at location of detector
            
            
            IF active_flag% = 1 AND check_pad% = 1 THEN
                error_no_zvertex& = error_no_zvertex& + 1
                GOTO 10
            END IF
            
            REM check that there has been a valid estimate made of the z position of the reaction vertex
 

            vend = SQR(vx# * vx# + vy# * vy# + vz# * vz#)

            REM check that velocity has not changed

            E4_meas# = gauss(E4#, sigma_E4)

            REM randomise measured energy because of detector resolution
            
            IF E4# < 0 THEN E4# = 0
            IF E4_meas# < E4_thresh_MeV THEN
                error_E4low& = error_E4low& + 1
                GOTO 10
            END IF

            REM check that light particle energy is above threshold after randomising energy
            

            
            zi& = z# * 1E5
            kz& = INT(zi& / pitch_z%)
            z_meas# = (kz& * pitch_z%) / 1.E5 + 0.5 * pitch_zmm / 1000.

            xi& = x# * 1E5
            kx& = INT(xi& / pitch_x%)
            x_meas# = (kx& * pitch_x%) / 1.E5 + 0.5 * pitch_xmm / 1000.

            REM bin measured z and x using Si array because of detector pitch
            
            
            Q_corr# = 0
            z_corr# = 0
                        
            iter_max% = 10
            IF active_flag% = 0 THEN
                iter_max% = 2
                z_vertex# = 0
                Delta_z_pad# = 0
                z_meas_pad# = 0
            END IF
                                    
            FOR iter% = 1 TO iter_max%
            
                z_corr_prev# = z_corr#
                
                IF active_flag% = 1 THEN
                    
                    IF iter% < 6 THEN z_corr# = ((AA# * E4_meas# + BB# - Q_corr#) / CC#) / 100.
                    
                    REM estimate of total distance between vertex and intersection of ejectile trajectory with beam axis.  Starting Q value is nominal

                    zeta_meas# = arccos#(z_corr# * charge_state * Bz / (z_factor# * SQR(E4_meas# * M4)))
                    rho_max_meas# = rho_max_factor# * SQR(E4_meas# * M4) * SIN(zeta_meas#) / (charge_state * Bz)
                    alpha# = 2 * arcsin#(R_meas_pad# / rho_max_meas#)
                    Delta_z_pad# = -z_corr# * alpha# / (2. * pi#)
                    z_vertex# = z_meas_pad# + Delta_z_pad#
                    
                    REM print #2, "iter "; iter%; "z_vertex "; z_vertex#
                    
                    REM this extrapolates the value of z at the pad inner radius to the vertex position, using PAB algorithm.
                                
                END IF
                                                                                                                    
                IF active_flag% = 0 THEN z_corr# = z_meas#
                
                R_meas# = SQR(x_meas# * x_meas# + y# * y#)


                FOR jj% = 1 TO 5
                    zeta_meas# = arccos#(z_corr# * charge_state * Bz / (z_factor# * SQR(E4_meas# * M4)))
                    rho_max_meas# = rho_max_factor# * SQR(E4_meas# * M4) * SIN(zeta_meas#) / (charge_state * Bz)
                    alpha# = 2 * arcsin#(R_meas# / rho_max_meas#)
                    Delta_z_Si# = z_corr# * alpha# / (2. * pi#)
                    z_corr# = z_meas# + Delta_z_Si# - z_vertex#
                    
                    REM print #2, "iter "; iter%; "z_corr "; z_corr#
                NEXT
                REM this corrects the measured z at the Si array to the true z, using PAB algorithm
                
            
                Q_corr# = AA# * E4_meas# + BB# - CC# * z_corr# * 100
 
                IF ABS(z_corr# - z_corr_prev#) < 0.0001 THEN GOTO 155
                REM check for convergence
            
            NEXT
            
            
            error_no_converge& = error_no_converge& + 1
            GOTO 10
            REM reject if does not converge after iter_max% iterations
            
            155
            
            z_change_padcm = Delta_z_pad# * 100.
            z_changecm = Delta_z_Si# * 100.
        
            rho_max_cm = rho_max_meas# * 100.
        
            Qvalue_uncorr = AA# * E4_meas# + BB# - CC# * (z_meas# - z_meas_pad#) * 100
            Qvalue_corr = Q_corr#
            
            
            IF ABS(Qvalue_uncorr - Qvalue_corr) > 2. THEN

                REM original criteria was that error flag = 1 if the absolute value of the difference > 10 MeV (August 2022)
                REM reduced to > 1 MeV later, but seems to cut off forward angles for (alpha,p) reaction.  Suggest (4th Jan 2023) 2 MeV as a compromise.

                error_poor_Q& = error_poor_Q& + 1
                GOTO 10
            END IF
            
            
            REM reject if calculate Q value that is clearly in error
            
            z_vertexcm = z_vertex# * 100
            
            phi_deg = phi * 180. / pi#
                       
            tns = t * 1E9
            tns_corr = tns * 2. * pi# / (2. * pi# - alpha#)
        

            IF ABS(Tcyc_ns - tns_corr) > 1 THEN
                error_slow& = error_slow& + 1
                GOTO 10
            END IF

            REM reject slow events at edge of detector
        

            ET_meas# = E0beam + Qvalue_corr
            cosangle = (E4_meas# / ET_meas# - A0# - C0#) / (2 * SQR(A0# * C0#))
            thetacm_meas = 180. - arccos#(cosangle) * 180. / pi#
            
            

            zcm = (z_meas# - z_meas_pad#) * 100
            z_corrcm = z_corr# * 100
            
            good_loops& = good_loops& + 1
            good_pass& = good_pass& + 1

            Q_store_pass(good_pass&) = Qvalue_corr
            Q_store_loops(good_loops&) = Qvalue_corr
            v_diff_sq_sum = (vend - v#) * (vend - v#) + v_diff_sq_sum
            v_sum = v_sum + v#
            
            
            E4_sum = E4_sum + E4#
            zeta_sum = zeta_sum + zeta#
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
            s_travel_tot_sum# = s_travel_tot_sum# + s_travel_tot#
            E4_loss_sum = E4_loss_sum + (E4_start - E4#)
            thetacm_meas_sum = thetacm_meas_sum + thetacm_meas
            
            REM
            REM  IF (thetacm > 10.99 AND thetacm < 11.01) THEN
            REM      PRINT #2, Q_store_pass(good_pass&); "    "; E4#; "    "; E4_meas#; "    "; z_corr#
            REM  END IF
            REM   For small theta_cm the large z correction is strongly correlated to the measured energy in the Si array, resulting in a smaller overall spread in Q values.
            REM
            
            10

        LOOP
        
        IF good_pass& > 0 THEN
 
            efficiency = good_pass& / n_pass& * 100.
            E4_av = E4_sum / good_pass&
            zeta_av = zeta_sum / good_pass&
            phi_deg_av = phi_deg_sum / good_pass&
            tns_av = tns_sum / good_pass&
            z_vertexcm_av = z_vertexcm_sum / good_pass&
            zcm_av = zcm_sum / good_pass&
            Qvalue_uncorr_av = Qvalue_uncorr_sum / good_pass&
            z_corrcm_av = z_corrcm_sum / good_pass&
            z_change_padcm_av = z_change_padcm_sum / good_pass&
            z_changecm_av = z_changecm_sum / good_pass&
            Qvalue_corr_av = Qvalue_corr_sum / good_pass&
            tns_corr_av = tns_corr_sum / good_pass&
            rho_max_cm_av = rho_max_cm_sum / good_pass&
            s_travelcm_tot_av# = s_travel_tot_sum# / good_pass& * 100.
            E4_loss_av = E4_loss_sum / good_pass&
            thetacm_meas_av = thetacm_meas_sum / good_pass&
        
            Q_pass_sum = 0.
 
    
            FOR n% = 1 TO good_pass&
                Q_pass_sum = Q_pass_sum + Q_store_pass(n%)
            NEXT


            Q_pass_average = Q_pass_sum / good_pass&
        
            sigma_pass_sq = 0.
            FOR n% = 1 TO good_pass&
                sigma_pass_sq = sigma_pass_sq + (Q_store_pass(n%) - Q_pass_average) * (Q_store_pass(n%) - Q_pass_average)
            NEXT

            sigma_pass = SQR(sigma_pass_sq / good_pass&) * 1000.
            fwhm_pass = 2.35 * sigma_pass
            
            radius_maxcm = radius_max * 100.

            
        
            PRINT #2, USING tmp$; thetacm; E4_av; zeta_av; tns_av; zcm_av; Qvalue_uncorr_av; z_change_padcm; z_changecm_av; z_vertexcm_av; tns_corr_av; z_corrcm_av; Qvalue_corr_av; fwhm_pass; rho_max_cm_av; radius_maxcm; s_travelcm_tot_av#; E4_loss_av; thetacm_meas_av; efficiency
        ELSE
        
            IF Rdetcm < 1 OR Rpadcm < 1 THEN
                radius_maxcm = radius_max * 100.
                form$ = "theta c. of m. =###,.# degrees    maximum radius =###,.## cm"
                PRINT #2, USING form$; thetacm; radius_maxcm
            END IF
            
        END IF
        
        thetacm = thetacm + thetacm_step
    LOOP
    11

    IF good_loops& = 0 THEN
        PRINT #2, "no good events"
        GOTO 101
    END IF

    v_diff_rms_percent = 100. * SQR(v_diff_sq_sum) / v_sum * SQR(good_loops&)

    Q_sum = 0.
    
       
    FOR n% = 1 TO good_loops&
        Q_sum = Q_sum + Q_store_loops(n%)
    NEXT


    Q_average = Q_sum / good_loops&
    

    PRINT #2,

    sigma_sq = 0.
    FOR n% = 1 TO good_loops&
        sigma_sq = sigma_sq + (Q_store_loops(n%) - Q_average) * (Q_store_loops(n%) - Q_average)
    NEXT
    
    

    PRINT #2,
    
    sigma = SQR(sigma_sq / good_loops&) * 1000.
    fwhm = 2.35 * sigma
    

    

    tmp$ = "average rms velocity error = #,.#### %"
    PRINT #2, USING tmp$; v_diff_rms_percent
    tmp$ = "average measured Q value = ##,.### MeV    standard deviation = ###,.# keV  FWHM = ###,.# keV"
    PRINT #2, USING tmp$; Q_average; sigma; fwhm

    PRINT #2, ""


    101
    

    PRINT #2, "total number of trajectories"; loops&; "  successful trajectories"; good_loops&
    PRINT #2, "low trajectory"; error_traj_low&; "  detector overrun"; error_det_overrun&; "  too slow"; error_slow&;
    PRINT #2, "  too large radius"; error_large_radius&; "  de/dx lookup error"; error_dedx&;
    PRINT #2, "  E4 low"; error_E4low&; "  miss detector"; error_miss_detector&;
    PRINT #2, "  no z vertex"; error_no_zvertex&; "  extrapolation fail"; error_no_converge&; "  poor Q value"; error_poor_Q&
    
    PRINT #2, ""
    PRINT #2, ""
    PRINT #2, ""
    

LOOP

IF plot% = 0 THEN STOP

END

FUNCTION arcsin# (x)
xsq = x * x
IF xsq < 1 THEN
    arcsin# = ATN(x / SQR(1 - (xsq)))
ELSE
    arcsin# = pi# / 2
END IF
END FUNCTION

FUNCTION arccos# (x) ' Inverse Cosine
xsq = x * x
IF xsq < 1 THEN
    arccos# = (2 * ATN(1)) - ATN(x / SQR(1 - xsq))
ELSE
    arccos# = 0
END IF
END FUNCTION

FUNCTION gauss (AM, S)
A = 0.0
FOR I% = 1 TO 12
    A = A + RND
NEXT
gauss = (A - 6.0) * S + AM
END FUNCTION

FUNCTION terpolate (x)

IST% = 0
FOR I% = 1 TO ii_max% - 1
    IF (energy4(I%) <= x AND energy4(I% + 1) >= x) THEN
        IST% = I%
        GOTO 2
    END IF
NEXT
2
IF IST% = 0 THEN
    terpolate = 0
    GOTO 120
END IF

ZM# = LOG(dedx4(IST% + 1) / dedx4(IST%)) / LOG(energy4(IST% + 1) / energy4(IST%))
C# = LOG(dedx4(IST%)) - ZM# * LOG(energy4(IST%))
Y# = ZM# * LOG(x) + C#
terpolate = EXP(Y#)
120
END FUNCTION

FUNCTION multi_scatter (Z, A, Q, T, D)
K# = .0393 / 1000 * Z * (Z + 1) * Q ^ 2 * D / A / T ^ 2
B# = 2.119 * Z ^ (-5 / 3) * (Z + 1) * D / A
LB# = LOG(B#) + 11.51
IF LB# > 0 THEN
    F = SQR(K# * LB#)
ELSE
    F = 0.
END IF
REM multi_scatter = 1.36 * F
REM original statement
multi_scatter = F
REM RMS scattering angle in radians, taken from Kantele's handbook

END FUNCTION



