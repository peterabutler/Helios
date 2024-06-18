# HELIOS simulations.  PA Butler  April 2014.  Modified for active gas volume mode June 2022.
# Python version created June 2024


import sys
import array
import math
import random
random.seed()

def multi_scatter(Z, A, Q, T, D):
 K = .0393 / 1000. * Z * (Z + 1.) * Q ** 2 * D / A / T ** 2
 B = 2.119 * Z ** (-5. / 3.) * (Z + 1.) * D / A
 LB = math.log(B) + 11.51
 if LB > 0.:
    F = math.sqrt(K * LB)
 else:
    F = 0.
 # multi_scatter = 1.36 * F
 # original statement
 return F
 # RMS scattering angle in radians, taken from Kantele's handbook
 
 
 

plot = 0
sintheta_weighting = 0

f = open ("reaction.txt", "r")
folder_name = (f.readline())
folder_name = folder_name[:-1] # -1 removes end of line

reaction_label = (f.readline())

f.close()

file_name = folder_name + "/input_" + reaction_label + ".txt" 
file_in = open(file_name,"r")
file_name = folder_name + "/output_" + reaction_label + ".txt"
file_out = open(file_name,"w")

target = 0.
Ez = 2500 / 0.15
# electric field if active target (V/m)

idum = 1
while idum != 0:
    dummy = (file_in.readline())
    if dummy == "":
     break 
    M1, M2, M4, EoA, Qvalue = [float(x) for x in (file_in.readline()).split(',')]
    dummy = (file_in.readline())
    Bz, charge_state, Rmagnetcm, inhomo_Bz_percent = [float(x) for x in (file_in.readline()).split(',')]
    dummy = (file_in.readline())
    Rdetcm, z1detcm, z2detcm = [float(x) for x in (file_in.readline()).split(',')]  
    dummy = (file_in.readline())
    pitch_zmm, pitch_xmm, fwhm_E4keV, E4_thresh_MeV = [float(x) for x in (file_in.readline()).split(',')]  
    dummy = (file_in.readline())
    target, M1_stop, M4_stop_file, beam_stop, eject_stop = [str(x) for x in (file_in.readline()).split(',')]
    dummy = (file_in.readline())
    ii_scat_max, Z_beam, ii_eject_scat, fwhm_beam_straggle_percent, fwhm_eject_straggle_percent = [float(x) for x in (file_in.readline()).split(',')]
    dummy = (file_in.readline())
    fwhm_beam_spot_mm, fwhm_E0_percent, fwhm_diverg_mrad = [float(x) for x in (file_in.readline()).split(',')]
    dummy = (file_in.readline())
    thetacm_step = (file_in.readline())

    
    
    target = float(target) 
    M1_stop = float(M1_stop) 
    beam_stop = int(beam_stop) 
    eject_stop = int(eject_stop)
    ii_scat_max = int(ii_scat_max)
    ii_eject_scat = int(ii_eject_scat)
    thetacm_step = float(thetacm_step)
    
    
    
    

    M4_stop_filename = "srim_files/" + M4_stop_file

    stop_file_in =open(M4_stop_filename,"r")

    Bx = 0
    By = 0
    Z_target = 0
    A_target = 0
    
    if M4_stop_file == "proton_srim.txt" or M4_stop_file == "deuterium_srim.txt" or M4_stop_file == "carbon_srim.txt" or M4_stop_file == "proton_on_isobutane_srim.txt":
        Z_target = 6.
        A_target = 12.
    
    
    if M4_stop_file == "proton_on_CD2_srim.txt" or M4_stop_file == "deuteron_on_CD2_srim.txt" or M4_stop_file == "carbon_on_carbon_srim.txt":
        Z_target = 6.
        A_target = 12.

    
    if M4_stop_file == "proton_on_deuterium_srim.txt": 
        Z_target = 1.
        A_target = 2.
  
    
    if M4_stop_file == "deuteron_on_deuterium_srim.txt" or M4_stop_file == "He-3_on_deuterium_srim.txt":
        Z_target = 1.
        A_target = 2.
   
    
    if M4_stop_file == "alpha_on_helium_srim.txt" or M4_stop_file == "alpha_on_He-CO2_srim.txt":
        Z_target = 2.
        A_target = 4.
 
    
    if M4_stop_file == "proton_on_helium_srim.txt" or M4_stop_file == "proton_on_He-CO2_srim.txt":
        Z_target = 2.
        A_target = 4.
 
    
    if Z_target == 0:
        print ("target not specified")
        sys.exit(0)
 
    
    scatter_factor_active = 1.8
    
    active_flag = 0
    if target > 1: active_flag = 1

    # if value of target thickness is greater than 1, then this value is assumed to be the gas pressure (in torr) of an active target
        
    if active_flag == 1:
        s_travel_increment = 0.01
        density_760 = 0.
        
        
        if ii_eject_scat < 0:
            scatter_factor_active = 1.0
            ii_eject_scat = 1
  
        # if scatter flag < 0 do not modify multiple scattering for gas target
         
          
    
        if Z_target == 1. and A_target == 2.:
            density_760 = 0.1702 / 1000
        
        # deuterium
    
        if Z_target == 2. and A_target == 4.:
            density_760 = 0.1662 / 1000
        
        # helium
    
        if Z_target == 6. and A_target == 12.:
            density_760 = 2.53 / 1000
            s_travel_increment = 0.001
        
        # isobutane
    
        if density_760 == 0.:
            print ("gas density not specified")
            sys.exit(0)
      
    
    

    # gas density in g/cm^3 at 15C, atmospheric pressure

    

        density = density_760 * target / 760.

    # actual gas density in g/cm^3
   
 
       
    ii = 0
    energy4 = []
    dedx4 = []
    with stop_file_in as inputfile:
     while True:
      line = (inputfile.readline())
      if line =='':
        break
      ii=ii + 1
      energy4_in, dedx4_in = line.split(',')
      energy4.append(float(energy4_in))
      dedx4.append(float(dedx4_in))
    ii_max = ii
    
    stop_file_in.close()
    
    def interpolate(x):
      IST = 0
      for I in range(0, ii_max):
       if energy4[I] <= x and energy4[I+1] >= x:
        IST = I
        break
      if IST == 0:
        tpl = 0.
        return tpl
      ZM = math.log(dedx4[IST + 1] / dedx4[IST]) / math.log(energy4[IST + 1] / energy4[IST])
      C = math.log(dedx4[IST]) - ZM * math.log(energy4[IST])
      Y = ZM * math.log(x) + C
      tpl = math.exp(Y)
      return tpl
    
    

    sigma_beam_straggle = fwhm_beam_straggle_percent / 235.
    sigma_beam_spot = fwhm_beam_spot_mm / 2350.
    sigma_eject_straggle = fwhm_eject_straggle_percent / 235.




    Rmagnet = Rmagnetcm / 100.

    #radius of magnet
    
    scatter_factor = 1.
    z_interact_startcm = 0.
    window_positioncm = 0.
    
    if active_flag == 1:
    
        scatter_factor = scatter_factor_active
    
        # this emperically adjusts the value of sigma for multiple scattering so that the Kantele subroutine
        # agrees with the Fano prescription, which reproduces the data of Kuhn et al. NIM B4 (1984) 332
    
        z_interact_startcm = z2detcm
     
        if (z1detcm < 0.):
            z2detcm = z1detcm - 50.
            window_positioncm = z1detcm - 5.
        else:
            z2detcm = z1detcm + 50.
            window_positioncm = -5.
        
        # window_positioncm is position of gas window, in cm
   
    
    # for active target the Si detector is assumed to be 50cm long
    # in this case the interaction region is between 0 and z_interact_startcm
    
    window_thickness = 0.1
    # for active target this is window foil thickness in mg/cm**2
    

   
    
    Rpadcm = 5.5
    cent = 0

    if cent == 0 and active_flag == 1:
        Rpadcm = Rdetcm
        Rdetcm = 3.
    
    
    Rdet = Rdetcm / 100.

    # radius of Si detector (m)
    
    Rpad = Rpadcm / 100.

    # radius of inner pads (m)
    
    FWHM_z_padmm = 0.  
    pitch_rpadmm = 0.    

    if active_flag == 1:
        FWHM_z_padmm = pitch_zmm
        pitch_rpadmm = pitch_xmm
        pitch_zmm = 0.95
        pitch_xmm = 2.0
   
    
    if FWHM_z_padmm < 0.02 and active_flag == 1:

        pitch_zmm = 0.01
        pitch_xmm = 0.01
        pitch_rpadmm = 0.01

   
    
    # for active target the Si detector radius, pitch in z and x direction are fixed
    # in this case the uncertainty in the z-value measured by the pad detectors is FWHM_z_padmm
    # if FWHM_z_padmm is very small (0.01mm) then it is assumed that pitch_zmm and pitch_xmm are also very small

    inhomo_Bz = inhomo_Bz_percent / 100. * Bz

    z1det = z1detcm / 100.
    z2det = z2detcm / 100.

    # front position and rear position of Si detector wrt target


    pitch_z = int(pitch_zmm * 100.)

    # pitch of z strips of Si detector


    pitch_x = int(pitch_xmm * 100.)

    # pitch of x strips of Si detector

    pitch_radius_pad = int(pitch_rpadmm * 100.)
    
    # pitch of pad sensors
    

    sigma_E4 = fwhm_E4keV / 2350.

    # uncertainty in ejectile energy measurement


    pi = math.pi
    z_factor = math.sqrt(2. * 1.60218E-13 / (1.66054E-27)) * 1.66054E-27 * 2. * pi / 1.60218E-19

    M3 = M1 + M2 - M4
    E0beam = EoA * M1
    
    ET0 = E0beam + Qvalue
    factor1 = (1. + M1 / M2 * Qvalue / ET0)
    factor2 = (M1 + M2) * (M3 + M4)
    A0 = M1 * M4 / factor2 * E0beam / ET0
    B0 = M1 * M3 / factor2 * E0beam / ET0
    C0 = M2 * M3 / factor2 * factor1
    D0 = M2 * M4 / factor2 * factor1
    
    VV = math.sqrt(E0beam ** 2 + 2 * E0beam * M1 * 931.36814) / (E0beam + M1 * 931.36814) * 2.9979246E+10
    XX = VV / (3.6E+08 * Z_beam ** .45)
    YY = 1 - .00119 * (Z_target - 6) * math.sqrt(XX) + .00001 * (Z_target - 6) ** 2 * XX
    charge_state_beam = Z_beam * (1 - math.exp(-1.25 * XX + .32 * XX ** 2 - .11 * XX ** 3)) * YY
    

    sigma_E0beam = fwhm_E0_percent / 235. * E0beam
    sigma_angle_deg = fwhm_diverg_mrad / 2.35 * (180. / (1000. * pi))


    Vcm = M1 / (M1 + M2) * 4.633981637 * math.sqrt(E0beam / M1) * 2997924.6
    factor4 = M4 * 1.6603145E-27 * (Vcm * Vcm) / 2. * 6.241509074E12
    rho_max_factor = 2. * math.sqrt(2) * math.sqrt(1.602176634E-13 * 1.6603145E-27) / 1.602176634E-19
    Tcyc_factor = 2 * pi * 1.6603145E-27 / 1.602176634E-19 * 1.00000E9
    Tcyc = Tcyc_factor * M4 * 1.00000E-9 / (charge_state * Bz)
    Helios_factor = Tcyc / (1.6603145E-27 * M4 * Vcm * 6.241509074E12)
    factor5 = 1 / Helios_factor / 100.
    AA = (M3 + M4) / M3
    BB = factor4 * AA - M2 * E0beam / (M1 + M2)
    CC = factor5 * AA
    
    

    Tcyc_ns = Tcyc * 1E9
    M1_str = str(int(M1))
    M2_str = str(int(M2))
    M3_str = str(int(M3))
    M4_str = str(int(M4))
    EoA_str = str(EoA)
    Qvalue_str = str(Qvalue)
    Bz_str = str(Bz)
    charge_state_str = str(int(charge_state))
    Tcyc_ns_str = str(round(Tcyc_ns,2))
    AA_str = str(round(AA,4))
    BB_str = str(round(BB,4))
    CC_str = str(round(CC,4))
    
    
    
    file_out.write("M1 = " + M1_str+ " M2 = " + M2_str + " M3 = " + M3_str + " M4 = " + M4_str + " Beam energy = " + EoA_str + " MeV/A" + "  Q value = " + Qvalue_str +" MeV\n" )
    file_out.write("Bz =" + Bz_str + " T" + "  q = " + charge_state_str + "  T cyc = " + Tcyc_ns_str + " ns\n")
    file_out.write("Qvalue (MeV) = " + AA_str + " *E4 (MeV) + " + BB_str + " - " + CC_str + " *z (cm) \n")
    
    mass = M4 * 1.66054E-27

    q = charge_state * 1.60218E-19 
    
    
    deltat = 1.0E-11

    # time interval for path integral (s)

   
    
    
    
    if active_flag == 0:
        file_out.write("radius Si detector = " + str(Rdetcm) + "cm" + "  detector between " + str(z1detcm) + " and " + str(z2detcm) + "cm" + "    radius magnet = " + str(Rmagnetcm) + "cm" + "    inhomogeneity in Bz = " + str(inhomo_Bz_percent) + "% \n")
        file_out.write("Si pitch z = " + str(pitch_zmm) + "mm" + "    pitch x = " + str(pitch_xmm) + "mm" + "     Si FWHM E4 = " + str(fwhm_E4keV) + " keV    E4 threshold = " + str(E4_thresh_MeV) + " MeV \n")
        file_out.write("target thickness = " + str(target) + " mg/cm**2 ")
        if ii_scat_max > 1:
            file_out.write("  no. of targets = " + str(ii_scat_max) +"\n")
        else:
            file_out.write("  no. of targets = 1 \n")
        file_out.write("beam stopping = " + str(M1_stop) + " MeV/mg/cm**2" + "     ejectile stopping power file: " + M4_stop_file )
        file_out.write("  beam stopping flag = " + str(beam_stop) + "  ejectile stopping flag = " + str(eject_stop) +"\n")

    else:
        file_out.write("inner radius pad detector = " + str(Rpadcm) + "cm" + "   radius Si detector = " + str(Rdetcm) + "cm" + "  detector between " + str(z1detcm) + "and " + str(z2detcm) + "cm" + "    radius magnet = " + str(Rmagnetcm) + "cm" + "    inhomogeneity in Bz = " + str(inhomo_Bz_percent) + "% \n")
        file_out.write ("Si pitch z = " + str(pitch_zmm) + "mm" + "    pitch x = " + str(pitch_xmm) + "mm" + "     pad FWHM z =" + str(FWHM_z_padmm) + "mm" + "    pitch pad = " + str(pitch_rpadmm) + "mm" + "     Si FWHM E4 =" + str(fwhm_E4keV) + " keV    E4 threshold =" + str(E4_thresh_MeV) + " MeV \n")
        file_out.write("gas pressure = " + str(target) + " torr" + "  between 0 and " + str(z_interact_startcm) + "cm \n")
        file_out.write("beam stopping = " + str(M1_stop) + " MeV/mg/cm**2" + "     ejectile stopping power file: " + M4_stop_file )
        file_out.write("  beam stopping flag = " + str(beam_stop) + "  ejectile stopping flag = " + str(eject_stop) +"\n")
    file_out.write("Z beam = " + str(int(Z_beam)))
    file_out.write("  beam charge state = " + str(int(charge_state_beam)))
    file_out.write("  Z target = " + str(int(Z_target)) + "  A target = " + str(int(A_target)))
    file_out.write("   ejectile multiple scattering flag = " + str(ii_eject_scat))
    if ii_scat_max > 0:
        file_out.write(" beam multiple scattering flag = 1 \n")
    else:
        file_out.write(" beam multiple scattering flag = 0 \n")
    file_out.write("multiple scattering factor = " + str(scatter_factor) +"\n")
    
    file_out.write("FWHM beam straggling = " + str(fwhm_beam_straggle_percent) + "% of energy loss" + "  FWHM ejectile straggling = " + str(fwhm_eject_straggle_percent) + "% of energy loss \n")
    file_out.write("FWHM beam spot = " + str(fwhm_beam_spot_mm) + "mm" + "  FWHM beam energy spread = " + str(fwhm_E0_percent) + "%" + "  FWHM beam divergence = " + str(fwhm_diverg_mrad) + "mrad \n")
    if active_flag == 1:
        file_out.write("window foil thickness = " + str(window_thickness) + " mg/cm^2" + "   window position = " + str(window_positioncm) + "cm \n")
    
    v_sum = 0.
    v_diff_sq_sum = 0.

    n_pass = 1
    if thetacm_step > 1. and plot == 0:
        n_pass = int(thetacm_step)
        thetacm_step = 0.2
    
    if plot == 0: thetacm = thetacm_step
    
    
    
    # if thetacm_step is between 0 and 1 then step over theta
    # if thetacm_step > 1 then make multiple passes (e.g. to determine efficiency)
    # if thetacm_step = 0 then plot trajectory for theta = 30 degrees
    # if thetacm_step < 0 then plot trajectories at theta = |thetacm_step|
    
    file_out.write(" \n")
    
    
    file_out.write("number of passes = " + str(n_pass) + "\n")
    
    file_out.write("theta cm     E4    zeta     t      z     Q value   Dz (pad) (Si)  vertex  corr. t       z     Q value  Delta Q rho_max rad_max distance E4 loss theta  efficiency \n")
    file_out.write(" (deg.)     (MeV) (deg.)   (ns)   (cm)    (MeV)     (cm)    (cm)    (cm)    (ns)       (cm)    (MeV)    (keV)    (cm)    (cm)    (cm)   (MeV)   (deg)   (%) \n")
    
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
    Q_store_loops = []
    
 
   
    while thetacm < 175.:
    
        thetacm_str = str(round(thetacm,1))
        print(thetacm_str, end="\r", flush=True) # this prints at the same location on the screen
        pass_loop = 0
        good_pass = 0
        Q_store_pass = []
        
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
        # Qvalue_corr = Qvalue
        


        
        
        while pass_loop < n_pass:
    
            pass_loop = pass_loop + 1
            loops = loops + 1
        
            

            

            E0 = random.gauss(E0beam, sigma_E0beam)

            # randomise beam energy because of energy spread
        

            if active_flag == 0:
                z_vertexcm = 0.
            
                tt = random.random() * target
                if beam_stop == 1:E0 = E0 - M1_stop * tt

                # beam loses energy randomly in the target


                if target > 0:
                    ii_target = int(random.random() * float(ii_scat_max) + 1.)
                    if ii_target > ii_scat_max: ii_target = ii_scat_max
                    D_thick = (float(ii_target) - 1.) * target + tt
                    sigma_beam_straggle_abs = sigma_beam_straggle * M1_stop * D_thick
                    E0 = random.gauss(E0, sigma_beam_straggle_abs)
               

                # select one of targets randomly and randomise beam energy because of straggling
        
        
            else:
        
                z_vertexcm = random.random() * z_interact_startcm
                if plot == 1 and thetacm_step == 180: z_vertexcm = 0.5 * z_interact_startcm
            
                # for active target reaction occurs beween 0 and z_vertex
            
                gas_lengthcm = z_vertexcm - window_positioncm
 
                # beam traverses this distance to interaction point
                       

 
                D_thick = gas_lengthcm * density * 1000. + window_thickness
 
                # effective gas thickness for beam stopping, in mg/cm^2
                # add nominal amount to take into account gas window foil
            
            
                if beam_stop == 1: E0 = E0 - M1_stop * D_thick

                # beam loses energy in the target
                # for active target should not correct for beam energy loss as this can be determined.


                sigma_beam_straggle_abs = sigma_beam_straggle * M1_stop * D_thick
                
                # print #2, M1_stop, gas_lengthcm, D_thick, sigma_beam_straggle_abs
                E0 = random.gauss(E0, sigma_beam_straggle_abs)

                # for active target estimate beam energy straggling.
        

            

            ET = E0 + Qvalue
            factor1 = (1. + M1 / M2 * Qvalue / ET)
            factor2 = (M1 + M2) * (M3 + M4)
            A = M1 * M4 / factor2 * E0 / ET
            B = M1 * M3 / factor2 * E0 / ET
            C = M2 * M3 / factor2 * factor1
            D = M2 * M4 / factor2 * factor1
    

            E4 = (A + C + 2 * math.sqrt(A * C) * math.cos((180. - thetacm) * pi / 180.)) * ET
            E3 = (B + D + 2 * math.sqrt(A * C) * math.cos(thetacm * pi / 180.)) * ET
            sinpsi = math.sin(thetacm * pi / 180.) / math.sqrt(E3 / (ET * D))
            psi = math.asin(sinpsi) * 180. / pi
            factor3 = E3 * (1. - M3 / M1) + E4 * (1. - M4 / M1) - Qvalue
            cos_zeta_plus_psi = M1 / (2 * math.sqrt(E3 * E4 * M3 * M4)) * factor3
            zeta0 = math.acos(cos_zeta_plus_psi) * 180. / pi - psi
            
            if z1detcm < 0 and zeta0 < 91:
                loops = loops - 1
                file_out.write("\n")
                file_out.write("zeta approaching 90 degrees, terminate calculation \n")
                thetacm = 175.
                break
            
            

            # 2-body kinematics
            
            E4_start = E4
            
            # initial energy of ejectile


            zeta1 = random.gauss(zeta0, sigma_angle_deg)

            #  randomise ejectile angle because of beam divergence


            zeta = zeta1
            
            
            rho_max_0 = rho_max_factor * math.sqrt(E4 * M4) * math.sin(zeta * pi / 180.) / (charge_state * Bz)
            if rho_max_0 > Rmagnet + 0.01:
                loops = loops - 1
                file_out.write("rho_max too large, terminate calculation \n")
                thetacm = 175.
                break
            
            
            # if trajectory radius of ejectile exceeds detector radius by 1cm discontinue theta loop
            

            if ii_scat_max > 0 and target > 0:
                sigma_scat_beam_deg = multi_scatter(Z_target, A_target, charge_state_beam, E0, D_thick) / scatter_factor * (180. / pi)
                sigma_scat_beam_foil_deg = 0.
                
                if active_flag == 1:
                    sigma_scat_beam_foil_deg = multi_scatter(6., 12., charge_state_beam, E0, 0.05) / scatter_factor * (180. / pi)
                

                # for active target take into account beam scattering in window foil

                sigma_scat_beam_deg_tot = math.sqrt(sigma_scat_beam_deg ** 2 + sigma_scat_beam_foil_deg ** 2)
                zeta = random.gauss(zeta1, sigma_scat_beam_deg_tot)
          
            
                                    

            # randomise ejectile angle because of multiple scattering of beam

            if E4 < 0: E4 = 0
            if E4 < E4_thresh_MeV:
                error_E4low = error_E4low + 1
                continue
           
            
            # check that light particle energy is above threshold


            M4_stop = interpolate(E4)
            
            if M4_stop == 0:
                error_dedx = error_dedx + 1
                continue
           
 
            # stopping power of ejectile
            
            if active_flag == 0:


                if target > 0:
                    if zeta > 90.:
                        D_thick = -tt / math.cos(zeta * pi / 180.)
                    else:
                        D_thick = (target - tt) / math.cos(zeta * pi / 180.)
                    
                    if ii_eject_scat == 1:
                        sigma_scat_deg = multi_scatter(Z_target, A_target, charge_state, E4, D_thick) / scatter_factor * (180. / pi)
                        zeta = random.gauss(zeta, sigma_scat_deg)
                   

                    # multiple scattering of ejectile in the target


                    sigma_eject_straggle_abs = sigma_eject_straggle * M4_stop * D_thick
                    E4 = random.gauss(E4, sigma_eject_straggle_abs)

                    # straggling of ejectile in target




                if eject_stop == 1 and target > 0. :
                    if zeta > 90. :
                        E4 = E4 + M4_stop * tt / math.cos(zeta * pi / 180.)
                    else:
                        E4 = E4 - M4_stop * (target - tt) / math.cos(zeta * pi / 180.)
                

                # energy loss of ejectile in the target


            if E4 < 0 : E4 = 0
            if E4 < E4_thresh_MeV:
                error_E4low = error_E4low + 1
                continue
           

            # check that light particle energy is above threshold after losing energy in target
        
        
            
            sigma_beam_spot_tot = sigma_beam_spot
 
            if active_flag == 1:
                scatter_displace = sigma_scat_beam_deg * pi / 180. * gas_lengthcm / 200.
                scatter_displace_foil = sigma_scat_beam_foil_deg * pi / 180. * gas_lengthcm / 100.
                sigma_beam_spot_tot = math.sqrt(sigma_beam_spot ** 2 + scatter_displace ** 2 + scatter_displace_foil ** 2)
            
            
            
            # for active target take into account displacement in x,y due to beam scattering.
            # For multiple scattering of beam in gas, reduce gas length by factor of 2 as this incrementally increases
            
            x = random.gauss(0., sigma_beam_spot_tot)
            y = random.gauss(0., sigma_beam_spot_tot)
            radius = math.sqrt(x * x + y * y)
            # IF radius > 2. * sigma_beam_spot_tot THEN GOTO 5
 
            # randomise starting x,y because of finite beam spot
            
            x_start = x * 1000.
            y_start = y * 1000.

            if active_flag == 0:
                z = 0.
            else:
        
                z_sigma_pad = FWHM_z_padmm / 2350.
                z = z_vertexcm / 100.

            v = math.sqrt(2. * E4 * 1.60218E-13 / mass)

            # velocity of ejectile+


            if plot == 0:
                phi = 2 * pi * random.random()
            else:
                # phi = pi#
                phi = 5. * pi / 4.
                #    (to compare with 46Ar(p,p') in Bradt et al.)
      

            # random azimuthal angle


            vz = v * math.cos(zeta * pi / 180.)
            vx = v * math.sin(zeta * pi / 180.) * math.cos(phi)
            vy = v * math.sin(zeta * pi / 180.) * math.sin(phi)
        

        
            check_pad = 1
            n = 0
            nprev = 0
            t = 0.
            s_travel = 0.
            s_travel_tot = 0.
        
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
            
            
            next_loop = 0
           
            while t < Tcyc * 0.5 or radius > Rdet:

                # loop until distance from beam axis < radius of detector. The first condition ensures that particle is returning to the beam axis


                if ((z1det < 0. and z < z1det) or (z1det > 0. and z > z1det)) and radius < Rdet:
                    error_traj_low = error_traj_low + 1
                    next_loop = 1
                    break
                
                # check that trajectory is not too low for finite size detector


                if ((z2det > 0. and z > z2det) or (z2det < 0. and z < z2det)):
                    error_det_overrun = error_det_overrun + 1
                    next_loop = 1
                    break

                # check that have not overran detector


                if t > Tcyc * 2.:
                    error_slow = error_slow + 1
                    next_loop =1
                    break            

                # check that trajectory time does not excessively exceed cyclotron time


                if radius > Rmagnet:
                    error_large_radius = error_large_radius + 1
                    next_loop = 1
                    break

                # check that radius of trajectory is within magnet radius
                
                
                z_interact_start = z_interact_startcm / 100.
                
                if radius > Rpad and check_pad == 1 and active_flag == 1 and ((z1detcm < 0. and z > (z1detcm + 5.)) or (z1detcm > 0. and z < (z1detcm - 5.))):
                
                    check_pad = 0
                    z_meas_pad = random.gauss(z, z_sigma_pad)
                    
                    # randomise z measurement at pad because of uncertainty of its determination in active target
                
                    radius_i = radius * 1.0E5
                    kradius = int(radius_i / pitch_radius_pad)
                    R_meas_pad = (kradius * pitch_radius_pad) / 1.E5 + 0.5 * pitch_rpadmm / 1000.
                    
                    # measured radius of pad is binned because of sensor pitch
                                                                  
                

                # as ejectile trajectory crosses the inner radius of the pad detector, record values of z and R
                                                                                                            
        
                s_iteration = v * deltat

                # distance travelled each iteration in m
            

                s_travel = s_travel + s_iteration
                s_travel_tot = s_travel_tot + s_iteration
                
                # total distance travelled by ejectile
            
                r_from_target = math.sqrt(x * x + y * y + z * z)
                Bzz = Bz - (r_from_target / 0.5) ** 2 * inhomo_Bz
                Bx = x * z * inhomo_Bz / 0.5 ** 2
                By = y * z * inhomo_Bz / 0.5 ** 2

                # include effect of field inhomogeneity
                

            
                      
                if active_flag == 1 and s_travel > s_travel_increment:
            
                    gas_thick = s_travel * 100. * density * 1000.

                    # thickness of gas travelled through in ~ 1cm increments (~ 1mm for isobutane), in mg/cm**2
                

                    gas_thick_tot = gas_thick_tot + gas_thick
            
                    E4_0 = E4
                    
                    M4_stop = interpolate(E4)
                    
                    if M4_stop == 0.:
                        error_dedx = error_dedx + 1
                        next_loop = 1
                        break
                   
            
                    if eject_stop == 1:
                        E4 = E4 - M4_stop * gas_thick
                    
                
                    # ejectile loses energy in the gas
                    # if active target should not correct for ejectile energy loss as this can be determined.
                    
                    if E4 < 0:  E4 = 0.
                    if E4 < E4_thresh_MeV:
                        error_E4low = error_E4low + 1
                        next_loop =1
                        break
                    
                    # check that light particle energy is above threshold after losing energy in gas
                    
                    
                    sigma_eject_straggle_abs = sigma_eject_straggle * M4_stop * gas_thick
                    E4 = random.gauss(E4, sigma_eject_straggle_abs)
                    # straggling of ejectile in gas

                    if E4 < 0:  E4 = 0
                    if E4 < E4_thresh_MeV:
                        error_E4low = error_E4low + 1
                        next_loop =1
                        break

                    # check that light particle energy is above threshold after straggling

            

                
            
            
                    factor = math.sqrt(E4 / E4_0)
                    vx = factor * vx
                    vy = factor * vy
                    vz = factor * vz
                

                    if ii_eject_scat == 1:
            
                        r_travel = math.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)
 
                        cos_theta_z = (z - z0) / r_travel
                        sin_theta_z = math.sqrt(1 - cos_theta_z * cos_theta_z)
                        sin_phi_z = (y - y0) / (r_travel * sin_theta_z)
                        cos_phi_z = (x - x0) / (r_travel * sin_theta_z)
                    
                        sigma_scat_rad = multi_scatter(Z_target, A_target, charge_state, E4, gas_thick) / scatter_factor
                    
                    
                        d_theta = random.gauss(0., sigma_scat_rad)
                                                                                          
                        sum_d_theta = sum_d_theta + d_theta
                        
                        Delta_x = (z - z0) * cos_phi_z * sum_d_theta
                        Delta_y = (z - z0) * sin_phi_z * sum_d_theta
                        Delta_z = -r_travel * sin_theta_z * sum_d_theta
                        
                    
                        x = x + Delta_x
                        y = y + Delta_y
                        z = z + Delta_z
                    
                        
                   
                        # multiple scattering of ejectile in the gas. The angular change has to be a running sum over the whole path
                    
                    
                        
                        Delta_x_tot = Delta_x_tot + Delta_x * 1000.
                        Delta_y_tot = Delta_y_tot + Delta_y * 1000.
                        Delta_z_tot = Delta_z_tot + Delta_z * 1000.
                        
                    
                        # tests
 
                   
                    
                    s_travel = 0
            
                    x0 = x
                    y0 = y
                    z0 = z
   
                       
                vx = -q * vy * Bzz * deltat / mass + q * vz * By * deltat / mass + vx
                vy = -q * vz * Bx * deltat / mass + q * vx * Bzz * deltat / mass + vy
                vz = -q * vx * By * deltat / mass + q * vy * Bx * deltat / mass + vz
                
                if active_flag == 1:
                    vz = vz - q * Ez * deltat / mass
              

                if n == nprev and plot == 1:

                    # plotting variables

                    x1 = x * xscale + 150.
                    y1 = y * yscale + 250.
                    z1 = z * zscale + z_plot_offset
                

                x = vx * deltat + x
                y = vy * deltat + y
                z = vz * deltat + z
                radius = math.sqrt(x * x + y * y)
                
                if radius > radius_max: radius_max = radius

                # distance from beam axis
  

                n = n + 1
                t = t + deltat

                # total time


                if n == nprev + 10 and plot == 1:
 
                    # plot x-y and z-y every 100 increments
  
                    x2 = x * xscale + 150.
                    y2 = y * yscale + 250.
                    z2 = z * zscale + z_plot_offset
                    # LINE (x1, y1)-(x2, y2), 4
                    # plots in x, y plane, colour red on white background
                    # LINE (z1, y1)-(z2, y2), 4
                    # plots in z, y plane, colour red on white background
                    #   plot trajectory
                    # can #ove the x,y or z,y plot for clarity
                    x_plot = x * 1000.
                    y_plot = y * 1000.
                    z_plot = z * 1000.
                    E4_out = E4
                    # PRINT #5, E4_out; "   "; x_plot; "   "; y_plot; "   "; z_plot
                    nprev = n
 
            if next_loop == 1:
                continue
            
            
            if (z1det < 0 and z > z1det) or (z1det < 0 and z < z2det) or (z1det > 0 and z < z1det) or (z1det > 0 and z > z2det):
                error_miss_detector = error_miss_detector + 1
                continue
        
            # check that intersection is at location of detector
            
            
            if active_flag == 1 and check_pad == 1:
                error_no_zvertex = error_no_zvertex + 1
                continue
            
            # check that there has been a valid estimate made of the z position of the reaction vertex
 

            vend = math.sqrt(vx * vx + vy * vy + vz * vz)

            # check that velocity has not changed

            E4_meas = random.gauss(E4, sigma_E4)

            # randomise measured energy because of detector resolution
            
            if E4 < 0 : E4 = 0
            if E4_meas < E4_thresh_MeV:
                error_E4low = error_E4low + 1
                continue

            # check that light particle energy is above threshold after randomising energy
            

            
            zi = z * 1E5
            kz = int(zi / pitch_z)
            z_meas = (kz * pitch_z) / 1.0E5 + 0.5 * pitch_zmm / 1000.

            xi = x * 1E5
            kx = int(xi / pitch_x)
            x_meas = (kx * pitch_x) / 1.0E5 + 0.5 * pitch_xmm / 1000.

            # bin measured z and x using Si array because of detector pitch
            
            
            Q_corr = 0
            z_corr = 0
                        
            iter_max = 10
            if active_flag == 0:
                iter_max = 2
                z_vertex = 0.
                Delta_z_pad = 0.
                z_meas_pad = 0.
  
            converge = 0                        
            for iter in range(1, iter_max+1):
            
                z_corr_prev = z_corr
                
                if active_flag == 1:
                    
                    if iter < 6 : z_corr = ((AA * E4_meas + BB - Q_corr) / CC) / 100.
                    
                    # estimate of total distance between vertex and intersection of ejectile trajectory with beam axis.  Starting Q value is nominal

                    zeta_meas = math.acos(z_corr * charge_state * Bz / (z_factor * math.sqrt(E4_meas * M4)))
                    rho_max_meas = rho_max_factor * math.sqrt(E4_meas * M4) * math.sin(zeta_meas) / (charge_state * Bz)
                    alpha = 2 * math.asin(R_meas_pad / rho_max_meas)
                    Delta_z_pad = -z_corr * alpha / (2. * pi)
                    z_vertex = z_meas_pad + Delta_z_pad
                    
                    # print #2, "iter "; iter%; "z_vertex "; z_vertex#
                    
                    # this extrapolates the value of z at the pad inner radius to the vertex position, using PAB algorithm.
                                
                
                                                                                                                    
                if active_flag == 0 : z_corr = z_meas
                
                R_meas = math.sqrt(x_meas * x_meas + y * y)


                
                for jj in range(1,6):
                    zeta_meas = math.acos(z_corr * charge_state * Bz / (z_factor * math.sqrt(E4_meas * M4)))
                    rho_max_meas = rho_max_factor * math.sqrt(E4_meas * M4) * math.sin(zeta_meas) / (charge_state * Bz)
                    alpha = 2. * math.asin(R_meas / rho_max_meas)
                    Delta_z_Si = z_corr * alpha / (2. * pi)
                    z_corr = z_meas + Delta_z_Si - z_vertex
                    
                    # print #2, "iter "; iter%; "z_corr "; z_corr#
           
                    # this corrects the measured z at the Si array to the true z, using PAB algorithm
                
            
                Q_corr = AA * E4_meas + BB - CC * z_corr * 100

                if abs(z_corr - z_corr_prev) < 0.0001:
                 converge = 1
                 break
                
            
            if converge == 0:
             error_no_converge = error_no_converge + 1
             continue
            # reject if does not converge after iter_max% iterations
            
        
            
            z_change_padcm = Delta_z_pad * 100.
            z_changecm = Delta_z_Si * 100.
        
            rho_max_cm = rho_max_meas * 100.
        
            Qvalue_uncorr = AA * E4_meas + BB - CC * (z_meas - z_meas_pad) * 100.
            Qvalue_corr = Q_corr
            
            
            if abs(Qvalue_uncorr - Qvalue_corr) > 2.:

                # original criteria was that error flag = 1 if the absolute value of the difference > 10 MeV (August 2022)
                # reduced to > 1 MeV later, but seems to cut off forward angles for (alpha,p) reaction.  Suggest (4th Jan 2023) 2 MeV as a compromise.

                error_poor_Q = error_poor_Q + 1
                continue
           
            
            # reject if calculate Q value that is clearly in error
            
            z_vertexcm = z_vertex * 100.
            
            phi_deg = phi * 180. / pi
                       
            tns = t * 1.0E9
            tns_corr = tns * 2. * pi / (2. * pi - alpha)
        

            if abs(Tcyc_ns - tns_corr) > 1:
                error_slow = error_slow + 1
                continue

            # reject slow events at edge of detector
                    
            ET_meas = E0beam + Qvalue_corr
            cosangle = (E4_meas / ET_meas - A0 - C0) / (2 * math.sqrt(A0 * C0))
            thetacm_meas = 180. - math.acos(cosangle) * 180. / pi
            
            

            zcm = (z_meas - z_meas_pad) * 100.
            z_corrcm = z_corr * 100.
            
            good_loops = good_loops + 1
            good_pass = good_pass + 1
            

            Q_store_pass.append(Qvalue_corr)
            
            Q_store_loops.append(Qvalue_corr)
            
            
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
            
            
        
        if good_pass > 0:
 
            efficiency = good_pass / n_pass * 100.
            E4_av = E4_sum / good_pass
            zeta_av = zeta_sum / good_pass
            phi_deg_av = phi_deg_sum / good_pass
            tns_av = tns_sum / good_pass
            z_vertexcm_av = z_vertexcm_sum / good_pass
            zcm_av = zcm_sum / good_pass
            Qvalue_uncorr_av = Qvalue_uncorr_sum / good_pass
            z_corrcm_av = z_corrcm_sum / good_pass
            z_change_padcm_av = z_change_padcm_sum / good_pass
            z_changecm_av = z_changecm_sum / good_pass
            Qvalue_corr_av = Qvalue_corr_sum / good_pass
            tns_corr_av = tns_corr_sum / good_pass
            rho_max_cm_av = rho_max_cm_sum / good_pass
            s_travelcm_tot_av = s_travel_tot_sum / good_pass * 100.
            E4_loss_av = E4_loss_sum / good_pass
            thetacm_meas_av = thetacm_meas_sum / good_pass
        
            Q_pass_sum = 0.
 
    
            for n in range(0, good_pass):
                Q_pass_sum = Q_pass_sum + Q_store_pass[n]
           

            Q_pass_average = Q_pass_sum / good_pass
        
            sigma_pass_sq = 0.
            for n in range(0, good_pass):
                sigma_pass_sq = sigma_pass_sq + (Q_store_pass[n] - Q_pass_average) * (Q_store_pass[n] - Q_pass_average)
           

            sigma_pass = math.sqrt(sigma_pass_sq / good_pass) * 1000.
            fwhm_pass = 2.35 * sigma_pass
            
            radius_maxcm = radius_max * 100.
            
            file_out.write("%7.2f %7.2f %7.2f %7.1f %7.2f %7.2f %8.2f "% (thetacm,E4_av,zeta_av,tns_av,zcm_av,Qvalue_uncorr_av,z_change_padcm))
            file_out.write("%7.2f %7.2f %8.2f %9.1f %9.3f "% (z_changecm_av,z_vertexcm_av,tns_corr_av,z_corrcm_av,Qvalue_corr_av))
            file_out.write("%7.1f %7.2f %7.2f %7.1f %7.3f "% (fwhm_pass,rho_max_cm_av,radius_maxcm,s_travelcm_tot_av,E4_loss_av))
            file_out.write("%7.2f %5.0f \n"% (thetacm_meas_av,efficiency))
           
         
        
 
        else:
        
            if Rdetcm < 1. or Rpadcm < 1.:
                radius_maxcm = radius_max * 100.
                form = "theta c. of m. = " + str(round(thetacm,1)) + "degrees    maximum radius = " + str(round(radius_maxcm,1)) + " cm \n"
                file_out.write(form)
                      
        
        
        thetacm = thetacm + thetacm_step
 
    if good_loops > 0:
     
     
  

     v_diff_rms_percent = 100. * math.sqrt(v_diff_sq_sum) / v_sum * math.sqrt(good_loops)

     Q_sum = 0.
    
       
     for n in range(0,good_loops): 
        Q_sum = Q_sum + Q_store_loops[n]
   


     Q_average = Q_sum / good_loops
    

     file_out.write("\n")

     sigma_sq = 0.
     for n in range(0,good_loops):
        sigma_sq = sigma_sq + (Q_store_loops[n] - Q_average) * (Q_store_loops[n] - Q_average)
    
    

     file_out.write("\n")
    
     sigma = math.sqrt(sigma_sq / good_loops) * 1000.
     fwhm = 2.35 * sigma
    

    

    
     file_out.write("average rms velocity error = "+ str(round(v_diff_rms_percent,4)) + "% \n")
    
     file_out.write("average measured Q value = " + str(round(Q_average,3)) + "MeV    standard deviation = " + str(round(sigma,1)) + " keV  FWHM = " + str(round(fwhm,1)) + "keV \n")

     file_out.write("\n")
    
    else:
     file_out.write("\n")
     file_out.write("no good events \n")
     file_out.write("\n")



    file_out.write("total number of trajectories " +str(loops) + "  successful trajectories " + str(good_loops) + "\n")
    file_out.write("low trajectory " + str(error_traj_low) + "  detector overrun " + str(error_det_overrun) + "  too slow " +str(error_slow))
    file_out.write("  too large radius " + str(error_large_radius) + "  de/dx lookup error " + str(error_dedx))
    file_out.write("  E4 low " + str(error_E4low) + "  miss detector " + str(error_miss_detector))
    file_out.write("  no z vertex " + str(error_no_zvertex) + "  extrapolation fail " + str(error_no_converge) + "  poor Q value " + str(error_poor_Q) + "\n")
    
    file_out.write("\n")
    file_out.write("\n")
    file_out.write("\n")
    
file_in.close()
file_out.close()




      
            
            
            
            
