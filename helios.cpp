// HELIOS simulations.  PA Butler  April 2014.  Modified for active gas volume mode June 2022.
// C++ version created June 2024

#include <iostream>
#include <fstream> 
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <random>

using namespace std;

int ii_max = 0, iter, jj, n;
int ii_target;
int beam_stop, eject_stop, ii_scat_max, ii_eject_scat;
int M1, M2, M3, M4, charge_state, Z_beam;

double pi = 4*atan(1);
double energy4[1000];
double dedx4[1000];
double Q_store_pass[100000];
double Q_store_loops[100000];
double z_vertexcm, tt, E4;
double EoA, Qvalue, Bz, Rmagnetcm, inhomo_Bz_percent,Rdetcm, z1detcm, z2detcm, pitch_zmm, pitch_xmm, fwhm_E4keV, E4_thresh_MeV;
double M1_stop, fwhm_beam_straggle_percent, fwhm_eject_straggle_percent, fwhm_beam_spot_mm, fwhm_E0_percent, fwhm_diverg_mrad, thetacm_step;
double thetacm, z, phi, D_thick, gas_lengthcm, zeta, z_sigma_pad, vz, M4_stop, v;
double gas_thick, gas_thick_tot, x1, yy1, z1, z_plot_offset, z_meas, z_vertex, Delta_z_pad, z_meas_pad, converge;
double z_corr_prev, zeta_meas, alpha_c, R_meas_pad, Delta_z_Si, rho_max_meas, sigma_scat_deg, R_meas, z_change_padcm, radius_maxcm;
double Q_corr, Qvalue_corr;
double sigma_scat_beam_deg, sigma_scat_beam_foil_deg, sigma_scat_beam_deg_tot, sigma_eject_straggle_abs;

string dummy, M4_stop_file, M4_stop_filename, line;
vector<string> row;

double xframe = 640.;
double yframe = 350.;
double xscale = xframe;
double yscale = xscale;
double zscale = xscale;

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()


double randomise(double rand_num_max)
{
    double sample, rand_num;

    std::uniform_real_distribution<> dis(0.0, rand_num_max); // Generate random number between 0 and rand_num_max  
	rand_num = dis(gen);
	
    return sample =  rand_num;
}

double gauss(double am, double ss, int check)
{   double sample;
	if (ss/am < 1.0E-8) {
		
		return sample = am;
	}   
    
    // instance of class std::normal_distribution with specific mean and stddev
        std::normal_distribution<double> d(am, ss); 

    // get random number with normal distribution using gen as random source
        return sample = d(gen); 

}

void read_input_line(string line) 
{
 // reading from comma delimited data files
	string word, fin; 
  
        row.clear(); 
  
        // used for breaking words 
        stringstream s(line); 
  
        // read every column data of a row and 
        // store it in a string variable, 'word' 
        while (getline(s, word, ',')) { 
  
            // add all the column data 
            // of a row to a vector 
            row.push_back(word); 
        } 
}

double multi_scatter(double Z, double A, double Q, double T, double D){
 double K, B, LB, F;
 
 // RMS scattering angle in radians, taken from Kantele's handbook
 
 K = .0393 / 1000. * Z * (Z + 1.) * pow(Q,2) * D / A / pow(T,2);
 B = 2.119 * pow(Z,(-5. / 3.)) * (Z + 1.) * D / A;
 LB = log(B) + 11.51;
 if (LB > 0.) {
    F = sqrt(K * LB);
 } else {
 
 F = 0.;
 }
 // multi_scatter = 1.36 * F
 // original statement

return F;
}
  
double interpolate(double x){
      int IST, i;
	  double ZM, C, Y, tpl;
	  
      IST = 0;
      for (i = 1; i < ii_max; i++) {
       if (energy4[i] <= x and energy4[i+1] >= x) {
        IST = i;
        break;
	   }
	  }	
      if (IST == 0){
        tpl = 0.;
        return tpl;
	  }
      ZM = log(dedx4[IST + 1] / dedx4[IST]) / log(energy4[IST + 1] / energy4[IST]);
      C = log(dedx4[IST]) - ZM * log(energy4[IST]);
      Y = ZM * log(x) + C;
      tpl = exp(Y);	  
      return tpl;
}
    
int main() { 
 
 int plot = 0;
 int sintheta_weighting = 0;

 string folder_name, reaction_label, file_name, line_in;
 
 ifstream f("reaction.txt");
 getline (f, folder_name); 
 getline (f, reaction_label);
 
 // folder_name = folder_name[:-1]; // -1 removes end of line

 f.close();

 file_name = folder_name + "/input_" + reaction_label + ".txt" ;
 ifstream file_in(file_name);
 file_name = folder_name + "/output_" + reaction_label + ".txt";
 ofstream file_out(file_name);

 double target = 0.;
 double Ez = 2500 / 0.15;
 // electric field if active target (V/m)

 int idum = 1; 
 
 
 while (file_in.peek() != EOF){
    getline(file_in, line);    	
	getline(file_in, line); read_input_line(line);
	M1 = std::stoi(row[0]); M2 = std::stoi(row[1]); M4 = std::stoi(row[2]); EoA = std::stod(row[3]); Qvalue = std::stod(row[4]);	
	getline(file_in, line); 
	getline(file_in, line); read_input_line(line);
	Bz = std::stod(row[0]); charge_state = std::stoi(row[1]); Rmagnetcm = std::stod(row[2]); inhomo_Bz_percent = std::stod(row[3]);
	getline(file_in, line); 
	getline(file_in, line); read_input_line(line);
	Rdetcm = std::stod(row[0]); z1detcm = std::stod(row[1]); z2detcm = std::stod(row[2]);
	getline(file_in, line); 
	getline(file_in, line); read_input_line(line);
	pitch_zmm = std::stod(row[0]); pitch_xmm = std::stod(row[1]); fwhm_E4keV = std::stod(row[2]); E4_thresh_MeV = std::stod(row[3]);
	getline(file_in, line); 
	getline(file_in, line); read_input_line(line);
	target = std::stod(row[0]); M1_stop = std::stod(row[1]); M4_stop_file = row[2]; beam_stop = std::stoi(row[3]); eject_stop = std::stoi(row[4]);
    getline(file_in, line); 
	getline(file_in, line); read_input_line(line);
	ii_scat_max = std::stoi(row[0]); Z_beam = std::stoi(row[1]); ii_eject_scat = std::stoi(row[2]); fwhm_beam_straggle_percent = std::stod(row[3]); fwhm_eject_straggle_percent = std::stod(row[4]);
    getline(file_in, line); 
	getline(file_in, line); read_input_line(line);
	fwhm_beam_spot_mm = std::stod(row[0]); fwhm_E0_percent = std::stod(row[1]); fwhm_diverg_mrad = std::stod(row[2]);
    getline(file_in, line); 
	getline(file_in, line); read_input_line(line);
	thetacm_step = std::stod(row[0]);
   
    string M4_stop_filename = "srim_files/" + M4_stop_file;
	

	ifstream stop_file_in(M4_stop_filename);

    double Bx = 0.;
    double By = 0.;
    int Z_target = 0;
    int A_target = 0;
    
    if (M4_stop_file == "proton_srim.txt" or M4_stop_file == "deuterium_srim.txt" or M4_stop_file == "carbon_srim.txt" or M4_stop_file == "proton_on_isobutane_srim.txt"){
        Z_target = 6;
	    A_target = 12;
	}
    
    if (M4_stop_file == "proton_on_CD2_srim.txt" or M4_stop_file == "deuteron_on_CD2_srim.txt" or M4_stop_file == "carbon_on_carbon_srim.txt"){
        Z_target = 6;
        A_target = 12;
    }
    
    if (M4_stop_file == "proton_on_deuterium_srim.txt"){ 
        Z_target = 1;
        A_target = 2;
    } 
    
    if (M4_stop_file == "deuteron_on_deuterium_srim.txt" or M4_stop_file == "He-3_on_deuterium_srim.txt"){
        Z_target = 1;
        A_target = 2;
    }   
    
    if (M4_stop_file == "alpha_on_helium_srim.txt" or M4_stop_file == "alpha_on_He-CO2_srim.txt"){
        Z_target = 2;
        A_target = 4;
    } 
    
    if (M4_stop_file == "proton_on_helium_srim.txt" or M4_stop_file == "proton_on_He-CO2_srim.txt"){
        Z_target = 2;
        A_target = 4;
    } 
    
    if (Z_target == 0){
        cout << "target not specified" << endl;
        exit(0);
    }
    
	
    double scatter_factor_active = 1.8;

	double s_travel_increment, density_760;   
    int active_flag = 0;
	
    if (target > 1){ 
	 active_flag = 1;
	}

    // if value of target thickness is greater than 1, then this value is assumed to be the gas pressure (in torr) of an active target
	
	density_760 = 0.;
    //  this statement has been inserted February 2025
        
    if (active_flag == 1){
        s_travel_increment = 0.01;
        density_760 = 0.;
    }      
        
        if (ii_eject_scat < 0){
            scatter_factor_active = 1.0;
            ii_eject_scat = 1;
        }
        // if scatter flag < 0 do not modify multiple scattering for gas target
         
          
    
        if (Z_target == 1 and A_target == 2){
            density_760 = 0.1702 / 1000;
	    }      
        // deuterium
    
        if (Z_target == 2 and A_target == 4){
            density_760 = 0.1662 / 1000;
	    }        
        // helium
    
        if (Z_target == 6 and A_target == 12){
            density_760 = 2.53 / 1000;
            s_travel_increment = 0.001;
        }       
        // isobutane
    
        if (density_760 == 0.){
            cout << "gas density not specified" << endl;
            exit(0);
        }
    

    // gas density in g/cm^3 at 15C, atmospheric pressure

    

        double density = density_760 * target / 760.;

    // actual gas density in g/cm^3
   
  
       
    int ii = 0;
	
    while (getline(stop_file_in, line)){ 
	  read_input_line(line);
      energy4[ii] =  std::stod(row[0]);
	  dedx4[ii] = std::stod(row[1]);
	

	  ii=ii+1;
	}
	
    ii_max = ii;
    
    stop_file_in.close();
    

    double sigma_beam_straggle = fwhm_beam_straggle_percent / 235.;
    double sigma_beam_spot = fwhm_beam_spot_mm / 2350.;
    double sigma_eject_straggle = fwhm_eject_straggle_percent / 235.;




    double Rmagnet = Rmagnetcm / 100.;

    //radius of magnet
    
    double scatter_factor = 1.;
    double z_interact_startcm = 0.;
    double window_positioncm = 0.;
    
    if (active_flag == 1){
    
        scatter_factor = scatter_factor_active;
    
        // this emperically adjusts the value of sigma for multiple scattering so that the Kantele subroutine
        // agrees with the Fano prescription, which reproduces the data of Kuhn et al. NIM B4 (1984) 332
    
        z_interact_startcm = z2detcm;
     
        if (z1detcm < 0.){
            z2detcm = z1detcm - 50.;
            window_positioncm = z1detcm - 5.;
        } else{
            z2detcm = z1detcm + 50.;
            window_positioncm = -5.;
        }
        // window_positioncm is position of gas window, in cm
   
    }
    // for active target the Si detector is assumed to be 50cm long
    // in this case the interaction region is between 0 and z_interact_startcm
    
    double window_thickness = 0.1;
    // for active target this is window foil thickness in mg/cm**2
    
    double pitch_rpadmm = 0.;
    //this line inserted February 2025
   
    
    double Rpadcm = 5.5;
    int cent = 0;

    if (cent == 0 and active_flag == 1){
        Rpadcm = Rdetcm;
        Rdetcm = 3.;
    }
    
    double Rdet = Rdetcm / 100.;

    // radius of Si detector (m)
    
    double Rpad = Rpadcm / 100.;

    // radius of inner pads (m)
    
    double FWHM_z_padmm = 0.;
    pitch_rpadmm = 0.;

    if (active_flag == 1){
        FWHM_z_padmm = pitch_zmm;
        pitch_rpadmm = pitch_xmm;
        pitch_zmm = 0.95;
        pitch_xmm = 2.0;
    }
    
    if (FWHM_z_padmm < 0.02 and active_flag == 1){

        pitch_zmm = 0.01;
        pitch_xmm = 0.01;
        pitch_rpadmm = 0.01;
    }
   
    
    // for active target the Si detector radius, pitch in z and x direction are fixed
    // in this case the uncertainty in the z-value measured by the pad detectors is FWHM_z_padmm
    // if FWHM_z_padmm is very small (0.01mm) then it is assumed that pitch_zmm and pitch_xmm are also very small

    double inhomo_Bz = inhomo_Bz_percent / 100. * Bz;

    double z1det = z1detcm / 100.;
    double z2det = z2detcm / 100.;

    // front position and rear position of Si detector wrt target


    int pitch_z = int(pitch_zmm * 100.);

    // pitch of z strips of Si detector


    int pitch_x = int(pitch_xmm * 100.);

    // pitch of x strips of Si detector

    int pitch_radius_pad = int(pitch_rpadmm * 100.);
    
    // pitch of pad sensors
    

    double sigma_E4 = fwhm_E4keV / 2350.;

    // uncertainty in ejectile energy measurement


    
    double z_factor = sqrt(2. * 1.60218E-13 / (1.66054E-27)) * 1.66054E-27 * 2. * pi / 1.60218E-19;

    M3 = M1 + M2 - M4;
    double E0beam = EoA * static_cast<double>(M1);
    
    double ET0 = E0beam + Qvalue;
    double factor1 = (1. + static_cast<double>(M1) / static_cast<double>(M2) * Qvalue / ET0);
    double factor2 = static_cast<double> ((M1 + M2) * (M3 + M4));
    double A0 = static_cast<double>(M1 * M4) / factor2 * E0beam / ET0;
    double B0 = static_cast<double>(M1 * M3) / factor2 * E0beam / ET0;
    double C0 = static_cast<double>(M2 * M3) / factor2 * factor1;
    double D0 = static_cast<double>(M2 * M4) / factor2 * factor1;
    
    double VV = sqrt(pow(E0beam,2) + 2 * E0beam * static_cast<double>(M1) * 931.36814) / (E0beam + static_cast<double>(M1) * 931.36814) * 2.9979246E+10;
    double XX = VV / (3.6E+08 * pow(static_cast<double>(Z_beam),.45));
    double YY = 1. - .00119 * (static_cast<double>(Z_target) - 6.) * sqrt(XX) + .00001 * pow((static_cast<double>(Z_target) - 6.),2) * XX;
    double charge_state_beam = static_cast<double>(Z_beam) * (1. - exp(-1.25 * XX + .32 * pow(XX,2) - .11 * pow(XX,3))) * YY;
    

    double sigma_E0beam = fwhm_E0_percent / 235. * E0beam;
    double sigma_angle_deg = fwhm_diverg_mrad / 2.35 * (180. / (1000. * pi));


    double Vcm = static_cast<double>(M1) / static_cast<double>(M1 + M2) * 4.633981637 * sqrt(E0beam / static_cast<double>(M1)) * 2997924.6;
    double factor4 = static_cast<double>(M4) * 1.6603145E-27 * (Vcm * Vcm) / 2. * 6.241509074E12;
    double rho_max_factor = 2. * sqrt(2) * sqrt(1.602176634E-13 * 1.6603145E-27) / 1.602176634E-19;
    double Tcyc_factor = 2 * pi * 1.6603145E-27 / 1.602176634E-19 * 1.00000E9;
    double Tcyc = Tcyc_factor * static_cast<double>(M4) * 1.00000E-9 / (static_cast<double>(charge_state) * Bz);
    double Helios_factor = Tcyc / (1.6603145E-27 * static_cast<double>(M4) * Vcm * 6.241509074E12);
    double factor5 = 1 / Helios_factor / 100.;
    double AA = static_cast<double>(M3 + M4) / static_cast<double>(M3);
    double BB = factor4 * AA - static_cast<double>(M2) * E0beam / static_cast<double>(M1 + M2);
    double CC = factor5 * AA;
    
    

    double Tcyc_ns = Tcyc * 1E9;
	
	    
    file_out << "M1 = " << M1 << " M2 = " << M2 << " M3 = " << M3 << " M4 = " << M4 << " Beam energy = " << fixed << setprecision(2) << EoA << " MeV/A" << "  Q value = " << setprecision(1) << Qvalue << " MeV" << endl;
    file_out << "Bz = " << Bz << " T"  << "  q = " << charge_state  << "  T cyc = " << setprecision(4) << Tcyc_ns << " ns" << endl;
    file_out <<  "Qvalue (MeV) = " << AA << " *E4 (MeV) + " << BB << " - " << CC << " *z (cm)" << endl;
	

    double mass = static_cast<double>(M4) * 1.66054E-27;

    double q = static_cast<double>(charge_state) * 1.60218E-19 ;
     
    
    double deltat = 1.0E-11;

    // time interval for path integral (s)

   
    
    
    
    if (active_flag == 0){
		file_out << setprecision(1) << "radius Si detector = " << Rdetcm << "cm" << "  detector between " << z1detcm << " and " << z2detcm << "cm" << "    radius magnet = " << Rmagnetcm << "cm" << "    inhomogeneity in Bz = " << inhomo_Bz_percent << "%" << endl;
		file_out << setprecision(2) << "Si pitch z = " << pitch_zmm << "mm" << "    pitch x = " << pitch_xmm << "mm" << "     Si FWHM E4 = " << fwhm_E4keV << " keV    E4 threshold = " << E4_thresh_MeV << " MeV" << endl;
		file_out << setprecision(1) << "target thickness = " << target << " mg/cm**2 ";
		if (ii_scat_max > 1){
            file_out << "  no. of targets = " << ii_scat_max << endl;
        } else{
            file_out << "  no. of targets = 1" << endl;
        file_out << setprecision(1) << "beam stopping = " << M1_stop << " MeV/mg/cm**2" << "     ejectile stopping power file: " << M4_stop_file;
        file_out << "  beam stopping flag = " << beam_stop << "  ejectile stopping flag = " << eject_stop << endl;
        }
    } else{
        file_out << setprecision(2) << "inner radius pad detector = " << Rpadcm << "cm" << "   radius Si detector = " << Rdetcm << "cm" << "  detector between " << z1detcm << "and " << z2detcm << "cm" << "    radius magnet = " << Rmagnetcm << "cm" << "    inhomogeneity in Bz = " << inhomo_Bz_percent << endl;
        file_out << setprecision(2) << "Si pitch z = " << pitch_zmm << "mm" << "    pitch x = " << pitch_xmm << "mm" << "     pad FWHM z =" << FWHM_z_padmm << "mm" << "    pitch pad = " << pitch_rpadmm << "mm" << "     Si FWHM E4 =" << fwhm_E4keV << " keV    E4 threshold =" << E4_thresh_MeV << " MeV" << endl;
        file_out << setprecision(1) << "gas pressure = " << target << " torr" << "  between 0 and " << z_interact_startcm << "cm" << endl;
        file_out << setprecision(1) << "beam stopping = " << M1_stop << " MeV/mg/cm**2" << "     ejectile stopping power file: " << M4_stop_file;
        file_out << "  beam stopping flag = " << beam_stop << "  ejectile stopping flag = " << eject_stop << endl;
    }
	file_out << "Z beam = " << Z_beam;
    file_out << "  beam charge state = " << charge_state_beam;
    file_out << "  Z target = " << Z_target << "  A target = " << A_target;
    file_out << "   ejectile multiple scattering flag = " << ii_eject_scat;
    if (ii_scat_max > 0){
        file_out << " beam multiple scattering flag = 1 " << endl;
    } else{
        file_out << " beam multiple scattering flag = 0 " << endl;
	}	
    file_out << "multiple scattering factor = " << setprecision(1) << scatter_factor << endl;
    
    file_out << "FWHM beam straggling = " << setprecision(1) << fwhm_beam_straggle_percent << "% of energy loss" << "  FWHM ejectile straggling = " << fwhm_eject_straggle_percent << "% of energy loss" << endl;
    file_out << "FWHM beam spot = " << fwhm_beam_spot_mm << "mm" << "  FWHM beam energy spread = " << fwhm_E0_percent << "%" << "  FWHM beam divergence = " << fwhm_diverg_mrad << "mrad" << endl;
    if (active_flag == 1){
        file_out << "window foil thickness = " << window_thickness << " mg/cm^2" << "   window position = " << window_positioncm << "cm" << endl;
    }
    double v_sum = 0.;
    double v_diff_sq_sum = 0.;
	
    
	
    int n_pass = 1;
    if (thetacm_step > 1. and plot == 0){
        n_pass = int(thetacm_step);
        thetacm_step = 0.2;
    }
    if (plot == 0){ 
	    thetacm = thetacm_step;
    }
    
    
    // if thetacm_step is between 0 and 1 then step over theta
    // if thetacm_step > 1 then make multiple passes (e.g. to determine efficiency)
    // if thetacm_step = 0 then plot trajectory for theta = 30 degrees
    // if thetacm_step < 0 then plot trajectories at theta = |thetacm_step|
    
    file_out << endl;
    
    
    file_out << "number of passes = " << n_pass << endl;
    
    file_out << "theta cm     E4    zeta     t      z     Q value   Dz (pad) (Si)  vertex  corr. t       z     Q value  Delta Q rho_max rad_max distance E4 loss theta  efficiency" << endl;
    file_out << " (deg.)     (MeV) (deg.)   (ns)   (cm)    (MeV)     (cm)    (cm)    (cm)    (ns)       (cm)    (MeV)    (keV)    (cm)    (cm)    (cm)   (MeV)   (deg)   (%)" << endl;
    
    int error_traj_low = 0;
    int error_det_overrun = 0;
    int error_slow = 0;
    int error_large_radius = 0;
    int error_dedx = 0;
    int error_E4low = 0;
    int error_miss_detector = 0;
    int error_no_zvertex = 0;
    int error_no_converge = 0;
    int error_poor_Q = 0;
    int loops = 0;
    int good_loops = 0;
	

    
    
 
   
    while (thetacm < 175.){
    
        
        
		cout << "\r"  << fixed << setprecision(1) << thetacm;
		
        int pass_loop = 0;
        int good_pass = 0;
        
        
        double E4_sum = 0.;
        double zeta_sum = 0.;
        double phi_deg_sum = 0.;
        double tns_sum = 0.;
        double z_vertexcm_sum = 0.;
        double zcm_sum = 0.;
        double Qvalue_uncorr_sum = 0.;
        double z_corrcm_sum = 0.;
        double z_changecm_sum = 0.;
        double z_change_padcm_sum = 0.;
        double Qvalue_corr_sum = 0.;
        double tns_corr_sum = 0.;
        double rho_max_cm_sum = 0.;
        double s_travel_tot_sum = 0.;
        double radius_max = 0.;
        double E4_loss_sum = 0.;
        double thetacm_meas_sum = 0.;
		double sigma_beam_straggle_abs;
        // Qvalue_corr = Qvalue
        


        
        
        while (pass_loop < n_pass){
    
            pass_loop = pass_loop + 1;
            loops = loops + 1;
        
            


            double E0 = gauss(E0beam, sigma_E0beam,1);

            // randomise beam energy because of energy spread

            if (active_flag == 0){
                z_vertexcm = 0.;
            
                tt = randomise(target);
                if (beam_stop == 1){ 
				 E0 = E0 - M1_stop * tt;
				}

                // beam loses energy randomly in the target


                if (target > 0){
                    ii_target = int(randomise(static_cast<double>(ii_scat_max)) + 1.);
                    if (ii_target > ii_scat_max) { ii_target = ii_scat_max;}
                    D_thick = (static_cast<double>(ii_target) - 1.) * target + tt;
					
					
                    sigma_beam_straggle_abs = sigma_beam_straggle * M1_stop * D_thick;
                    if(sigma_beam_straggle_abs/E0 < 1.0E-6) {					
                       
					}
                    E0 = gauss(E0, sigma_beam_straggle_abs, 101);
		
                }

                // select one of targets randomly and randomise beam energy because of straggling
        
        
            } else{
        
             
				z_vertexcm = randomise(1.0)* z_interact_startcm;
				
                if (plot == 1 and thetacm_step == 180){ 
				 z_vertexcm = 0.5 * z_interact_startcm;
                }
                // for active target reaction occurs beween 0 and z_vertex
            
                gas_lengthcm = z_vertexcm - window_positioncm;
 
                // beam traverses this distance to interaction point
                       

 
                D_thick = gas_lengthcm * density * 1000. + window_thickness;
 
                // effective gas thickness for beam stopping, in mg/cm^2
                // add nominal amount to take into account gas window foil
            
            
                if (beam_stop == 1){ 
				 E0 = E0 - M1_stop * D_thick;
                }
                // beam loses energy in the target
                // for active target should not correct for beam energy loss as this can be determined.


                sigma_beam_straggle_abs = sigma_beam_straggle * M1_stop * D_thick;
                
                
                E0 = gauss(E0, sigma_beam_straggle_abs, 2);

                // for active target estimate beam energy straggling.
        
            }
            

            double ET = E0 + Qvalue;
            factor1 = (1. + static_cast<double>(M1) / static_cast<double>(M2) * Qvalue / ET);
            factor2 = static_cast<double>((M1 + M2) * (M3 + M4));
            double A = static_cast<double>(M1 * M4) / factor2 * E0 / ET;
            double B = static_cast<double>(M1 * M3) / factor2 * E0 / ET;
            double C = static_cast<double>(M2 * M3) / factor2 * factor1;
            double D = static_cast<double>(M2 * M4) / factor2 * factor1;
    

            E4 = (A + C + 2 * sqrt(A * C) * cos((180. - thetacm) * pi / 180.)) * ET;
            double E3 = (B + D + 2 * sqrt(A * C) * cos(thetacm * pi / 180.)) * ET;
            double sinpsi = sin(thetacm * pi / 180.) / sqrt(E3 / (ET * D));
            double psi = asin(sinpsi) * 180. / pi;
            double factor3 = E3 * (1. - static_cast<double>(M3) / static_cast<double>(M1)) + E4 * (1. - static_cast<double>(M4) / static_cast<double>(M1)) - Qvalue;
            double cos_zeta_plus_psi = static_cast<double>(M1) / (2 * sqrt(E3 * E4 * static_cast<double>(M3 * M4))) * factor3;
            double zeta0 = acos(cos_zeta_plus_psi) * 180. / pi - psi;
            
            if (z1detcm < 0 and zeta0 < 91){
                loops = loops - 1;
                file_out << endl;
                file_out <<"zeta approaching 90 degrees, terminate calculation" << endl;
                thetacm = 175.;
                break;
            
            }

            // 2-body kinematics
            
            double E4_start = E4;
            
            // initial energy of ejectile


            double zeta1 = gauss(zeta0, sigma_angle_deg, 3);

            //  randomise ejectile angle because of beam divergence


            zeta = zeta1;
            
            
            double rho_max_0 = rho_max_factor * sqrt(E4 * M4) * sin(zeta * pi / 180.) / (static_cast<double>(charge_state) * Bz);
            if (rho_max_0 > Rmagnet + 0.01){
                loops = loops - 1;
                file_out <<"rho_max too large, terminate calculation" << endl;
                thetacm = 175.;
                break;
            }
            
            // if trajectory radius of ejectile exceeds detector radius by 1cm discontinue theta loop
			
            
            sigma_scat_beam_deg = 0.;
			sigma_scat_beam_foil_deg = 0.;
            if (ii_scat_max > 0 and target > 0){
                sigma_scat_beam_deg = multi_scatter(Z_target, A_target, charge_state_beam, E0, D_thick) / scatter_factor * (180. / pi);
                
                
                if (active_flag == 1) {
                    sigma_scat_beam_foil_deg = multi_scatter(6., 12., charge_state_beam, E0, 0.05) / scatter_factor * (180. / pi);
                }

                // for active target take into account beam scattering in window foil

                sigma_scat_beam_deg_tot = sqrt(pow(sigma_scat_beam_deg,2) + pow(sigma_scat_beam_foil_deg,2));
				if(sigma_scat_beam_deg_tot/zeta1 < 1.0E-6) {					
                       
					}
                zeta = gauss(zeta1, sigma_scat_beam_deg_tot, 4);
          
            }
                                    

            // randomise ejectile angle because of multiple scattering of beam

            if (E4 < 0){ 
			 E4 = 0;
			}
            if (E4 < E4_thresh_MeV){
                error_E4low = error_E4low + 1;
                continue;
            }
            
            // check that light particle energy is above threshold


            M4_stop = interpolate(E4);
            
            if (M4_stop == 0){
                error_dedx = error_dedx + 1;
                continue;
            }
 
            // stopping power of ejectile
            
            if (active_flag == 0) {


                if (target > 0){
                    if (zeta > 90.){
                        D_thick = -tt / cos(zeta * pi / 180.);
                    } else{
                        D_thick = (target - tt) / cos(zeta * pi / 180.);
                    }
                    if (ii_eject_scat == 1) {
                        sigma_scat_deg = multi_scatter(Z_target, A_target, static_cast<double>(charge_state), E4, D_thick) / scatter_factor * (180. / pi);
                        zeta = gauss(zeta, sigma_scat_deg, 5);
				    }

                    // multiple scattering of ejectile in the target


                    sigma_eject_straggle_abs = sigma_eject_straggle * M4_stop * D_thick;
                    E4 = gauss(E4, sigma_eject_straggle_abs, 105);

                    // straggling of ejectile in target
			    }



                if (eject_stop == 1 and target > 0.) {
                    if (zeta > 90.) {
                        E4 = E4 + M4_stop * tt / cos(zeta * pi / 180.);
                    } else{
                        E4 = E4 - M4_stop * (target - tt) / cos(zeta * pi / 180.);
					}	
                }

                // energy loss of ejectile in the target

            }
            if (E4 < 0) { 
			 E4 = 0;
			}
            if (E4 < E4_thresh_MeV){
                error_E4low = error_E4low + 1;
                continue;
            }

            // check that light particle energy is above threshold after losing energy in target
        
        
            double scatter_displace, scatter_displace_foil;
            double sigma_beam_spot_tot = sigma_beam_spot;
 
            if (active_flag == 1){
                scatter_displace = sigma_scat_beam_deg * pi / 180. * gas_lengthcm / 200.;
                scatter_displace_foil = sigma_scat_beam_foil_deg * pi / 180. * gas_lengthcm / 100.;
                sigma_beam_spot_tot = sqrt(pow(sigma_beam_spot,2) + pow(scatter_displace,2) + pow(scatter_displace_foil,2));
            }
            
            
            // for active target take into account displacement in x,y due to beam scattering.
            // For multiple scattering of beam in gas, reduce gas length by factor of 2 as this incrementally increases
            
            
			double x = gauss(0., sigma_beam_spot_tot, 6);
            double y = gauss(0., sigma_beam_spot_tot, 7);
            double radius = sqrt(x * x + y * y);
            // IF radius > 2. * sigma_beam_spot_tot THEN GOTO 5
 
            // randomise starting x,y because of finite beam spot
            
            double x_start = x * 1000.;
            double y_start = y * 1000.;

            if (active_flag == 0){
                z = 0.;
            } else{
        
                z_sigma_pad = FWHM_z_padmm / 2350.;
                z = z_vertexcm / 100.;
            }
            v = sqrt(2. * E4 * 1.60218E-13 / mass);

            // velocity of ejectile+


            if (plot == 0) {
                phi =  randomise(2 * pi);
            } else{
                // phi = pi//
                phi = 5. * pi / 4.;
                //    (to compare with 46Ar(p,p') in Bradt et al.)
            }

            // random azimuthal angle


            double vz = v * cos(zeta * pi / 180.);
            double vx = v * sin(zeta * pi / 180.) * cos(phi);
            double vy = v * sin(zeta * pi / 180.) * sin(phi);
        

            
            int check_pad = 1;
            int n = 0;
            int nprev = 0;
            double t = 0.;
            double s_travel = 0.;
            double s_travel_tot = 0.;
        
            double x0 = x;
            double y0 = y;
            double z0 = z;
            double sum_d_theta = 0.;
            double sum_d_phi = 0.;

            
            gas_thick_tot = 0.;
            double Delta_x_tot = 0.;
            double Delta_y_tot = 0;
            double Delta_z_tot = 0;
            double Delta_s_tot = 0;
            
            
            int next_loop = 0;
           
            while (t < Tcyc * 0.5 or radius > Rdet) {

                // loop until distance from beam axis < radius of detector. The first condition ensures that particle is returning to the beam axis


                if (((z1det < 0. and z < z1det) or (z1det > 0. and z > z1det)) and radius < Rdet) {
                    error_traj_low = error_traj_low + 1;
                    next_loop = 1;
                    break;
                }
                // check that trajectory is not too low for finite size detector


                if (((z2det > 0. and z > z2det) or (z2det < 0. and z < z2det))) {
                    error_det_overrun = error_det_overrun + 1;
                    next_loop = 1;
                    break;
                }
                // check that have not overran detector


                if (t > Tcyc * 2.){
                    error_slow = error_slow + 1;
                    next_loop =1;
					
                    break            ;
                }
                // check that trajectory time does not excessively exceed cyclotron time


                if (radius > Rmagnet) {
                    error_large_radius = error_large_radius + 1;
                    next_loop = 1;
                    break;
                }
                // check that radius of trajectory is within magnet radius
                
                
                double z_interact_start = z_interact_startcm / 100.;
                
                if (radius > Rpad and check_pad == 1 and active_flag == 1 and ((z1detcm < 0. and z > (z1detcm + 5.)) or (z1detcm > 0. and z < (z1detcm - 5.)))) {
                
                    check_pad = 0;
                    z_meas_pad = gauss(z, z_sigma_pad, 8);
                    
                    // randomise z measurement at pad because of uncertainty of its determination in active target
                
                    double radius_i = radius * 1.0E5;
                    double kradius = int(radius_i / pitch_radius_pad);
                    R_meas_pad = (kradius * pitch_radius_pad) / 1.E5 + 0.5 * pitch_rpadmm / 1000.;
                    
                    // measured radius of pad is binned because of sensor pitch
                                                                  
                }

                // as ejectile trajectory crosses the inner radius of the pad detector, record values of z and R
                                                                                                            
        
                double s_iteration = v * deltat;

                // distance travelled each iteration in m
            

                s_travel = s_travel + s_iteration;
                s_travel_tot = s_travel_tot + s_iteration;
                
                // total distance travelled by ejectile
            
                double r_from_target = sqrt(x * x + y * y + z * z);
                double Bzz = Bz - pow((r_from_target / 0.5),2) * inhomo_Bz;
                double Bx = x * z * inhomo_Bz / pow(0.5,2);
                double By = y * z * inhomo_Bz / pow(0.5,2);

                // include effect of field inhomogeneity
                

            
                      
                if (active_flag == 1 and s_travel > s_travel_increment) {
            
                    gas_thick = s_travel * 100. * density * 1000.;

                    // thickness of gas travelled through in ~ 1cm increments (~ 1mm for isobutane), in mg/cm**2
                

                    gas_thick_tot = gas_thick_tot + gas_thick;
            
                    double E4_0 = E4;
                    
                    M4_stop = interpolate(E4);
                    
                    if (M4_stop == 0.) {
                        error_dedx = error_dedx + 1;
                        next_loop = 1;
                        break;
                    }
            
                    if (eject_stop == 1) {
                        E4 = E4 - M4_stop * gas_thick;
                    }
                
                    // ejectile loses energy in the gas
                    // if active target should not correct for ejectile energy loss as this can be determined.
                    
                    if (E4 < 0) {  
					 E4 = 0.;
					}
                    if (E4 < E4_thresh_MeV) {
                        error_E4low = error_E4low + 1;
                        next_loop =1;
                        break;
                    }
                    // check that light particle energy is above threshold after losing energy in gas
                    
                    
                    sigma_eject_straggle_abs = sigma_eject_straggle * M4_stop * gas_thick;
                    E4 = gauss(E4, sigma_eject_straggle_abs, 9);
                    // straggling of ejectile in gas

                    if (E4 < 0) {  
					 E4 = 0;
					}
                    if (E4 < E4_thresh_MeV) {
                        error_E4low = error_E4low + 1;
                        next_loop =1;
                        break;
                    }
                    // check that light particle energy is above threshold after straggling

            

                
            
            
                    double factor = sqrt(E4 / E4_0);
                    vx = factor * vx;
                    vy = factor * vy;
                    vz = factor * vz;
                

                    if (ii_eject_scat == 1){
            
                        double r_travel = sqrt(pow((x - x0),2) + pow((y - y0),2) + pow((z - z0),2));
 
                        double cos_theta_z = (z - z0) / r_travel;
                        double sin_theta_z = sqrt(1 - cos_theta_z * cos_theta_z);
                        double sin_phi_z = (y - y0) / (r_travel * sin_theta_z);
                        double cos_phi_z = (x - x0) / (r_travel * sin_theta_z);
                    
                        double sigma_scat_rad = multi_scatter(Z_target, A_target, static_cast<double>(charge_state), E4, gas_thick) / scatter_factor;
                    
                    
                        double d_theta = gauss(0., sigma_scat_rad, 10);
                                                                                          
                        sum_d_theta = sum_d_theta + d_theta;
                        
                        double Delta_x = (z - z0) * cos_phi_z * sum_d_theta;
                        double Delta_y = (z - z0) * sin_phi_z * sum_d_theta;
                        double Delta_z = -r_travel * sin_theta_z * sum_d_theta;
                        
                    
                        x = x + Delta_x;
                        y = y + Delta_y;
                        z = z + Delta_z;
                    
                        
                   
                        // multiple scattering of ejectile in the gas. The angular change has to be a running sum over the whole path
                    
                    
                        
                        Delta_x_tot = Delta_x_tot + Delta_x * 1000.;
                        Delta_y_tot = Delta_y_tot + Delta_y * 1000.;
                        Delta_z_tot = Delta_z_tot + Delta_z * 1000.;
                        
                    
                        // tests
                    }
                   
                    
                    s_travel = 0;
            
                    x0 = x;
                    y0 = y;
                    z0 = z;
				}	
                       
                vx = -q * vy * Bzz * deltat / mass + q * vz * By * deltat / mass + vx;
                vy = -q * vz * Bx * deltat / mass + q * vx * Bzz * deltat / mass + vy;
                vz = -q * vx * By * deltat / mass + q * vy * Bx * deltat / mass + vz;
                
                if (active_flag == 1) {
                    vz = vz - q * Ez * deltat / mass;
                }

                if (n == nprev and plot == 1) {

                    // plotting variables

                    x1 = x * xscale + 150.;
                    yy1 = y * yscale + 250.;
                    z1 = z * zscale + z_plot_offset;
                }

                x = vx * deltat + x;
                y = vy * deltat + y;
                z = vz * deltat + z;
                radius = sqrt(x * x + y * y);
                
                if (radius > radius_max) { 
				 radius_max = radius;
				}

                // distance from beam axis
  

                n = n + 1;
                t = t + deltat;

                // total time


                if (n == nprev + 10 and plot == 1) {
 
                    // plot x-y and z-y every 100 increments
  
                    double x2 = x * xscale + 150.;
                    double y2 = y * yscale + 250.;
                    double z2 = z * zscale + z_plot_offset;
                    // LINE (x1, yy1)-(x2, y2), 4
                    // plots in x, y plane, colour red on white background
                    // LINE (z1, yy1)-(z2, y2), 4
                    // plots in z, y plane, colour red on white background
                    //   plot trajectory
                    // can remove the x,y or z,y plot for clarity
                    double x_plot = x * 1000.;
                    double y_plot = y * 1000.;
                    double z_plot = z * 1000.;
                    double E4_out = E4;
                    // PRINT #5, E4_out; "   "; x_plot; "   "; y_plot; "   "; z_plot
                    nprev = n;
				}	
            }
            if (next_loop == 1) {
                continue;
            }
            
            if ((z1det < 0 and z > z1det) or (z1det < 0 and z < z2det) or (z1det > 0 and z < z1det) or (z1det > 0 and z > z2det)) {
                error_miss_detector = error_miss_detector + 1;
                continue;
            }
            // check that intersection is at location of detector
            
            
            if (active_flag == 1 and check_pad == 1) {
                error_no_zvertex = error_no_zvertex + 1;
                continue;
            }
            // check that there has been a valid estimate made of the z position of the reaction vertex
 

            double vend = sqrt(vx * vx + vy * vy + vz * vz);

            // check that velocity has not changed

            double E4_meas = gauss(E4, sigma_E4, 11);

            // randomise measured energy because of detector resolution
            
            if (E4 < 0) { 
			 E4 = 0;
			}
            if (E4_meas < E4_thresh_MeV) {
                error_E4low = error_E4low + 1;
                continue;
            }
            // check that light particle energy is above threshold after randomising energy
            

            
            long zi = int(z * 1.0E5);
            long kz = int(zi / pitch_z);
            z_meas = static_cast<double>(kz * pitch_z) / 1.0E5 + 0.5 * pitch_zmm / 1000.;

            long xi = int(x * 1.0E5);
            long kx = int(xi / pitch_x);
            double x_meas = static_cast<double>(kx * pitch_x) / 1.0E5 + 0.5 * pitch_xmm / 1000.;

            // bin measured z and x using Si array because of detector pitch
            
            
            Q_corr = 0;
            double z_corr = 0;
                        
            int iter_max = 10;
            if (active_flag == 0) {
                iter_max = 2;
                z_vertex = 0.;
                Delta_z_pad = 0.;
                z_meas_pad = 0.;
            }
            converge = 0; 
            		
            for (iter = 1; iter < iter_max+1; iter++){
            
                z_corr_prev = z_corr;
                
                if (active_flag == 1) {
                    
                    if (iter < 6) { 
					 z_corr = ((AA * E4_meas + BB - Q_corr) / CC) / 100.;
                    }
                    // estimate of total distance between vertex and intersection of ejectile trajectory with beam axis.  Starting Q value is nominal

                    zeta_meas = acos(z_corr * static_cast<double>(charge_state) * Bz / (z_factor * sqrt(E4_meas * M4)));
                    rho_max_meas = rho_max_factor * sqrt(E4_meas * M4) * sin(zeta_meas) / (static_cast<double>(charge_state) * Bz);
                    alpha_c = 2 * asin(R_meas_pad / rho_max_meas);
                    Delta_z_pad = -z_corr * alpha_c / (2. * pi);
                    z_vertex = z_meas_pad + Delta_z_pad;
					
                    
                    // print #2, "iter "; iter%; "z_vertex "; z_vertex#
                    
                    // this extrapolates the value of z at the pad inner radius to the vertex position, using PAB algorithm.
                                
                }
                                                                                                                    
                if (active_flag == 0) { 
				 z_corr = z_meas;
				}
                
                R_meas = sqrt(x_meas * x_meas + y * y);


                
                for (jj=1; jj < 6; jj++) {
                    zeta_meas = acos(z_corr * static_cast<double>(charge_state) * Bz / (z_factor * sqrt(E4_meas * M4)));
                    rho_max_meas = rho_max_factor * sqrt(E4_meas * M4) * sin(zeta_meas) / (static_cast<double>(charge_state) * Bz);
                    alpha_c = 2. * asin(R_meas / rho_max_meas);
                    Delta_z_Si = z_corr * alpha_c / (2. * pi);
                    z_corr = z_meas + Delta_z_Si - z_vertex;
					
					
                    
                    // print #2, "iter "; iter%; "z_corr "; z_corr#
           
                    // this corrects the measured z at the Si array to the true z, using PAB algorithm
                }
            
                Q_corr = AA * E4_meas + BB - CC * z_corr * 100;

                if (fabs(z_corr - z_corr_prev) < 0.0001) {
                 converge = 1;
				 
				 
                 break;
				}
            }    
            
            if (converge == 0) {
             error_no_converge = error_no_converge + 1;
             continue;
			} 
            // reject if does not converge after iter_max% iterations
                    
            
            z_change_padcm = Delta_z_pad * 100.;
            double z_changecm = Delta_z_Si * 100.;
			
			
        
            double rho_max_cm = rho_max_meas * 100.;
        
            double Qvalue_uncorr = AA * E4_meas + BB - CC * (z_meas - z_meas_pad) * 100.;
			
			
			
            Qvalue_corr = Q_corr;
            
            
            if (fabs(Qvalue_uncorr - Qvalue_corr) > 2.) {

                // original criteria was that error flag = 1 if the absolute value of the difference > 10 MeV (August 2022)
                // reduced to > 1 MeV later, but seems to cut off forward angles for (alpha,p) reaction.  Suggest (4th Jan 2023) 2 MeV as a compromise.

                error_poor_Q = error_poor_Q + 1;
                continue;
            }
            
            // reject if calculate Q value that is clearly in error
            
            z_vertexcm = z_vertex * 100.;
            
            double phi_deg = phi * 180. / pi;
                       
            double tns = t * 1.0E9;
            double tns_corr = tns * 2. * pi / (2. * pi - alpha_c);
        

            if (fabs(Tcyc_ns - tns_corr) > 1.) {
                error_slow = error_slow + 1;
				
                continue;
            }
            // reject slow events at edge of detector
                    
            double ET_meas = E0beam + Qvalue_corr;
            double cosangle = (E4_meas / ET_meas - A0 - C0) / (2 * sqrt(A0 * C0));
            double thetacm_meas = 180. - acos(cosangle) * 180. / pi;
            
            

            double zcm = (z_meas - z_meas_pad) * 100.;
            double z_corrcm = z_corr * 100.;
            
            good_loops = good_loops + 1;
            good_pass = good_pass + 1;
            
            Q_store_pass[good_pass] = Qvalue_corr;
            Q_store_loops[good_loops] = Qvalue_corr;

            
            v_diff_sq_sum = (vend - v) * (vend - v) + v_diff_sq_sum;
            v_sum = v_sum + v;
            
            
            E4_sum = E4_sum + E4;
            zeta_sum = zeta_sum + zeta;
            phi_deg_sum = phi_deg_sum + phi_deg;
            tns_sum = tns_sum + tns;
            z_vertexcm_sum = z_vertexcm_sum + z_vertexcm;
            zcm_sum = zcm_sum + zcm;
            Qvalue_uncorr_sum = Qvalue_uncorr_sum + Qvalue_uncorr;
            z_corrcm_sum = z_corrcm_sum + z_corrcm;
            z_change_padcm_sum = z_change_padcm_sum + z_change_padcm;
            z_changecm_sum = z_changecm_sum + z_changecm;
            Qvalue_corr_sum = Qvalue_corr_sum + Qvalue_corr;
            tns_corr_sum = tns_corr_sum + tns_corr;
            rho_max_cm_sum = rho_max_cm_sum + rho_max_cm;
            s_travel_tot_sum = s_travel_tot_sum + s_travel_tot;
            E4_loss_sum = E4_loss_sum + (E4_start - E4);
            thetacm_meas_sum = thetacm_meas_sum + thetacm_meas;
        }    
            
        
        if (good_pass > 0){
 
            double efficiency = static_cast<double>(good_pass) / static_cast<double>(n_pass) * 100.;
            double E4_av = E4_sum / static_cast<double>(good_pass);
            double zeta_av = zeta_sum / static_cast<double>(good_pass);
            double phi_deg_av = phi_deg_sum / static_cast<double>(good_pass);
            double tns_av = tns_sum / static_cast<double>(good_pass);
            double z_vertexcm_av = z_vertexcm_sum / static_cast<double>(good_pass);
            double zcm_av = zcm_sum / static_cast<double>(good_pass);
            double Qvalue_uncorr_av = Qvalue_uncorr_sum / static_cast<double>(good_pass);
			
            double z_corrcm_av = z_corrcm_sum / static_cast<double>(good_pass);
            double z_change_padcm_av = z_change_padcm_sum / static_cast<double>(good_pass);
            double z_changecm_av = z_changecm_sum / static_cast<double>(good_pass);
            double Qvalue_corr_av = Qvalue_corr_sum / static_cast<double>(good_pass);
            double tns_corr_av = tns_corr_sum / static_cast<double>(good_pass);
            double rho_max_cm_av = rho_max_cm_sum / static_cast<double>(good_pass);
            double s_travelcm_tot_av = s_travel_tot_sum / static_cast<double>(good_pass) * 100.;
            double E4_loss_av = E4_loss_sum / static_cast<double>(good_pass);
            double thetacm_meas_av = thetacm_meas_sum / static_cast<double>(good_pass);
        
            double Q_pass_sum = 0.;
 
    
            for (n=1; n < good_pass + 1; n++){
                Q_pass_sum = Q_pass_sum + Q_store_pass[n];
            }

            double Q_pass_average = Q_pass_sum / good_pass;
        
            double sigma_pass_sq = 0.;
            for (n=1; n < good_pass + 1; n++){
                sigma_pass_sq = sigma_pass_sq + (Q_store_pass[n] - Q_pass_average) * (Q_store_pass[n] - Q_pass_average);
            }

            double sigma_pass = sqrt(sigma_pass_sq / static_cast<double>(good_pass)) * 1000.;
            double fwhm_pass = 2.35 * sigma_pass;
            
            radius_maxcm = radius_max * 100.;
            
			file_out << fixed << setprecision(2) << setw(8) << thetacm << setw(8) << E4_av << setw(8) << zeta_av  <<  setprecision(1) << setw(8) << tns_av << setprecision(2) << setw(8) << zcm_av << setw(8) << Qvalue_uncorr_av << setw(8) << z_change_padcm;
			file_out << fixed << setprecision(2) << setw(8) << z_changecm_av << setw(8) << z_vertexcm_av << setw(9) << tns_corr_av  << setprecision(1) << setw(10) << z_corrcm_av << setprecision(3) << setw(10) << Qvalue_corr_av;
			file_out << fixed << setprecision(1) << setw(8) << fwhm_pass << setprecision(2) << setw(8) << rho_max_cm_av << setw(8) << radius_maxcm << setprecision(1) << setw(8) << s_travelcm_tot_av << setprecision(3) << setw(8) << E4_loss_av;
			file_out << fixed << setprecision(2) << setw(8) << thetacm_meas_av << setprecision(0) << setw(6) << efficiency << endl;
			
			
        
 
        }else{
        
            if (Rdetcm < 1. or Rpadcm < 1.) {
                radius_maxcm = radius_max * 100.;

                file_out << "theta c. of m. = " << setprecision(1) << thetacm << "degrees    maximum radius = " << radius_maxcm << " cm" << endl;
            }          
        
        }
        thetacm = thetacm + thetacm_step;
    }
    if (good_loops > 0) {
     
     
  

     double v_diff_rms_percent = 100. * sqrt(v_diff_sq_sum) / v_sum * sqrt(static_cast<double>(good_loops));

     double Q_sum = 0.;
    
       
     for (n =1; n < good_loops+1; n++){ 
        Q_sum = Q_sum + Q_store_loops[n];
   
     }

     double Q_average = Q_sum / static_cast<double>(good_loops);
    

     file_out << endl;

     double sigma_sq = 0.;
     for (n=1; n < good_loops+1; n++){
        sigma_sq = sigma_sq + (Q_store_loops[n] - Q_average) * (Q_store_loops[n] - Q_average);
     }
    

     file_out << endl;
    
     double sigma = sqrt(sigma_sq / static_cast<double>(good_loops)) * 1000.;
     double fwhm = 2.35 * sigma;
    

    

    
     file_out << "average rms velocity error = " << setprecision(4) <<  v_diff_rms_percent << "%" << endl;
    
     file_out << "average measured Q value = " << setprecision(3) << Q_average <<  "MeV    standard deviation = " << setprecision(1) << sigma << " keV  FWHM = " << fwhm << "keV" << endl;

     file_out << endl;
    
    } else{
     file_out << endl;
     file_out << "no good events" << endl;
     file_out << endl;

    }

    file_out << "total number of trajectories " << loops << "  successful trajectories " << good_loops << endl;
    file_out << "low trajectory " << error_traj_low << "  detector overrun " << error_det_overrun << "  too slow " << error_slow;
    file_out << "  too large radius " << error_large_radius << "  de/dx lookup error " << error_dedx;
    file_out << "  E4 low " << error_E4low << "  miss detector " << error_miss_detector;
    file_out << "  no z vertex " << error_no_zvertex << "  extrapolation fail " << error_no_converge << "  poor Q value " << error_poor_Q << endl;
    
    file_out << endl;
    file_out << endl;
    file_out << endl;
 }
 file_in.close();
 file_out.close();
return 0;
} 




      
            
            
            
            
