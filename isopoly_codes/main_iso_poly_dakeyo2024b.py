# Importation required to run this code
import isopoly_solar_wind_solve_and_plot as ipsw
import streamline_calc_dakeyo2024a as stream 
#########################################
# Inputs of the model : Isopoly 
#########################################

# Length of the output model
N = 1e3
L = 1.496e11      # set to 1au by default

# Polytropic indexes
gamma_p_max = 1.52
gamma_e_max = 1.23

# Coronal temperature
Tp0 = 1.63e6
Te0 = 0.71e6

# Isothermal radius (in solar radii)
r_iso_p = 13.6 
r_iso_e = 10.3

# Expansion factor parameters
fm = 2
r_exp = 1.9          # in solar radii
sig_exp = 0.1       # in solar radii
#########################################
# Plotting option 
plot_f = False
plot_gamma = False

plot_unT = True
plot_energy = False
#########################################

###############################################################
# Running the main function
(r, n, u, Tp, Te, gamma_p, gamma_e, f, bol_super) = ipsw.solve_isopoly(
                                        N, L, gamma_p_max, gamma_e_max, 
                                        Tp0, Te0, r_iso_p, r_iso_e,
                                        fm, r_exp, sig_exp, plot_f, 
                                        plot_gamma, plot_unT, 
                                        plot_energy)
###############################################################



#########################################
# Streamline tracing 
#########################################
stream_calc = True
plot_streamline = True
# Probe location for streamline tracing
phi_sat = 10     # in degrees

# Streamline calculation
if(stream_calc):
    (r_phi, phi, v_alf, u_phi) = stream.streamline_calc(r, n, u, phi_sat, plot_streamline)  
###############################################################






