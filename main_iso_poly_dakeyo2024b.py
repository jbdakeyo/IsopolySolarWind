# Importation required to run this code
import isopoly_solar_wind_solve_and_plot as ipsw

#########################################
# Inputs of the model 
#########################################

# Length of the output model
N = 5e4
L = 1.496e11      # set to 1au by default

# Polytropic indexes
gamma_p_max = 1.45
gamma_e_max = 1.25

# Coronal temperature
Tp0 = 1.1e6
Te0 = 0.8e6

# Isothermal radius (in solar radii)
r_iso_p = 9 
r_iso_e = 15 

# Expansion factor parameters
fm = 30
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









