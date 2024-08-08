# Package used is the code
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from function_to_run_isopoly_dakeyo2024b import function_to_run_isopoly as func_isopoly


def solve_isopoly(N, L, gamma_p_max, gamma_e_max, Tp0, Te0, r_iso_p, r_iso_e
                                      ,fm, r_exp, sig_exp, plot_f, plot_gamma, plot_unT
                                      , plot_energy):
    
    
    #########################################
    # Initialization of physical quantities
    mp = 1.67e-27
    M = 1.99e30
    G = 6.67e-11
    k = 1.38e-23
    r0 = 6.96e8
    L_min = 10 * r0
    
    nb_vent_ref = 5
    #########################################
    
    def differential(x):
        N = len(x)
        dxdt = np.zeros_like(x)
        for i in range(N-1):
            dxdt[i] = (x[i+1] - x[i])
            dxdt[-1] = dxdt[-2]            
        return dxdt
    
    
    def derive_cen(t,x):
        N = len(x)
        dxdt = np.zeros_like(x)
        for i in range(1, N-1):
            dxdt[i] = (x[i+1] - x[i-1]) / (t[i+1] - t[i-1])
            dxdt[0] = dxdt[1]
            dxdt[-1] = dxdt[-2]
            
        return dxdt
    
    
    
    # Set the dimension of the input quantities that required it
    bol_iso = True
    if(r_iso_p == float('inf')): 
        r_iso_p = L
        bol_iso = True
    if(r_iso_e == float('inf')): 
        r_iso_e = L
        bol_iso = True
    
    r_iso_p = r_iso_p * r0
    r_iso_e = r_iso_e * r0
    r_exp = r_exp * r0
    sig_exp = sig_exp * r0
    
    if(L < L_min):
        print('------------------------------------------------------------------------')
        print('WARNING : L too small --> L < 10*r0 not considered')
        print('--> Calculation done with L = 10*r0')
        print('------------------------------------------------------------------------')
        
        
        
    

    
    
    
    
    
    #########################################
    # Graphic parameters
    ep_trait = 2
    pol = 14
    
    color_reg = [ 'red', 'cyan', 'blue' ]
    xvalues = np.array([1, 10, 100, 200])
    #########################################
    
    N = int(N)
    r = np.geomspace(r0, L, N) 
    
    ind_r_iso_p = np.argmin( abs(r - r_iso_p) ) 
    ind_r_iso_e = np.argmin( abs(r - r_iso_e) ) 
    
    r_iso_min = np.min([r_iso_e, r_iso_p])
    r_iso_max = np.max([r_iso_e, r_iso_p])
    
    
    if(fm<=0): 
        print('------------------------------------------------------------------------')
        print('ERROR : fm <= 0 --> unphysical value')
        print('------------------------------------------------------------------------')
        sys.exit()
    elif(fm<1): 
        print('------------------------------------------------------------------------')
        print('WARNING : fm < 1 --> sub-spherical expansion')
        print('------------------------------------------------------------------------')
    elif( sig_exp<=5e-3*r0 ): 
        print('------------------------------------------------------------------------')
        print('ERROR : sig_exp too small --> increase sig_exp value')
        print('------------------------------------------------------------------------')
        sys.exit()
    elif( r_exp<r0 ): 
        print('------------------------------------------------------------------------')
        print('ERROR : r_exp <= r0 --> unphysical value')
        print('------------------------------------------------------------------------')
        sys.exit()
        
        
    # Calculation of the expansion factor profile (Kopp & Holzer 1976)
    f1 = 1 - (fm -1) * np.exp( (r0 - r_exp)/sig_exp )
    f = ( fm + f1 * np.exp( - ( r - r_exp )/ sig_exp ) ) / ( 1 + np.exp( - (r - r_exp)/sig_exp) )
    
    
    
    
    ##############################################################
    # Computation of the isopoly solutions
    
    (r, rc_iso, u_h, n_h, Tp, Te, cs_T, ind_rc_iso, 
         gamma_p, gamma_e, bol_supersonic) = func_isopoly(r0, Te0, Tp0,
                                          gamma_p_max, gamma_e_max, ind_r_iso_p, ind_r_iso_e,
                                          L ,N, f, bol_iso)
    ##############################################################
    
    
    
    
    ##########################################################"
    # Chargement des mesures
    
    os.chdir( r'C:\Users\jdakeyo\Documents\These_LESIA_Toulouse\Codes_Python\Analyse_donnee_Parker_Solar_Probe_encounter_p_et_e' )
    
    data = np.loadtxt('HELIOS_PSP_med_profile_' + str(nb_vent_ref) + '_vents_r_v_Tp_Te_ne.txt')
    
    
    rmoy_bin = data[:, 0]
    vp_family = data[:, 1 : nb_vent_ref + 1]
    tp_family = data[:, 1*nb_vent_ref +1: 2*nb_vent_ref +1]
    te_family = data[:, 2*nb_vent_ref +1: 3*nb_vent_ref +1]
    np_family = data[:, 3*nb_vent_ref +1: 4*nb_vent_ref +1]
    
    
    tp_family[ vp_family == 0 ] = float('nan')
    te_family[ vp_family == 0 ] = float('nan')
    np_family[ vp_family == 0 ] = float('nan')
    vp_family[ vp_family == 0 ] = float('nan')
    
    
    np_med_1au_dakeyo2022 = np.array([9.46, 9.59, 6.99, 6.14, 5.37])
    
    
    # On ajuste n_h_all du modele au mesures complete
    ##################################################
    ind_r1au = np.argmin( abs(r - L) )
    u_1au_pop = np.array([ vp_family[-1,:] ])
    num_vent = np.argmin( abs(u_1au_pop - u_h[ind_r1au]/1e3) )
    
    #ind_np = np.argwhere( np_family[:,num_vent] >0 )[:,0]
    #np_med_1au = np.median(np_family[ind_np, num_vent] * (rmoy_bin[ind_np] /(L/r0) )**2 )
    n_h = n_h / n_h[-1] * np_med_1au_dakeyo2022[num_vent]
    ##########################################################"
    
    
    
    u_h = u_h / 1e3
    cs_T = cs_T / 1e3
    r = r / r0
    r_iso_p = r_iso_p / r0
    r_iso_e = r_iso_e / r0
    rc_iso = rc_iso / r0
    
    
    
    ind_sel_reg = [0] * 3
    ind_sel_reg[0] = np.argwhere( r < r_iso_min/r0 )
    ind_sel_reg[1] = np.argwhere( (r > r_iso_min/r0) & (r < r_iso_max/r0) )
    ind_sel_reg[2] = np.argwhere( r > r_iso_max/r0 )
    
    
    

    # Plot velocity
    ######################################
    if(plot_unT):
        plt.figure(figsize=(20,5))
        plt.subplot(1,3,1)
        for i in range(len(ind_sel_reg)):
            plt.plot(r[ind_sel_reg[i]], u_h[ind_sel_reg[i]], color = color_reg[i], linewidth = ep_trait )
        plt.scatter(rc_iso, cs_T[ind_rc_iso], s=50, marker='o', color='black', label='Sonic point', zorder=2)
        plt.xlabel('Radius ($r \: / \: r_\\odot$)', fontsize=pol)
        plt.ylabel('u (km/s)', fontsize=pol)
        plt.title('Velocity : $u_0$ = ' + str(np.round(u_h[0],2)) + 
                  ' km/s | $u_{1au}$ = '
                  + str(int(u_h[ind_r1au])) + ' km/s', fontsize=0.9*pol)
        plt.xscale('log')
        plt.xlim([ r[0], 1.05*r[-1] ])
        plt.ylim([0, 1.05*np.max(u_h) ])
        ax = plt.gca()
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        plt.xticks(fontsize= pol)
        plt.yticks(fontsize= pol)
        ax = plt.gca()
        plt.axvspan(1, 2.5, facecolor='silver', alpha=0.5)
        plt.grid()
        plt.legend(fontsize=0.85*pol)
    
    
    
        # Plot density
        ########################################
        plt.subplot(1,3,2)
        for i in range(len(ind_sel_reg)):
            plt.plot(r[ind_sel_reg[i]], n_h[ind_sel_reg[i]], color = color_reg[i], linewidth = ep_trait )
        plt.xlabel('Radius ($r \: / \: r_\\odot$)', fontsize=pol)
        plt.ylabel('n (#/$cm^{-3}$)', fontsize=pol)
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Density : $n_0$ = ' + '%.2E' % int(n_h[0]) + ' #/$cm^{-3}$ | $n_{1au}$ = '
                 + str(int(n_h[ind_r1au])) + ' #/$cm^{-3}$' , fontsize=0.9*pol)
        plt.xlim([ r[0], 1.05*r[-1] ])
        plt.ylim([ 0.5 * np.min(n_h), 1.5 * np.max(n_h) ])
        ax = plt.gca()
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        plt.xticks(fontsize= pol)
        plt.yticks(fontsize= pol)
        ax = plt.gca()
        plt.axvspan(1, 2.5, facecolor='silver', alpha=0.5)
        plt.grid()
        
    
    
        # Plot temperature
        ###############################################
        plt.subplot(1,3,3)
        plt.plot(r[r<r_iso_p], Tp[r<r_iso_p], color = color_reg[0], linewidth = ep_trait )
        plt.plot(r[r>=r_iso_p], Tp[r>=r_iso_p], color = color_reg[-1], linewidth = ep_trait )
        plt.plot(r[r<r_iso_e], Te[r<r_iso_e], '--', color = color_reg[0], linewidth = ep_trait )
        plt.plot(r[r>=r_iso_e], Te[r>=r_iso_e], '--', color = color_reg[-1], linewidth = ep_trait )
        plt.xlabel('Radius ($r \: / \: r_\\odot$)', fontsize=pol)
        plt.ylabel('Température (MK)', fontsize=pol)
        plt.title('Temperature : $T_{p|1au}$ = ' + '%.2E' % int(Tp[ind_r1au]) + 
                  'K | $T_{e|1au}$ = ' + '%.2E' % int(Te[ind_r1au]) + 'K ', fontsize=0.9*pol)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim([ r[0], 1.05*r[-1] ])
        plt.ylim([ 0.75*np.min([np.min(Tp), np.min(Te)]),1.5 * np.max([np.max(Tp), np.max(Te)]) ])
        ax = plt.gca()
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        plt.xticks(fontsize= pol)
        plt.yticks(fontsize= pol)
        ax = plt.gca()
        plt.axvspan(1, 2.5, facecolor='silver', alpha=0.5)
        plt.grid()
    
    
    
    # Plot f
    #########################################
    if(plot_f):
        plt.figure()
        plt.plot(r, f, color='black', linewidth=1.7)
        plt.xlabel('Radius ($r \: / \: r_0$)', fontsize=pol)
        plt.ylabel('f(r)', fontsize=pol)
        plt.title('Expansion factor profile', fontsize = 0.9*pol)
        #plt.xscale('log')
        plt.yscale('log')
        plt.xlim([1, 5])
        ax = plt.gca()
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        plt.xticks(fontsize= pol)
        plt.yticks(fontsize= pol)
        plt.grid()
    
    
    # Plot gamma 
    #########################################
    if(plot_gamma):
        plt.figure( figsize=(7,5))
        plt.plot(r[r<=r_iso_p], gamma_p[r<=r_iso_p], color = color_reg[0], linewidth = ep_trait )
        plt.plot(r[r>=0.99*r_iso_p], gamma_p[r>=0.99*r_iso_p], color = color_reg[-1], linewidth = ep_trait, label='$\\gamma_{p|max}$' )
        plt.plot(r[r<=r_iso_e], gamma_e[r<=r_iso_e], '--', color = color_reg[0], linewidth = ep_trait )
        plt.plot(r[r>=0.99*r_iso_e], gamma_e[r>=0.99*r_iso_e], '--', color = color_reg[-1], linewidth = ep_trait, label='$\\gamma_{e|max}$' )
        plt.xlabel('Radius ($r \: / \: r_\\odot$)', fontsize=pol)
        plt.ylabel('$\\gamma$', fontsize=pol)
        plt.xscale('log')
        plt.title('Polytropic indices profiles', fontsize = 0.9*pol)
        plt.xlim([ r[0], 1.05*r[-1] ])
        plt.ylim([0.9 , 1.05 * np.max([gamma_p_max , gamma_e_max])])
        ax = plt.gca()
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        plt.xticks(fontsize= pol)
        plt.yticks(fontsize= pol)
        ax = plt.gca()
        plt.axvspan(1, 2.5, facecolor='silver', alpha=0.5)
        plt.grid()
        plt.legend()
        
    
    # Calculation of the mechanical energy
    ##################################################
    r = r * r0
    u_h = u_h * 1e3
    dh = differential(r) 
    
    
    E_cin = mp / 2. * u_h**2
    E_grav = - G * M * mp / r
    E_press = np.zeros_like(E_cin)
    E_cin_int = np.zeros_like(E_cin)
    
    # Work of pressure force
    n_tilde = n_h / n_h[ind_rc_iso] 
    
    n_tilde_p = n_h / n_h[ind_r_iso_p]
    n_tilde_e = n_h / n_h[ind_r_iso_e]  
    
    F_press = np.zeros_like(E_press)
    F_press_p = np.zeros_like(E_press)
    F_press_e = np.zeros_like(E_press)
    
    
    i = 0
    F_press_p[i] = ( k / n_tilde_p[i]  ) \
        * ( Tp[i+1] * n_tilde_p[i+1] - Tp[i] * n_tilde_p[i] ) / (r[i+1] - r[i]) 
        
    F_press_e[i] = ( k / n_tilde_e[i] ) \
        * ( Te[i+1] * n_tilde_e[i+1] - Te[i] * n_tilde_e[i] ) / (r[i+1] - r[i])
    for i in range(1, N-1):
        F_press_p[i] = ( k / n_tilde_p[i]  ) \
            * ( Tp[i+1] * n_tilde_p[i+1] - Tp[i-1] * n_tilde_p[i-1] ) / (r[i+1] - r[i-1]) 
            
        F_press_e[i] = ( k / n_tilde_e[i] ) \
            * ( Te[i+1] * n_tilde_e[i+1] - Te[i-1] * n_tilde_e[i-1] ) / (r[i+1] - r[i-1])
    
    F_press_p[-1] = ( k / n_tilde_p[i]  ) \
        * ( Tp[i] * n_tilde_p[i] - Tp[i-1] * n_tilde_p[i-1] ) / (r[i] - r[i-1]) 
        
    F_press_e[-1] = ( k / n_tilde_e[i] ) \
        * ( Te[i] * n_tilde_e[i] - Te[i-1] * n_tilde_e[i-1] ) / (r[i] - r[i-1])
                
            
    E_press_p = np.cumsum(F_press_p * dh)  
    E_press_e = np.cumsum(F_press_e * dh) 
    E_press = E_press_p + E_press_e
    
    

    # Adjustment of the thermal energy to be positive at least of the order of gravitational potential 
    # There not influence on the modeled speed, since it's defined to a constant, however the total energy
    # is positive for existing wind solution
    diff = - E_grav[-1]
    E_press = E_press + diff + abs(E_press[-1])

    
    
    
    
    # Calcul par equation de Bernouilli du chauffage hors cas adiabatique
    
    F_press_adiab_p = np.zeros_like(E_press)
    F_press_adiab_e = np.zeros_like(E_press)
    Pp = k * n_h *1e6 * Tp
    Pe = k * n_h *1e6 * Te
    
    
    i = 0
    F_press_adiab_p[i] = Pp[0]/( n_h[0]**(5/3) ) * ( n_h[i+1]**(5/3) - n_h[i]**(5/3)  ) / (r[i+1] - r[i]) * 1/(1e6*n_h[i])
    F_press_adiab_e[i] = Pe[0]/( n_h[0]**(5/3) ) * ( n_h[i+1]**(5/3) - n_h[i]**(5/3)  ) / (r[i+1] - r[i]) * 1/(1e6*n_h[i])
    
    for i in range(1, N-1):
           
        F_press_adiab_p[i] = Pp[0]/( n_h[0]**(5/3) ) * ( n_h[i+1]**(5/3) - n_h[i-1]**(5/3)  ) / (r[i+1] - r[i-1]) * 1/(1e6*n_h[i])
        F_press_adiab_e[i] = Pe[0]/( n_h[0]**(5/3) ) * ( n_h[i+1]**(5/3) - n_h[i-1]**(5/3)  ) / (r[i+1] - r[i-1]) * 1/(1e6*n_h[i])
    
    F_press_adiab_p[-1] = Pp[0]/( n_h[0]**(5/3) ) * ( n_h[i]**(5/3) - n_h[i-1]**(5/3)  ) / (r[i] - r[i-1]) * 1/(1e6*n_h[i])
    F_press_adiab_e[-1] = Pe[0]/( n_h[0]**(5/3) ) * ( n_h[i]**(5/3) - n_h[i-1]**(5/3)  ) / (r[i] - r[i-1]) * 1/(1e6*n_h[i])
    
    E_press_adiab_p = np.cumsum(F_press_adiab_p * dh)  # k * 5/2 * Tp + k * 5/2 * Tp
    E_press_adiab_e = np.cumsum(F_press_adiab_e * dh)
    
    
    
    # Reajustement des energies par rapport aux valeurs modeliser limite (offset a cause de densite)
    
    E_press_adiab_max_p = k * 5/2 * Tp[-1] 
    E_press_adiab_max_e = k * 5/2 * Te[-1]
    
    diff_p_adiab = E_press_adiab_p[-1] - E_press_adiab_max_p
    diff_e_adiab = E_press_adiab_e[-1] - E_press_adiab_max_e
    
    E_press_adiab_p = E_press_adiab_p - diff_p_adiab
    E_press_adiab_e = E_press_adiab_e - diff_e_adiab
    
    E_press_adiab = E_press_adiab_p + E_press_adiab_e
    #E_press_adiab = E_press_adiab - E_press_adiab[-1] + E_press[-1] 
    
    
    E_cin = E_cin /(1.6e-19*1e3)
    E_press = E_press /(1.6e-19*1e3) #+ 4 
    E_press_p = E_press_p /(1.6e-19*1e3) #+ 4 
    E_press_e = E_press_e /(1.6e-19*1e3) #+ 4 
    E_press_adiab = E_press_adiab /(1.6e-19*1e3) #+ 4
    E_press_adiab_p = E_press_adiab_p /(1.6e-19*1e3) #+ 4 
    E_press_adiab_e = E_press_adiab_e /(1.6e-19*1e3) #+ 4  
    E_grav = E_grav /(1.6e-19*1e3)
    #E_tot = E_tot /(1.6e-19*1e3)
    
    E_tot = E_cin + E_press + E_grav
    delta_E = E_press - E_press_adiab
    
    r = r / r0
    
    
    xvalues = np.array([1, 10, 100, 200])
    
    if(plot_energy):
        plt.figure( figsize=(7,5))
        plt.plot(r, E_cin, linewidth=2 , color='blue', label = '$E_c$')
        plt.plot(r, E_grav, linewidth=2 , color='green', label = '$E_g$')
        plt.plot(r, E_press, linewidth=2, color='red' , label = '$E_{th}$')
        #plt.plot(r, E_press_adiab, ':', linewidth=2, color='gray' , label = '$E_{th|5/3}$')
        #plt.plot(r, E_press_adiab_p , label = 'E_press_adiab_p')
        #plt.plot(r, E_press_adiab_e , label = 'E_press_adiab_e')
        #plt.plot(r, E_press - E_press_adiab,'--' , linewidth=2, color='red' , label = '$\Delta$E')
        plt.plot(r, E_tot, '--', linewidth=2, color='black' , label = '$E_{tot}$')
        plt.xlim([ r[0], 1.05*r[-1] ])
        plt.xscale('log')
        plt.ylim([1.15*np.min([E_cin, E_grav, E_press, E_tot]), 1.07*np.max([E_cin, E_grav, E_press, E_tot]) ])
        #plt.title('E_mec')
        plt.xlabel('Radius ($r \: / \: r_\\odot$)', fontsize=0.9*pol)
        plt.ylabel('Energy (KeV)', fontsize=0.9*pol)
        plt.xticks(fontsize= 0.9*pol)
        plt.yticks(fontsize= 0.9*pol)
        plt.grid()
        plt.legend(loc = 1, fontsize= 0.9*pol)
        ax = plt.gca()
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        plt.axvspan(1, 2.5, facecolor='silver', alpha=0.5)
        
    
   
    return(r, n_h, u_h, Tp, Te, gamma_p, gamma_e, f, bol_supersonic)
        
        
        
        
        
        
