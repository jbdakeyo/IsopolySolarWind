# Package used is the code
import math
import numpy as np
import sys
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

###############################################################
###### Main function calculating the solar wind solution ######
###############################################################
def function_to_run_isopoly(r0, Te0, Tp0, gamma_p_max, gamma_e_max, ind_r_iso_p, ind_r_iso_e,
                                                                L ,N, f, bol_iso):
    
    
    
    if((gamma_p_max <= 0) | (gamma_e_max <= 0)): 
        print('-------------------------------------------------------')
        print('ERROR : gamma_p or gamma_e <= 0 --> Undefined thermal regime')
        print('-------------------------------------------------------')
        sys.exit()
    elif((Tp0 < 0) | (Te0 < 0)): 
        print('-------------------------------------------------------')
        print('ERROR : Tp0 or Te0 < 0 --> Unphysical quantities')
        print('-------------------------------------------------------')
        sys.exit()
    elif((ind_r_iso_p == 0) | (ind_r_iso_e == 0)): 
        print('-------------------------------------------------------')
        print('EROR : r_iso_p or r_iso_e <= 0')
        print('--> Misdefined quantities or Undefined solution')
        print('-------------------------------------------------------')
        sys.exit()
        
    if(bol_iso != 1):
        if((ind_r_iso_p >= N-1) | (ind_r_iso_e >= N-1)): 
            print('-------------------------------------------------------')
            print('WARNING : r_iso_p or r_iso_e > the size of the domain')
            print('--> Decrease r_iso_p or r_iso_e')
            print('-------------------------------------------------------')
            
        
    
        

    # Centered derivative function
    def derive_cen(t,x):
        N = len(x)
        dxdt = np.zeros_like(x)
        for i in range(1, N-1):
            dxdt[i] = (x[i+1] - x[i-1]) / (t[i+1] - t[i-1])
            dxdt[0] = dxdt[1]
            dxdt[-1] = dxdt[-2]
            
        return dxdt
    
    # Differential function
    def differential(x):
        N = len(x)
        dxdt = np.zeros_like(x)
        for i in range(N-1):
            dxdt[i] = (x[i+1] - x[i])
            dxdt[-1] = dxdt[-2]            
        return dxdt
    
    # Local sound speed
    cs_T_back = []
    cs_T_for = []
    
    ### Percentage around rc and uc where the equations are linearized to treat the sonic point calculation
    lim_pourc = 0.001 
    
    
    ################################################
    ### Forward calculation in finite difference ###
    ################################################
    def grad_forward_polytr(r, u_h, ind_rc, gamma_p, gamma_e, N, cs_poly_p, cs_poly_e) : 

        for i in range(ind_rc, N-1):
            deriv_csp2 = ( cs_poly_p[i]**2 - cs_poly_p[i-1]**2 ) / (r[i] - r[i-1])
            deriv_cse2 = ( cs_poly_e[i]**2 - cs_poly_e[i-1]**2 ) / (r[i] - r[i-1])
            
            gamma_p_i = gamma_p[i]
            gamma_e_i = gamma_e[i]
            
            Cn = cs_iso * rc_iso**2
            
            ########################################################################
            # Separate the initialization of quantities in isothermal and polytropic 
            # region for electrons and protons
            
            if(gamma_p_i !=1):
                if(u_h[ind_r_iso_p] == 0): break
                x_p = ( (u_h[ind_r_iso_p] * r[ind_r_iso_p]**2 * f[ind_r_iso_p])/(u_h[i] * r[i]**2 * f[i]) )**( gamma_p_i - 1) 
                n_tilde_iso_p = Cn / (u_h[ind_r_iso_p] * r[ind_r_iso_p]**2)
            else : 
                x_p = 1.
                # Value below is set to 1 to avoid division by zero, but since they are multiplied by zero there is no influence 
                n_tilde_iso_p = 1.
                
                
            if(gamma_e_i !=1): 
                if(u_h[ind_r_iso_e] == 0): break
                x_e = ( (u_h[ind_r_iso_e] * r[ind_r_iso_e]**2 * f[ind_r_iso_e])/(u_h[i] * r[i]**2 * f[i]) )**( gamma_e_i - 1 ) 
                n_tilde_iso_e = Cn / (u_h[ind_r_iso_e] * r[ind_r_iso_e]**2)
            else : 
                x_e = 1.
                # Values below is set to 1 to avoid division by zero, but since they are multiplied by zero there is no influence 
                n_tilde_iso_e = 1.
            ########################################################################
         
            
            # Bi-fluid sound speed 
            u_c = np.sqrt( cs_poly_p[i]**2 * x_p + cs_poly_e[i]**2 * x_e )
 
            # Treatment of the sonic point calculation forward+
            if( abs(r[ind_rc] - r[i])/r[i] <= lim_pourc or abs(u_c - u_h[i])/u_h[i] <= lim_pourc  ): 
                signe = +1
                pourc = (100 + signe*lim_pourc)/100
                alpha = ( 2 *  u_h[i]/r[i] ) * (1 - 1/pourc) / ( pourc**2 - 1 )
   
            # Resolution far outside the neighborhood of the sonic point
            else:
                AA = 1 - ( cs_poly_p[i]**2 * x_p + cs_poly_e[i]**2 * x_e ) / u_h[i]**2   
                BB = 2 * ( 1 - 0.5*f[i]*r[i] * d_inv_f[i] ) * ( cs_poly_p[i]**2 * x_p + cs_poly_e[i]**2 * x_e ) - G*M/r[i]
                if( (n_tilde_iso_p ==1) | (n_tilde_iso_e ==1) ):
                    CC = 0
                else:
                    CC =  - ( x_p * deriv_csp2 * np.log( 1./n_tilde_iso_p * Cn/(u_h[i]  * r[i]**2) ) * ( 1. / u_h[i])  \
                        + x_e * deriv_cse2 * np.log( 1./n_tilde_iso_e * Cn/(u_h[i]  * r[i]**2) ) * ( 1. / u_h[i])  ) 
                alpha = ( BB/(u_h[i] * r[i]) + CC) / AA 
    
    
            # Speed forward
            u_h[i+1] = u_h[i] + (r[i+1] - r[i]) * alpha
            
            # Local sound speed calculation
            Tp_i = Tp0 * x_p 
            Te_i = Te0 * x_e 
            cs_p_i = np.sqrt( k * gamma_p_i * Tp_i /mp )
            cs_e_i = np.sqrt( k * gamma_e_i * Te_i /mp )
            #print(cs_p_i)
            cs_T_for.append( np.sqrt(cs_p_i**2 *x_p + cs_e_i**2 *x_e) )
            
            
    #################################################
    ### Backward calculation in finite difference ###
    #################################################
    def grad_backward_polytr(r, u_h, ind_rc, gamma_p, gamma_e, cs_poly_p, cs_poly_e) : 
    
        for i in range(ind_rc, 0, -1):
            deriv_csp2 = ( cs_poly_p[i]**2 - cs_poly_p[i-1]**2 ) / (r[i] - r[i-1])
            deriv_cse2 = ( cs_poly_e[i]**2 - cs_poly_e[i-1]**2 ) / (r[i] - r[i-1])
            
            gamma_p_i = gamma_p[i]
            gamma_e_i = gamma_e[i]
            
            Cn = cs_iso * rc_iso**2

            ########################################################################
            # Separate the initialization of quantities in isothermal and polytropic 
            # region for electrons and protons
            
            if(gamma_p_i !=1):
                if(u_h[ind_r_iso_p] == 0): break
                x_p = ( (u_h[ind_r_iso_p] * r[ind_r_iso_p]**2 * f[ind_r_iso_p])/(u_h[i] * r[i]**2 * f[i]) )**( gamma_p_i - 1) 
                n_tilde_iso_p = Cn / (u_h[ind_r_iso_p] * r[ind_r_iso_p]**2)
            else : 
                x_p = 1.
                # Value below is set to 1 to avoid division by zero, but since they are multiplied by zero there is no influence 
                n_tilde_iso_p = 1.
                
                
            if(gamma_e_i !=1): 
                if(u_h[ind_r_iso_e] == 0): break
                x_e = ( (u_h[ind_r_iso_e] * r[ind_r_iso_e]**2 * f[ind_r_iso_e])/(u_h[i] * r[i]**2 * f[i]) )**( gamma_e_i - 1 ) 
                n_tilde_iso_e = Cn / (u_h[ind_r_iso_e] * r[ind_r_iso_e]**2)
            else : 
                x_e = 1.
                # Value below is set to 1 to avoid division by zero, but since they are multiplied by zero there is no influence 
                n_tilde_iso_e = 1.


            # Bi-fluid sound speed 
            u_c = np.sqrt( cs_poly_p[i]**2 * x_p + cs_poly_e[i]**2 * x_e )
  
 
            # Treatment of the sonic point calculation backward
            if( abs(r[ind_rc] - r[i])/r[i] <= lim_pourc or abs(u_c - u_h[i])/u_h[i] <= lim_pourc  ):
                signe = -1
                pourc = (100 + signe*lim_pourc)/100
                alpha = ( 2 *  u_h[i]/r[i] ) * (1 - 1/pourc) / ( pourc**2 - 1 )
            # Resolution far outside the neighborhood of the sonic point
            else:
                AA = 1 - ( cs_poly_p[i]**2 * x_p + cs_poly_e[i]**2 * x_e ) / u_h[i]**2   
                BB = 2 * ( 1 - 0.5*f[i]*r[i] * d_inv_f[i] ) * ( cs_poly_p[i]**2 * x_p + cs_poly_e[i]**2 * x_e ) - G*M/r[i]
                if( (n_tilde_iso_p ==1) | (n_tilde_iso_e ==1) ):
                    CC = 0
                else:
                    CC =  - ( x_p * deriv_csp2 * np.log( 1./n_tilde_iso_p * Cn/(u_h[i]  * r[i]**2) ) * ( 1. / u_h[i])  \
                        + x_e * deriv_cse2 * np.log( 1./n_tilde_iso_e * Cn/(u_h[i]  * r[i]**2) ) * ( 1. / u_h[i])  ) 
                alpha = ( BB/(u_h[i] * r[i]) + CC) / AA 
    
    
            # Speed backward
            u_h[i-1] = u_h[i] - (r[i] - r[i-1]) * alpha
            
            # Local sound speed calculation
            Tp_i = Tp0 * x_p 
            Te_i = Te0 * x_e 
            cs_p_i = np.sqrt( k * gamma_p_i * Tp_i /mp )
            cs_e_i = np.sqrt( k * gamma_e_i * Te_i /mp )
            cs_T_back.append( np.sqrt(cs_p_i**2 *x_p + cs_e_i**2 *x_e) )

    # Initialization of physical quantities
    mp = 1.67e-27
    M = 1.99e30
    G = 6.67e-11
    k = 1.38e-23
    r0 = 6.96e8

    
    # Sonic point from coronal input temperature
    cs_iso = np.sqrt( k * (Te0 + Tp0) / mp )
    rc_iso = G*M / (2 * cs_iso**2)
    
    
    if(cs_iso == 0): 
        print('ERROR -- Invalid input coronal temperatures')
        return
    
    
    # The radial distance vector is defined logarithmically to reduce time calculation at large distance (slow quantities variation)
    r = np.geomspace(r0, L, N) 
    r_iso_p = r[ind_r_iso_p]
    r_iso_e = r[ind_r_iso_e]

    
    ###########################################################
    # Numerical search of rc, assuming rc is located in the isothermal region (and so xs=1)
    d_inv_f = derive_cen(r, 1/f)
    
    # Function to minimize to find rc
    func_rc = cs_iso**2 * (2 - r*f * d_inv_f) - G*M/r
    # Derivative of func_rc, used for the oversampling
    func_rc_prim = derive_cen(r, func_rc)
    
    
    # Re-sampling of r vector depending super radial expansion derivative (to better treat large gradient with finite difference)
    # Larger variation --> smaller dh
    ######################################################################
    dfdr = derive_cen(r, f)
    derive_lim_f = 0.05/(r[1]-r[0])  # Value of df/dr over which we oversample
    ratio_f = dfdr/derive_lim_f
    derive_lim_func = 100           # Value of func_rc over which we oversample
    ratio_func = abs(func_rc_prim)/derive_lim_func
    ind_fprim = np.argwhere( (ratio_func > 1) | (ratio_f >1) )[:,0] 
    r_new = []

    # Oversampling of r vector 
    for pp in ind_fprim:
        nb_pts = int(ratio_f[pp]) + int(ratio_func[pp])
        for mm in range(1, nb_pts):
            r_new.append(  (1-mm/nb_pts)*r[pp] + (mm/nb_pts)*r[pp+1]  )
    r_conc = np.concatenate((r, np.array(r_new) ),axis=-1)
    r_conc.sort()
    
    r_old = r.copy()
    f_old = f.copy()
    # Interpolation of f to match with r vector
    func_interp_fexp = interp1d( r , f, kind='linear',axis=0)
    f = func_interp_fexp(r_conc)
    r = r_conc.copy()
    
    ind_r_iso_p = np.argmin( abs(r - r_iso_p) )
    ind_r_iso_e = np.argmin( abs(r - r_iso_e) )
    
    
    
    #################################################
    # Illustrating the localized oversampling
    '''
    pol = 14
    xvalues = np.array([1, 10, 100, 200])
    dh_old = differential(r_old)
    dh_opt = differential(r)

    plt.figure()
    plt.plot(r_old/r0, dh_old, label='$dr_{(init)}$')
    plt.plot(r/r0, dh_opt, label='$dr_{(optim)}$')
    plt.xlabel('Radius ($r \: / \: r_0$)', fontsize=pol)
    plt.ylabel('Spatial step of the resolution', fontsize=pol)
    plt.xscale('log')
    plt.yscale('log')
    ax = plt.gca()
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    plt.xticks(xvalues, fontsize= pol)
    plt.yticks(fontsize= pol)
    plt.grid()
    plt.legend(fontsize= 0.9*pol)
    '''
    #################################################
    
    
    

    ###########################################################
    # Numerical search of rc --> UPDATED with oversampled r vector
    d_inv_f = derive_cen(r, 1/f)
    func_rc = cs_iso**2 * (2 - r*f * d_inv_f) - G*M/r
    peaks, _ = find_peaks( - abs(func_rc))
    ind_rc = peaks[-1]

    
    # Case of rc associated with f-supersonic type solution
    ind_rc_super=0
    if( any(func_rc > 0)):
        ind_plus = np.argwhere(func_rc > 0)[0,0]
        if(ind_plus > 0 ):
            dist = np.array([ func_rc[ind_plus] , func_rc[ind_plus-1] ])
            ind_rc_super = ind_plus - np.argmin(abs(dist))
    else: ind_rc_super = 0
    ind_rc_new = np.array([ind_rc_super, ind_rc])
    
    
    



    
    ###########################################################
    # Redefining r_iso for oversampled vector
    r_iso_p = r[ind_r_iso_p]
    r_iso_e = r[ind_r_iso_e]
    
    #########################################
    # Builing gamma_p and gamma_e vectors
    gamma_p_0 = 1
    gamma_e_0 = 1
    
    # Gamma isotherme --> gamma = 1
    gamma_iso_p = np.linspace(gamma_p_0, gamma_p_0, ind_r_iso_p  ) 
    gamma_iso_e = np.linspace(gamma_e_0, gamma_e_0, ind_r_iso_e  ) 
    
    # Gamma polytropique --> gamma > 1        
    ind_poly_p = len(r) - ( ind_r_iso_p ) 
    ind_poly_e = len(r) - ( ind_r_iso_e )
    
    gamma_poly_p = np.geomspace(gamma_p_max, gamma_p_max, ind_poly_p  )  
    gamma_poly_e = np.geomspace(gamma_e_max, gamma_e_max, ind_poly_e  )

    gamma_p = np.concatenate( (gamma_iso_p, gamma_poly_p), axis=0)
    gamma_e = np.concatenate( (gamma_iso_e, gamma_poly_e), axis=0)

    cs_poly_p = np.sqrt( gamma_p * k * Tp0 / mp )
    cs_poly_e = np.sqrt( gamma_e * k * Te0 / mp )

    
    # Velocity 
    u_h = np.zeros_like(r)
    u_h[ind_rc] = cs_iso    # Velocity at the soic point


    # Determination of the wind solution
    ############################################################
    bol_resol = False
    bol_supersonic = 0
    cs_T_back.append(cs_iso)
    
    # Isothermal part of the solution (backward)
    grad_backward_polytr(r, u_h, ind_rc, gamma_p, gamma_e, cs_poly_p, cs_poly_e)
    
    # Exception in case of f-subsonic solution is not verified (subosnic before rc)
    # --> f-supersonic solution
    if( (np.max(u_h[:ind_rc]) > cs_iso) | (math.isnan(u_h[0])) | (math.isnan(u_h[-1])) 
           | (r_iso_p <= 1.1*rc_iso) | (r_iso_e <= 1.1*rc_iso) ):
        
        if(ind_rc_new[0] == 0): 
            print('-------------------------------------------------------')
            print('ERROR : rc < r0 --> No transonic solution')
            print('-------------------------------------------------------')
            sys.exit()
            
            
        rc_iso = r[ind_rc_new[0]]  
        ind_rc = ind_rc_new[0]  
        u_h[ind_rc] = cs_iso       

        cs_T_back = []
        cs_T_back.append(cs_iso)
        grad_backward_polytr(r, u_h, ind_rc, gamma_p, gamma_e, cs_poly_p, cs_poly_e)
        '''
        print('-----------    -----------')
        print('SOLUTION TYPE : Supersonic before rc')
        print('-----------    -----------')
        '''
        bol_supersonic = True
        
    # Both isothermal and polytropic part of the solution (forward)
    grad_forward_polytr(r, u_h, ind_rc, gamma_p, gamma_e, len(r), cs_poly_p, cs_poly_e)
    
    # If it fails to compute correctly the forward step 
    if( u_h[ind_r_iso_p] == 0 ): 
        print('------------------------------------------------------------------------')
        print('ERROR : r_iso_p or/and r_iso_e < rc --> No solution')
        print('------------------------------------------------------------------------')
        sys.exit()


    # Local sound speed
    cs_T = np.concatenate(( np.array(cs_T_back[::-1]), np.array(cs_T_for) ))
    

    # Re-sampling tot fit to the initial size of r 
    func_interp_u_h = interp1d( r, u_h, kind='linear', axis=-1)
    func_interp_gp = interp1d( r, gamma_p, kind='linear', axis=-1)
    func_interp_ge = interp1d( r, gamma_e, kind='linear', axis=-1)
    func_interp_cs = interp1d( r, cs_T, kind='linear', axis=-1)
    u_h = func_interp_u_h( r_old ) 
    f = f_old.copy()
    gamma_p = func_interp_gp( r_old )
    gamma_e = func_interp_ge( r_old )
    cs_T = func_interp_cs( r_old )
    
    ind_r_iso_p = np.argmin( abs(r_old - r_iso_p) ) 
    ind_r_iso_e = np.argmin( abs(r_old - r_iso_e) ) 
    ind_rc = np.argmin( abs(r_old - rc_iso) ) 
    r = r_old
    
    
    # If velocity go down cs after the sonic point
    if( any(u_h[ind_rc:] < cs_T[ind_rc:]) ): 
        print('------------------------------------------------------------------------')
        print('ERROR : N too small     OR    (r_iso_p , r_iso_e) < rc  -->  No solution')
        print('--> Increase N, or (Tp0 ,Tpe), or (r_iso_p, r_iso_e)')
        print('------------------------------------------------------------------------')
        sys.exit()
    
    
    
    # Density calculation from mass flux conservation
    n_h = ( u_h[ind_rc] * r[ind_rc]**2 * f[ind_rc] ) / ( u_h * r**2 * f )
    
    # Temperature calculation from polytropic relation
    Tp_riso = Tp0
    Te_riso = Te0
    
    Tp = Tp_riso * ( n_h * 1./n_h[ind_r_iso_p] )**(gamma_p - 1)  
    Te = Te_riso * ( n_h * 1./n_h[ind_r_iso_e] )**(gamma_e - 1) 
    
     
    

    return(r, rc_iso, u_h, n_h, Tp, Te, cs_T, ind_rc, gamma_p, gamma_e, bol_supersonic)












