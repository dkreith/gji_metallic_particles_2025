import numpy as np
import cmath
import scipy.optimize as opt
import scipy.interpolate as interp
import scipy.signal as sgnl
import matplotlib.pyplot as plt

def mg(f,sig0,nu):
    sig = (1+2*nu*f)/(1-nu*f)
    return sig

def calc_nu(p1,p2,pk,R,a,L):
    V_K = 4/3*np.pi*a**3
    V_Z = np.pi*(p2-p1)*R**2
    if -p1 < a:
        V_Sl = np.pi/3*(a+p1)**2*(2*a-p1)
    else:
        V_Sl = 0
    if p2 < a:
        V_Sr = np.pi/3*(a-p2)**2*(2*a+p2)
    else:
        V_Sr = 0
    V_eff = V_K-V_Sl-V_Sr
    nu = V_eff/V_Z
    return nu

def reg(nu,sig,sig0,nu0):
    f = 1/nu*(sig-1)/(sig+2)
    sigr = mg(f,sig0,nu0)
    return sigr

def gf(d,R):
    k = np.pi*R**2/d
    return k

def load_data(folder, filename, p1, p2, pk, R, a, L, mod="none", nu0=0.01,
              ind=1):
    
    k = gf((p2-p1),R)
    nu = calc_nu(p1,p2,pk,R,a,L)
    
    dat_U = np.loadtxt(folder+"/"+filename+"_U.txt", dtype=complex)
    dat_I = np.loadtxt(folder+"/"+filename+"_I.txt", dtype=complex)
    
    w = dat_U[:,0]
    U = dat_U[:,ind]
    I = dat_I[:,ind]

    c0 = 1
    mu = 5e-8
    F = 96485.332
    sigma0 = 2*c0*mu*F
    sigma_eff = I/(k*U)/sigma0
    
    if mod == "reg":
        sigma_eff_r = reg(nu, sigma_eff, sigma0, nu0)
    elif mod == "none":
        sigma_eff_r = sigma_eff
    
    return [w.real, sigma_eff_r, filename, sigma0, nu]

def chargeability(sig):
    sig0 = sig.real[0]
    sig1 = sig.real[len(sig)-1]
    m = (sig1-sig0)/sig1
    return m

def wong(omega, c0, c2, T, epsr, mu, a, alpha, beta, nu, mod="norm"):
    sigma_n = []
    
    # Constants
    eps0 = 8.85e-12
    F = 96485
    kB = 8.6173e-5
    
    sig0 = 2*F*mu*c0
    
    for w in omega:
    
        # Define parameters
        D = mu*kB*T
        k = np.sqrt(2*c0*F*mu/eps0/epsr/D)
    
        lambda1 = cmath.sqrt(k**2+1j*w/D)
        lambda2 = cmath.sqrt(1j*w/D)
        
        f1 = (lambda1**2*a**2+2*lambda1*a+2)/((lambda1*a +1)*k**2)*(1j*w/D)
        f2 = (lambda1**2*a**2+2*lambda1*a+2)/(lambda1*a +1)
        f3 = (lambda2*a+1)/(lambda2**2*a**2+2*lambda2*a+2)
    
        f = (1 + ((3*(1+(beta*a*f3/D))+3*c2/(c2-2*c0)*((alpha/mu)-1))/
                   (c2/(c2-2*c0)*(f1+(alpha/mu*(f2-2))+
                                (beta*a*lambda1**2/D/k**2)+2)
                   -(2+f1)*(1+(beta*a/D*f3)))) )
    
        sigma_n.append((1+2*nu*f)/(1-nu*f))
    
    if mod=="norm":
        return np.array(sigma_n)
    elif mod=="sig":
        return np.array(sigma_n)*sig0

def coleColePar(param, w):
    rho = param[0]*(1 - param[1]*(1 - 1/(1+(1j*w*param[2])**param[3])))
    sigma = 1/rho
    return sigma

def findmax(w, sigma, high=0):
    
    # Interpolation
    w_int = np.logspace(np.log10(np.amin(w)),
                        np.log10(np.amax(w)),1000)
    sig_int = interp.interp1d(w, sigma, kind='cubic', fill_value='extrapolate')
    
    sig_int_arr = sig_int(w_int)
    
    peaks, _ = sgnl.find_peaks(sig_int_arr)
    if len(peaks)==0:
        sig_max = np.nan
        tau = np.nan
    else:
        if high == 0:
            index = peaks[0]
        elif high == 1:
            index = peaks[len(peaks)-1]
        w_max = w_int[index]
        tau = 1/w_max
    
        sig_max = sig_int_arr[index]
    
    return [sig_max, tau]

def coleColeFit(w, sig_dat):
    
    def residual(param):
        sig_cc = coleColePar(param, w)
        return (((sig_dat.real-sig_cc.real))+((sig_dat.imag-sig_cc.imag)))
    
    rho0_est = np.real(1/sig_dat[0])
    m_est = chargeability(sig_dat)
    sig_max, tau_est = findmax(w, sig_dat.imag)
    
    sol = opt.least_squares(residual, [rho0_est, m_est, tau_est, 1])
    
    return sol.x

def fit_c3(c0, T, epsr, mu, a, alpha, beta, nu, w_dat, sig_dat,
           mod=None):
 
    sigmax_dat, tau_dat = findmax(w_dat, sig_dat.imag)
    w0 = np.logspace(-8, 8, 1600)
    c3 = np.logspace(-10, np.log10(c0), 100)
    sigmax = []
    
    for c3_temp in c3:
        sig_temp = wong(w0, c0, c3_temp, T, epsr, mu, a, alpha, beta, nu,
                        mod="sig")
        sigmax_temp, tau_temp = findmax(w0, sig_temp.imag)
        sigmax.append(sigmax_temp)
    
    c3_int = np.logspace(-10, np.log10(c0),5000)
    sigmax_int = interp.interp1d(c3, sigmax, kind='cubic',
                             fill_value='extrapolate')

    sigmax_int_arr = sigmax_int(c3_int)
    
    c3_opt = 0
    
    for ii in range(len(c3_int)-1):
        if (sigmax_int_arr[ii] >= sigmax_dat and 
            sigmax_int_arr[ii+1] < sigmax_dat):
            c3_opt = c3_int[ii]
        
    if mod == "plot":
        fig, ax = plt.subplots()
        ax.semilogx(c3_int, sigmax_int_arr)
        ax.plot([1e-10, c3_opt], [sigmax_dat, sigmax_dat], ":k")
        ax.plot([c3_opt, c3_opt], [sigmax_int_arr[-1], sigmax_dat], ":k")
    
    return c3_opt

def wong_eff(omega, c0, c2, T, epsr, mu, a, alpha, beta, nu, w_dat, sig_dat,
             mod="sig"):
    
    _, tau_dat = findmax(w_dat, sig_dat.imag)
    
    wong_0 = wong(omega, c0, c2, T, epsr, mu, a, alpha, beta, nu, mod=mod)
    _, tau_0 = findmax(omega, wong_0.imag)
    
    mu_fac = tau_0/tau_dat

    wong_1 = wong(omega, c0, c2, T, epsr, mu*mu_fac, a, alpha, beta, nu,
                  mod=mod)/mu_fac
    return wong_0, wong_1