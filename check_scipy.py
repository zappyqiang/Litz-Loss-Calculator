"""Cross-check of core.js: skin factor and strand proximity loss vs SciPy and direct integration.
Compare with: node -e \"const {strandAC}=require('./core.js');console.log(strandAC(0.1e-3,1e6,1.724e-8))\"
Columns: f, F_skin, k_p [W/m per (A/m rms)^2], proximity / low-frequency limit."""
import numpy as np
from scipy.special import iv
from scipy.integrate import quad
mu0=4e-7*np.pi; rho=1.724e-8; sig=1/rho; d=0.1e-3; r=d/2
for f in [1e5,1e6,5e6]:
    w=2*np.pi*f; delta=np.sqrt(2*rho/(w*mu0)); g=(1+1j)/delta; z=g*r
    FR=(z*iv(0,z)/(2*iv(1,z))).real
    # proximity by direct integration of |J|^2/(2 sigma), H0=1 A/m peak
    C=2*mu0/(g*iv(0,z))
    integ=quad(lambda p: abs(iv(1,g*p))**2*p,0,r,limit=200)[0]
    Pp=(w**2*sig*abs(C)**2*np.pi/2)*integ
    Plf=sig*np.pi*w**2*mu0**2*r**4/8
    print(f, FR, 2*Pp, Pp/Plf)
