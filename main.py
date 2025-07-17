from dip_api import dip_frame
import matplotlib.pyplot as plt
from dip_api import PenetrationAnalisys
   
D = 50
N1 = 210 
N2 = 210 
sigma1 = 0.1   
sigma2 = 0.1 
p1 = 1.0  
p2 = 10.0
ps = 1.0   
chi1 = 0.0   
chi2 = 0.0  
tau1 = 0.0  
tau2 = 0.0  
tau_s= 0.0   
Ns = 1    
eta = 0.02   
nfree = 1000    
p = [D, N1, N2, sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, tau_s, Ns, eta, nfree]
pro = dip_frame(p)

z = pro[0] - 0.5
phi1 = pro[2]
phi2 = pro[3]

PA = PenetrationAnalisys(z, phi1, phi2)

print(PA.z_m)

plt.plot(z, phi1, color="red")
plt.plot(z, phi2, color="blue")
plt.ylabel(r"$\varphi(z)$")
plt.xlabel("$z$")
plt.show()