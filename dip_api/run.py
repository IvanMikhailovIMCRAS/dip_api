import numpy as np
import matplotlib.pyplot as plt
import os
from typing import List
from penetration import PenetrationAnalisys

def dip_frame(parameters: List[int|float]) -> np.ndarray:
    """_summary_

    Args:
        parameters (List[int | float]): 
            D - distance between grafting surfaces
            N1 - polymerization degree of the left brush
            N2 - polymerization degree of the right brush
            sigma1 - grafting density of the left brush  
            sigma2 - grafting density of the right brush
            p1 - the Kuhn leght of the left brush polymer chain
            p2 - the Kuhn leght of the right brush polymer chain
            ps - the Kuhn leght of the solvent molecule   
            chi1 - Flory-Huggins parameter for the left polymer-solvent
            chi2 - Flory-Huggins parameter for the right polymer-solvent
            tau1 - dipole moment for the left polymer chain segment
            tau2 - dipole moment for the right polymer chain segment 
            tau_s - dipole moment for the solvent molecule
            Ns - polymerization degree of the solvent molecule
            eta - learning rate
            nfree - number of the free steps 

    Raises:
        TimeoutError: _description_
        FileNotFoundError: _description_

    Returns:
        np.ndarray: _description_
    """
    with open(file='INPUT', mode='w') as file: 
        file.write('\n'.join(list(map(str, parameters))))  
    os.system("./dip.exe") 
    with open(file="INFO", mode='r') as file:
        text = file.read() 
    if "Awesome" not in text:
        raise TimeoutError("dip.exe shows an error")
    if not os.path.exists("profile.out"):
        raise FileNotFoundError("file profile.out is not found")
    pro = np.loadtxt("profile.out", skiprows=1).T
    return pro
        
        
if __name__ == '__main__':
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
    end1 = pro[4]
    end2 = pro[5]
    
    # PA = PenetrationAnalisys(z, phi1, phi2)
    # delta = 0.5 * (PA.delta1 + PA.delta2)
    # z_new = np.linspace(0, D, 1000)
    # phi_A = 0.5 * PA.phi_m * (1.0 - np.tanh((z_new - PA.z_m)/delta))
    # phi_B = 0.5 * PA.phi_m * (1.0 + np.tanh((z_new - PA.z_m)/delta))
   
    # print((PA.delta1 - PA.delta2)/ delta)
    # print(2* (PA.Sigma1 - PA.Sigma2)/ (PA.Sigma1 + PA.Sigma2))
    
    # plt.plot(z, phi1, color='red')
    # plt.plot(z_new, phi_A, '--', color='red')
    # plt.plot(z, phi2, color='blue')
    # plt.plot(z_new, phi_B, '--', color='blue')  
    # plt.ylabel(r'$\varphi_{1,2}(z)$')
    
    # plt.plot(z, phi1*phi2)
    
    plt.plot(z, end1, color='red') 
    plt.plot(z, end2, color='blue')   
    # plt.ylabel(r'$u_{1,2}(z)$')
    
    # plt.xlim(15, 25)    
    plt.xlabel('$z$')
    plt.ylabel('$g(z)$')
    plt.show()
    
    # plt.xlabel('$z$')
    # plt.ylabel(r'$U_{1,2}(z)$')
    # plt.xlim(0, np.max(z)+0.5)
    # plt.plot(z, u1)
    # plt.plot(z, u2)
    # plt.show()