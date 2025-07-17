import numpy as np
import os
from typing import List

# HOME_DIR = f"{os.path.dirname(os.path.realpath(__file__))}"[:-7]
# print(HOME_DIR)
HOME_DIR = os.getcwd()

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
    with open(file=os.path.join(HOME_DIR,"INPUT"), mode='w') as file:
        file.write('\n'.join(list(map(str, parameters))))  
    os.system(os.path.join(HOME_DIR,"dip.exe")) 
    with open(file=os.path.join(HOME_DIR,"INFO"), mode='r') as file:
        text = file.read() 
    if "Awesome" not in text:
        raise TimeoutError("dip.exe shows an error")
    if not os.path.exists(os.path.join(HOME_DIR,"profile.out")):
        raise FileNotFoundError("file profile.out is not found")
    pro = np.loadtxt(os.path.join(HOME_DIR,"profile.out"), skiprows=1).T
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