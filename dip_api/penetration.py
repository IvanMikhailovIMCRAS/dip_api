import warnings

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import brentq, fsolve

warnings.filterwarnings("ignore")


class PenetrationAnalisys:

    def __init__(self, z: np.ndarray, phi1: np.ndarray, phi2: np.ndarray) -> None:
        """class initialization

        Args:
            z (np.ndarray): distance from left grafting surface
            phi1 (np.ndarray): volume fraction of left brush
            phi2 (np.ndarray): volume fraction of right brush
        """
        self.z = np.hstack([0.0, z, z[-1] + z[0]])

        self.D = np.max(self.z)

        self.phi1 = interp1d(
            self.z, np.hstack([1e-9, phi1, phi1[-1]]), kind="quadratic"
        )

        self.phi2 = interp1d(self.z, np.hstack([phi2[0], phi2, 1e-9]), kind="quadratic")

        self.z_m = brentq(
            lambda x: self.phi1(x) - self.phi2(x), z.min(), z.max(), xtol=1e-12
        )
        self.phi_m = self.phi1(self.z_m) + self.phi2(self.z_m)

        self.Sigma1, _ = quad(self.phi1, self.z_m, self.D)
        self.Sigma2, _ = quad(self.phi2, 0.0, self.z_m)

        self.delta1 = fsolve(self._target1, 1.0)[0]
        self.delta2 = fsolve(self._target2, 1.0)[0]

        self.friction_force = np.sum(phi1**2 * phi2**2 / (phi1**2 + phi2**2))

    @property
    def Pressure(self) -> float:
        """Calculation of disjoining pressure

        Args:
            phi_m (float): value of midpoint volume fraction at equilibrium point of the system

        Returns:
            float: pressure value
        """
        return -self.phi_m - np.log(1 - self.phi_m)

    @property
    def Gamma1(self) -> float:
        """Overlap integral

        Returns:
            float: left brush overlap integral
        """
        return quad(lambda z: self.phi1(z) * self.phi2(z), self.z_m, self.D)[0]

    @property
    def Gamma2(self) -> float:
        """Overlap integral

        Returns:
            float: right brush overlap integral
        """
        return quad(lambda z: self.phi1(z) * self.phi2(z), 0.0, self.z_m)[0]

    @property
    def L1(self) -> float:
        """first moment (Leermakers)

        Returns:
            float: left brush inperpenetration zone
        """
        numerator, _ = quad(lambda z: self.phi1(z) * (z - self.z_m), self.z_m, self.D)
        return numerator / self.Sigma1

    @property
    def L2(self) -> float:
        """first moment (Leermakers)

        Returns:
            float: right brush inperpenetration zone
        """
        numerator, _ = quad(lambda z: self.phi2(z) * (self.z_m - z), 0.0, self.z_m)
        return numerator / self.Sigma2

    def _target1(self, delta):
        norm, _ = quad(
            lambda z: self.phi_m / 2 * (1.0 - np.tanh((z - self.z_m) / delta)),
            self.z_m,
            self.D,
        )
        return (self.Sigma1 - norm) ** 2

    def _target2(self, delta):
        norm, _ = quad(
            lambda z: self.phi_m / 2 * (1.0 + np.tanh((z - self.z_m) / delta)),
            0.0,
            self.z_m,
        )
        return (self.Sigma2 - norm) ** 2

    def rho1(self, z):
        return self.phi_m / 2 * (1.0 - np.tanh((z - self.z_m) / self.delta1))

    def rho2(self, z):
        return self.phi_m / 2 * (1.0 + np.tanh((z - self.z_m) / self.delta2))
