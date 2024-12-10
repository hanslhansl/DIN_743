import math as m
from abc import abstractmethod
from dataclasses import dataclass
from DIN_743_3 import *



@dataclass
class Kerbe:
    d : float

    def sigma_zd(self, F_zd):
        return F_zd / (m.pi / 4 * self.d**2)
    def sigma_b(self, M_b):
        return M_b / (m.pi / 32 * self.d**3) * 1000
    def tau_t(self, M_t):
        return M_t / (m.pi / 16 * self.d**3) * 1000

    @property
    def D(self):
        return self.d + 2 * self.t
    pass


class K_F_nach_Welle_Nabe:
    """Tabelle 1"""
    def K_Fsigma(self, Rz, sigma_B_d_eff):
        return 1
    def K_Ftau(self, Rz, sigma_B_d_eff):
        return 1
class K_F_nach_Formel:
    """Glg 18 & 19"""
    def _K_Fsigma(self, Rz, sigma_B_d_eff):
        return 1 - 0.22 * m.log(Rz, 10) * (m.log(sigma_B_d_eff / 20, 10) - 1)
    def _K_Ftau(self, Rz, sigma_B_d_eff):
        return 0.575 * self._K_Fsigma(Rz, sigma_B_d_eff) + 0.425

    def K_Fsigma(self, Rz, sigma_B_d_eff):
        if hasattr(self, "Rz_B"):
           return self._K_Fsigma(Rz, sigma_B_d_eff) / self._K_Fsigma(self.Rz_B, sigma_B_d_eff)
        return self._K_Fsigma(Rz, sigma_B_d_eff)
    def K_Ftau(self, Rz, sigma_B_d_eff):
        if hasattr(self, "Rz_B"):
           return self._K_Ftau(Rz, sigma_B_d_eff) / self._K_Ftau(self.Rz_B, sigma_B_d_eff)
        return self._K_Ftau(Rz, sigma_B_d_eff)


class ExperimentelleKerbwirkungszahlen(Kerbe):
    def beta_sigma_zd(self, **kwargs):
        beta_d_BK = self.beta_sigma_zd_d_BK(**kwargs)
        K_3_d = K_3(self.d, beta_d_BK)
        K_3_d_BK = K_3(self.d_BK, beta_d_BK)
        print(f"\tbeta_sigma_zd_d_BK = {beta_d_BK}")
        print(f"\tK_3_zd_d_BK = {K_3_d_BK}")
        print(f"\tK_3_zd_d = {K_3_d}")
        return beta_d_BK * K_3_d_BK / K_3_d
    def beta_sigma_b(self, **kwargs):
        beta_d_BK = self.beta_sigma_b_d_BK(**kwargs)
        K_3_d = K_3(self.d, beta_d_BK) 
        K_3_d_BK = K_3(self.d_BK, beta_d_BK)
        print(f"\tbeta_sigma_b_d_BK = {beta_d_BK}")
        print(f"\tK_3_b_d_BK = {K_3_d_BK}")
        print(f"\tK_3_b_d = {K_3_d}")
        return beta_d_BK * K_3_d_BK / K_3_d
    def beta_tau(self, **kwargs):
        beta_d_BK = self.beta_tau_d_BK(**kwargs)
        K_3_d = K_3(self.d, beta_d_BK)
        K_3_d_BK = K_3(self.d_BK, beta_d_BK)
        print(f"\tbeta_tau_d_BK = {beta_d_BK}")
        print(f"\tK_3_t_d_BK = {K_3_d_BK}")
        print(f"\tK_3_t_d = {K_3_d}")
        return beta_d_BK * K_3_d_BK / K_3_d


@dataclass
class Passfeder(ExperimentelleKerbwirkungszahlen, K_F_nach_Welle_Nabe):
    """Tabelle 1"""
    umdrehungskerbe = False
    d_BK = 40

    i : int

    def __post_init__(self):
        assert self.i in (1, 2)
        
    def beta_sigma_zd_d_BK(self, sigma_B_d: float, **kwargs):
        return (3 * (sigma_B_d / 1000)**0.38) * (1.15 if self.i == 2 else 1)
    def beta_sigma_b_d_BK(self, sigma_B_d: float, **kwargs):
        return self.beta_sigma_zd_d_BK(sigma_B_d)
    def beta_tau_d_BK(self, sigma_B_d: float, **kwargs):
        return (0.56 * 3 * (sigma_B_d / 1000)**0.38 + 0.1) * (1.15 if self.i == 2 else 1)
    
@dataclass
class Presssitz(ExperimentelleKerbwirkungszahlen, K_F_nach_Welle_Nabe):
    """Tabelle 1"""
    umdrehungskerbe = True
    d_BK = 40

    def beta_sigma_zd_d_BK(self, sigma_B_d: float, **kwargs):
        return 2.7 * (sigma_B_d / 1000)**0.43
    def beta_sigma_b_d_BK(self, sigma_B_d: float, **kwargs):
        return self.beta_sigma_zd_d_BK(sigma_B_d)
    def beta_tau_d_BK(self, sigma_B_d: float, **kwargs):
        return 0.65 * self.beta_sigma_d_BK(sigma_B_d)

    def sigma_zd(self, F_zd):
        return F_zd / (m.pi / 4 * self.d**2)
    def sigma_b(self, M_b):
        return M_b / (m.pi / 32 * self.d**3)
    def sigma_t(self, M_t):
        return M_t / (m.pi / 16 * self.d**3)
    
@dataclass
class Spitzkerbe(ExperimentelleKerbwirkungszahlen, K_F_nach_Formel):
    """Abschnitt 4.2.3"""
    umdrehungskerbe = True
    r = 0.1
    d_BK = 15
    Rz_B = 20

    t : float

    def beta_sigma_zd_d_BK(self, sigma_B_d: float, **kwargs):
        return 0.109 * sigma_B_d / 100 + 1.074
    def beta_sigma_b_d_BK(self, sigma_B_d: float, **kwargs):
        return 0.0923 * sigma_B_d / 100 + 0.985
    def beta_tau_d_BK(self, sigma_B_d: float, **kwargs):
        return 0.8 * self.beta_sigma_b_d_BK(sigma_B_d)

    def sigma_zd(self, F_zd):
        assert 0.05 <= self.t/self.d <= 0.20
        return super().sigma_zd(F_zd, self.d)
    def sigma_b(self, M_b):
        assert 0.05 <= self.t/self.d <= 0.20
        return super().sigma_b(M_b, self.d)
    def sigma_t(self, M_t):
        assert 0.05 <= self.t/self.d <= 0.20
        return super().sigma_t(M_t, self.d)
    
@dataclass
class UmlaufendeRechtecknut(ExperimentelleKerbwirkungszahlen, K_F_nach_Formel):
    """Abschnitt 4.2.4"""
    umdrehungskerbe = True
    d_BK = 30
    Rz_B = 20

    t : float
    r : float

    def __post_init__(self):
        raise NotImplementedError

    def rho_s(self, sigma_S_d):
        return 10 ** -(0.514 + 0.00152 - sigma_S_d)

    def beta_sigma_zd_d_BK(self, sigma_S_d: float, **kwargs):
        r_f = self.r + 2.9 * self.rho_s(sigma_S_d)
        return min(0.9 * (1.27 + 1.17 * m.sqrt(self.t / r_f)), 4)
    def beta_sigma_b_d_BK(self, sigma_S_d: float, **kwargs):
        r_f = self.r + 2.9 * self.rho_s(sigma_S_d)
        return min(0.9 * (1.14 + 1.08 * m.sqrt(self.t / r_f)), 4)
    def beta_tau_d_BK(self, sigma_S_d: float, **kwargs):
        r_f = self.r + self.rho_s(sigma_S_d)
        return min(1.48 + 0.45 * m.sqrt(self.t / r_f), 2.5)


class BekannteFormzahl(Kerbe, K_F_nach_Formel):
    """4.3.1"""
    def beta_sigma_zd(self, **kwargs):
        alpha = self.alpha_sigma_zd()
        n = self.n_zd(**kwargs)
        print(f"\talpha_sigma_zd = {alpha}")
        print(f"\tn_zd = {n}")
        return alpha / n
    def beta_sigma_b(self, **kwargs):
        alpha = self.alpha_sigma_b()
        n = self.n_b(**kwargs)
        print(f"\talpha_sigma_b = {alpha}")
        print(f"\tn_b = {n}")
        return alpha / n
    def beta_tau(self, **kwargs):
        alpha = self.alpha_tau()
        n = self.n_t(**kwargs)
        print(f"\talpha_tau = {alpha}")
        print(f"\tn_t = {n}")
        return alpha / n
    
    def n_zd(self, sigma_S_d, harte_randschicht, **kwargs):
        G_s = self.G_zd_s(**kwargs)
        print(f"\tG_zd_s = {G_s}")
        if harte_randschicht:
            return 1 + m.sqrt(G_s) * 10 ** -0.7
        else:
            return 1 + m.sqrt(G_s) * 10 ** -(0.33 + sigma_S_d / 712)
    def n_b(self, sigma_S_d, harte_randschicht, **kwargs):
        G_s = self.G_b_s(**kwargs)
        print(f"\tG_b_s = {G_s}")
        if harte_randschicht:
            return 1 + m.sqrt(G_s) * 10 ** -0.7
        else:
            return 1 + m.sqrt(G_s) * 10 ** -(0.33 + sigma_S_d / 712)
    def n_t(self, sigma_S_d, harte_randschicht, **kwargs):
        G_s = self.G_t_s(**kwargs)
        print(f"\tG_t_s = {G_s}")
        if harte_randschicht:
            return 1 + m.sqrt(G_s) * 10 ** -0.7
        else:
            return 1 + m.sqrt(G_s) * 10 ** -(0.33 + sigma_S_d / 712)

class AbsatzUndRundnut(BekannteFormzahl):
    umdrehungskerbe = True

    def _alpha(self, A, B, C, z):
        alpha = 1 + 1 / m.sqrt(A * self.r / self.t + 2 * B * self.r / self.d * (1 + 2 * self.r / self.d)**2 + C * (self.r / self.t)**z * self.d / self.D)
        assert self.r / self.t >= 0.03
        assert self.d / self.D <= 0.98
        assert alpha <= 6
        return alpha
    def alpha_sigma_zd(self):
        return self._alpha(self.A_zd(), self.B_zd(), self.C_zd(), self.z_zd())
    def alpha_sigma_b(self):
        return self._alpha(self.A_b(), self.B_b(), self.C_b(), self.z_b())
    def alpha_tau(self):
        return self._alpha(self.A_t(), self.B_t(), self.C_t(), self.z_t())
    
@dataclass
class Rundnut(AbsatzUndRundnut):
    """Tabelle 2"""
    
    r : float
    t : float
        
    def phi(self):
        if self.d / self.D(self.d) > 0.67 and self.r > 0:
            return 1 / (4 * m.sqrt(self.t / self.r) + 2)
        return 0

    def G_zd_s(self, **kwargs):
        print(f"\tϕ = {self.phi(**kwargs)}")
        return 2 * (1 + self.phi(**kwargs)) / self.r
    def G_b_s(self, **kwargs):
        print(f"\tϕ = {self.phi(**kwargs)}")
        return 2 * (1 + self.phi(**kwargs)) / self.r
    def G_t_s(self, **kwargs):
        print(f"\tϕ = {self.phi(**kwargs)}")
        return 1 / self.r
    
    def A_zd(self):
        return 0.22
    def A_b(self):
        return 0.2
    def A_t(self):
        return 0.7

    def B_zd(self):
        return 1.37
    def B_b(self):
        return 2.75
    def B_t(self):
        return 10.3
    
    def C_zd(self):
        return 0
    def C_b(self):
        return 0
    def C_t(self):
        return 0

    def z_zd(self):
        return 0
    def z_b(self):
        return 0
    def z_t(self):
        return 0
    
@dataclass
class Absatz(AbsatzUndRundnut):
    """Tabelle 2"""
    
    r : float
    t : float

    def phi(self, **kwargs):
        if self.d / self.D > 0.67 and self.r > 0:
            return 1 / (4 * m.sqrt(self.t / self.r) + 2)
        return 0

    def G_zd_s(self, **kwargs):
        print(f"\tϕ = {self.phi(**kwargs)}")
        return 2.3 * (1 + self.phi(**kwargs)) / self.r
    def G_b_s(self, **kwargs):
        print(f"\tϕ = {self.phi(**kwargs)}")
        return 2.3 * (1 + self.phi(**kwargs)) / self.r
    def G_t_s(self, **kwargs):
        print(f"\tϕ = {self.phi(**kwargs)}")
        return 1.15 / self.r

    def A_zd(self):
        return 0.62
    def A_b(self):
        return 0.62
    def A_t(self):
        return 3.4

    def B_zd(self):
        return 3.5
    def B_b(self):
        return 5.8
    def B_t(self):
        return 19
    
    def C_zd(self):
        return 0
    def C_b(self):
        return 0.2
    def C_t(self):
        return 1

    def z_zd(self):
        return 0
    def z_b(self):
        return 3
    def z_t(self):
        return 2
    
@dataclass
class Freistrich(Absatz):
    pass

    
@dataclass
class Querbohrung(BekannteFormzahl):
    """Abschnitt 5.2.3 """
    umdrehungskerbe = False
    
    r : float

    def alpha_sigma_zd(self):
        return 3 - (2 * self.r / self.d)
    def alpha_sigma_b(self):
        return 3 + 1.4 * (2 * self.r / self.d) - 2.8 * m.sqrt(2 * self.r / self.d)
    def alpha_tau(self):
        return 2.023 - 1.125 * m.sqrt(2 * self.r / self.d)

    def G_zd_s(self, **kwargs):
        return 2.3 / self.r
    def G_b_s(self, **kwargs):
        return 2.3 / self.r + 2 / self.d
    def G_t_s(self, **kwargs):
        return 1.15 / self.r + 2 / self.d

    def sigma_zd(self, F_zd):
        return F_zd / (m.pi * self.d**2 / 4 - 2 * self.r * self.d)
    def sigma_b(self, M_b):
        return M_b / (m.pi * self.d**3 / 32 - self.r * self.d**2 / 3) * 1000
    def tau_t(self, M_t):
        return M_t / (m.pi * self.d**3 / 16 - self.r * self.d**2 / 3) * 1000


def K_1(werkstoff: Werkstoff, d_eff, zugfestigkeit : bool):
    """
    Glg 10-14
    zugfestigkeit: True für sigma_B, False für sigma_S
    """
    streckgrenze = not zugfestigkeit
    if werkstoff.art == Werkstoff.Art.Nitrierstahl or zugfestigkeit and werkstoff.art == Werkstoff.Art.Baustahl:
        if d_eff <= 100:
            return 1.
        elif d_eff < 300:
            return 1 - 0.23 * m.log(d_eff / 100, 10)
        elif d_eff <= 500:
            return 0.89
    elif streckgrenze and werkstoff.art == Werkstoff.Art.Baustahl:
        if d_eff <= 32:
            return 1.
        elif d_eff < 300:
            d_B = 16
            return 1 - 0.26 * m.log(d_eff / 2 / d_B, 10)
        elif d_eff <= 500:
            return 0.75
    elif werkstoff.art == Werkstoff.Art.CrNiMoEinsatzstahl or zugfestigkeit and Werkstoff.Art.vergüteterStahl:
        d_B = 16
        if d_eff <= 16:
            return 1.
        elif d_eff < 300:
            return 1 - 0.26 * m.log(d_eff / d_B, 10)
        elif d_eff <= 500:
            return 0.67
    elif werkstoff.art == Werkstoff.Art.andererEinsatzstahl:
        if d_eff <= 16:
            return 1.
        elif d_eff < 150:
            d_B = 16
            return 1 - 0.41 * m.log(d_eff / d_B, 10)
        elif d_eff <= 500:
            return 0.6
    elif werkstoff.art == Werkstoff.Art.vergüteterStahl:
        d_B = 16
        if d_eff <= 16:
            return 1.
        elif d_eff < 300:
            return 1 - 0.34 * m.log(d_eff / d_B, 10)
        elif d_eff <= 500:
            return 0.57

    raise NotImplementedError

def K_2_zd(d):
    """Glg 15"""
    return 1
def K_2_b(d):
    """Glg 16"""
    if 7.5 <= d < 150:
        return 1 - 0.2 * m.log(d / 7.5, 10) / m.log(20, 10)
    elif d >= 150:
        return 0.8
    raise NotImplementedError
def K_2_t(d):
    """Glg 16"""
    return K_2_b(d)

def K_3(d, alpha):
    if 7.5 <= d < 150:
        return 1 - 0.2 * m.log(alpha, 10) * m.log(d / 7.5, 10) / m.log(20, 10)
    elif d >= 150:
        return 1 - 0.2 * m.log(alpha)
    raise NotImplementedError