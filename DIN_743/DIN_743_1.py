import math as m
from typing import Literal, Optional
from . import DIN_743_2, DIN_743_3



def _safe_division(a : float, b : float):
    if b == 0:
        if a == 0:
            return 0
        return a * float("inf")
    return a / b

def _check_safety(val, min_or_interv, name, full_name, _print):
    is_safe = True
    if val == float("inf"):
        _print("\033[33m", full_name, " nicht berechnet\033[0m", sep="")
    else:
        if isinstance(min_or_interv, tuple):
            res = min_or_interv[0] <= val <= min_or_interv[1]
            if res:
                _print("\033[32m", name, " = ", val, " in ", min_or_interv, "\033[0m", sep="")
            else:
                _print("\033[31m", name, " = ", val, " not in ", min_or_interv, ", ", full_name, " ist nicht erfüllt\033[0m", sep="")
                is_safe = False
        else:
            res = min_or_interv <= val
            if res:
                _print("\033[32m", name, " = ", val, " >= ", min_or_interv, "\033[0m", sep="")
            else:
                _print("\033[31m", name, " = ", val, " < ", min_or_interv, ", ", full_name, " ist nicht erfüllt\033[0m", sep="")
                is_safe = False
    return is_safe


def K_sigmazd(beta_sigmazd : float, K_2zd_d : float, K_Fsigma : float, K_V : float):
    """Glg 8"""
    return (beta_sigmazd / K_2zd_d + 1 / K_Fsigma - 1) / K_V
def K_sigmab(beta_sigmab : float, K_2b_d : float, K_Fsigma : float, K_V : float):
    """Glg 8"""
    return K_sigmazd(beta_sigmab, K_2b_d, K_Fsigma, K_V)
def K_tau(beta_tau : float, K_2t_d : float, K_Ftau : float, K_V : float):
    """Glg 9"""
    return K_sigmazd(beta_tau, K_2t_d, K_Ftau, K_V)

def sigma_zdWK(sigma_zdW_d_B : float, K_1B_d_eff : float, K_sigmazd : float):
    """Glg 5"""
    return sigma_zdW_d_B * K_1B_d_eff / K_sigmazd
def sigma_bWK(sigma_bW_d_B : float, K_1B_d_eff : float, K_sigmab : float):
    """Glg 6"""
    return sigma_zdWK(sigma_bW_d_B, K_1B_d_eff, K_sigmab)
def tau_tWK(tau_tW_d_B : float, K_1B_d_eff : float, K_tau : float):
    """Glg 7"""
    return sigma_zdWK(tau_tW_d_B, K_1B_d_eff, K_tau)

def sigma_zd_bFK(K_1S_d_eff : float, K_2Fzd : float, gamma_Fzd : float, sigma_S_d_B : float):
    """Glg 31"""
    return K_1S_d_eff * K_2Fzd * gamma_Fzd * sigma_S_d_B
def tau_tFK(K_1S_d_eff : float, K_2Ft : float, gamma_Ft : float, sigma_S_d_B : float):
    """Glg 32"""
    return sigma_zd_bFK(K_1S_d_eff, K_2Ft, gamma_Ft, sigma_S_d_B) / m.sqrt(3)

def _ADK(fall : Literal[1, 2], mv : float, a : float, FK : float, WK : float, psi : float, _print):
    if fall == 1:
        l = mv
        r = (FK - WK) / (1 - psi)
        if l <= r:
            res = WK - psi * mv        # Glg 10 - 12
        else:
            res = FK - mv              # Glg 13 & 14
    else:
        l = _safe_division(mv, a)
        r = (FK - WK) / (WK - FK * psi)
        if l <= r:
            res = WK / (1 + psi * l)   # Glg 15 - 17
        else:
            res = FK / (1 + l)         # Glg 18 & 19

    _print("Fall ", fall, ",", l, "<=" if l <= r else ">", r)
    return res
def ADK(fall : Literal[1, 2],
        sigma_zdm : float, sigma_bm : float, tau_tm : float,
        sigma_zda : float, sigma_ba : float, tau_ta : float,
        sigma_zdFK : float, sigma_bFK : float, tau_tFK : float,
        sigma_zdWK : float, sigma_bWK : float, tau_tWK : float,
        psi_zdsigmaK : float, psi_bsigmaK : float, psi_tauK : float,
        _print):
    """Glg 10 - 19, 23 & 24
    Returns σ_mv, τ_mv, σ_zdADK, σ_bADK, τ_tADK"""

    if sigma_zdm + sigma_bm < 0:
        H = m.pow(sigma_bm + sigma_zdm, 3) / abs(sigma_bm + sigma_zdm) + 3 * tau_tm**2
        sigma_mv = H / abs(H) * m.sqrt(abs(H))
    else:
        sigma_mv = m.sqrt((sigma_zdm + sigma_bm)**2 + 3 * tau_tm**2)    # Glg 23

    if sigma_mv < 0:
        tau_mv = 0.
    else:
        tau_mv = sigma_mv / m.sqrt(3)   # Glg 24
    
    return (sigma_mv, tau_mv,
            _ADK(fall, sigma_mv, sigma_zda, sigma_zdFK, sigma_zdWK, psi_zdsigmaK, _print),
            _ADK(fall, sigma_mv, sigma_ba, sigma_bFK, sigma_bWK, psi_bsigmaK, _print),
            _ADK(fall, tau_mv, tau_ta, tau_tFK, tau_tWK, psi_tauK, _print))

def psi_zdsigmaK(sigma_zdWK, K_1B_d_eff, sigma_B_d_B):
    """Glg 20"""
    return sigma_zdWK / (2 * K_1B_d_eff * sigma_B_d_B - sigma_zdWK)
def psi_bsigmaK(sigma_bWK, K_1B_d_eff, sigma_B_d_B):
    """Glg 21"""
    return psi_zdsigmaK(sigma_bWK, K_1B_d_eff, sigma_B_d_B)
def psi_tauK(tau_tWK, K_1B_d_eff, sigma_B_d_B):
    """Glg 22"""
    return psi_zdsigmaK(tau_tWK, K_1B_d_eff, sigma_B_d_B)

def gamma_Fzd(alpha : float, umdrehungskerbe : bool, harte_randschicht : bool):
    """Tabelle 2"""
    if not umdrehungskerbe or harte_randschicht:
        return 1.
    else:
        if alpha < 1.5:
            return 1.
        elif alpha < 2:
            return 1.05
        elif alpha < 3:
            return 1.1
        else:
            return 1.15
def gamma_Fb(alpha : float, umdrehungskerbe : bool, harte_randschicht : bool):
    """Tabelle 2"""
    return gamma_Fzd(alpha, umdrehungskerbe, harte_randschicht)
def gamma_Ft():
    """Tabelle 2"""
    return 1.

def K_2Fzd():
    """Tabelle 3"""
    return 1.
def K_2Fb(harte_randschicht : bool, hohlwelle : bool):
    """Tabelle 3"""
    if harte_randschicht:
        return 1.
    elif hohlwelle:
        return 1.1
    else:
       return 1.2
def K_2Ft(harte_randschicht : bool, hohlwelle : bool):
    """Tabelle 3"""
    if harte_randschicht:
        return 1.
    elif hohlwelle:
        return 1.
    else:
       return 1.2

def S_Dauer(sigma_zda : float, sigma_ba : float, tau_ta : float, sigma_zdADK : float, sigma_bADK : float, tau_tADK : float):
    sqrt = m.sqrt((_safe_division(sigma_zda, sigma_zdADK) + _safe_division(sigma_ba, sigma_bADK))**2 + _safe_division(tau_ta, tau_tADK)**2)
    return _safe_division(1, sqrt)
def S_Verform(sigma_zdmax : float, sigma_bmax : float, tau_tmax : float, sigma_zdFK : float, sigma_bFK : float, tau_tFK : float):
    return S_Dauer(sigma_zdmax, sigma_bmax, tau_tmax, sigma_zdFK, sigma_bFK, tau_tFK)

class Festigkeit:
    """Implementiert die Berechnungen aus DIN 743 Teil 1"""
    def __init__(self,
                fall : Literal[1, 2],
                werkstoff : DIN_743_3.Werkstoff,
                kerbe : DIN_743_2.Kerbe,
                d_eff : float,
                F_zdm : float, F_zda : float, F_zdmax : float, M_bm : float, M_ba : float, M_bmax : float, M_tm : float, M_ta : float, M_tmax : float,
                Rz : float,
                K_V : float,
                harte_randschicht : bool,
                hohlwelle : bool = False,

                K_1B_d_eff : Optional[float] = None,
                K_1S_d_eff : Optional[float] = None,
                K_2zd_d : Optional[float] = None,
                K_2b_d : Optional[float] = None,
                K_2t_d : Optional[float] = None,
                K_Fsigma : Optional[float] = None,
                K_Ftau : Optional[float] = None,
                beta_sigmazd : Optional[float] = None,
                beta_sigmab : Optional[float] = None,
                beta_tau : Optional[float] = None,
                K_sigmazd : Optional[float] = None,
                K_sigmab : Optional[float] = None,
                K_tau : Optional[float] = None,
                S_min : float | tuple[float, float] = 1.2,

                _print = print,
                _assert = False,
                print_all = False):
        """
        Im Druckbereich sind sigma_zdm und sigma_bm negativ
        - d_eff: für die Wärmebehandlung maßgebender Durchmesser
        """
        
        assert(fall in (1, 2))

        self.fall = fall
        self.werkstoff = werkstoff
        self.kerbe = kerbe
        self.d_eff = d_eff
        self.Rz = Rz
        self.K_V = K_V
        self.harte_randschicht = harte_randschicht
        self.hohlwelle = hohlwelle

        self.K_1B_d_eff = K_1B_d_eff
        self.K_1S_d_eff = K_1S_d_eff
        self.K_2zd_d = K_2zd_d
        self.K_2b_d = K_2b_d
        self.K_2t_d = K_2t_d
        self.K_Fsigma = K_Fsigma
        self.K_Ftau = K_Ftau
        self.beta_sigmazd = beta_sigmazd
        self.beta_sigmab = beta_sigmab
        self.beta_tau = beta_tau
        self.K_sigmazd = K_sigmazd
        self.K_sigmab = K_sigmab
        self.K_tau = K_tau

        [_print(f"{key} = {value}") for key, value in vars(self).items() if value != None]

        self.sigma_zdm = self.kerbe.sigma_zd(F_zdm)
        self.sigma_zda = self.kerbe.sigma_zd(F_zda)
        self.sigma_zdmax = self.kerbe.sigma_zd(F_zdmax)
        _print(f"σ_zdm = {self.sigma_zdm}")
        _print(f"σ_zda = {self.sigma_zda}")
        _print(f"σ_zdmax = {self.sigma_zdmax}")

        self.sigma_bm = self.kerbe.sigma_b(M_bm)
        self.sigma_ba = self.kerbe.sigma_b(M_ba)
        self.sigma_bmax = self.kerbe.sigma_b(M_bmax)
        _print(f"σ_bm = {self.sigma_bm}")
        _print(f"σ_ba = {self.sigma_ba}")
        _print(f"σ_bmax = {self.sigma_bmax}")

        self.tau_tm = self.kerbe.tau_t(M_tm)
        self.tau_ta = self.kerbe.tau_t(M_ta)
        self.tau_tmax = self.kerbe.tau_t(M_tmax)
        _print(f"τ_tm = {self.tau_tm}")
        _print(f"τ_ta = {self.tau_ta}")
        _print(f"τ_tmax = {self.tau_tmax}")

        self.sigma_B_d_B = self.werkstoff.sigma_B_d_B
        self.sigma_S_d_B = self.werkstoff.sigma_S_d_B
        _print(f"σ_B(d_B) = {self.sigma_B_d_B}")
        _print(f"σ_S(d_B) = {self.sigma_S_d_B}")

        zda = self.sigma_zda != 0 or print_all
        ba = self.sigma_ba != 0 or print_all
        ta = self.tau_ta != 0 or print_all
        zdmax = self.sigma_zdmax != 0 or print_all
        bmax = self.sigma_bmax != 0 or print_all
        tmax = self.tau_tmax != 0 or print_all

        self.sigma_zdW_d_B = self.werkstoff.sigma_zdW_d_B
        self.sigma_bW_d_B = self.werkstoff.sigma_bW_d_B
        self.tau_tW_d_B = self.werkstoff.tau_tW_d_B
        if zda:
            _print(f"σ_zdW(d_B) = {self.sigma_zdW_d_B}")
        if ba:
            _print(f"σ_bW(d_B) = {self.sigma_bW_d_B}")
        if ta:
            _print(f"τ_tW(d_B) = {self.tau_tW_d_B}")

        if self.K_1B_d_eff == None:
            self.K_1B_d_eff = DIN_743_2.K_1(werkstoff=self.werkstoff, d_eff=self.d_eff, zugfestigkeit=True)
        if self.K_1S_d_eff == None:
            self.K_1S_d_eff = DIN_743_2.K_1(werkstoff=self.werkstoff, d_eff=self.d_eff, zugfestigkeit=False)
        _print(f"K_1B(d_eff) = {self.K_1B_d_eff}")
        _print(f"K_1S(d_eff) = {self.K_1S_d_eff}")

        if self.K_2zd_d == None:
            self.K_2zd_d = DIN_743_2.K_2_zd(d=self.kerbe.d)
        if self.K_2b_d == None:
            self.K_2b_d = DIN_743_2.K_2_b(d=self.kerbe.d)
        if self.K_2t_d == None:
            self.K_2t_d = DIN_743_2.K_2_t(d=self.kerbe.d)
        if zda:
            _print(f"K_2zd(d) = {self.K_2zd_d}")
        if ba:
            _print(f"K_2b(d) = {self.K_2b_d}")
        if ta:
            _print(f"K_2t(d) = {self.K_2t_d}")

        # 743-2 Abschnitt 7
        self.sigma_B_d_eff = self.K_1B_d_eff * self.sigma_B_d_B
        _print(f"σ_B(d_eff) = {self.sigma_B_d_eff}")
        # keine Ahnung
        self.sigma_B_d = self.sigma_B_d_B * self.K_1B_d_eff
        _print(f"σ_B(d) = {self.sigma_B_d}")
        # Abschnitt 5.2
        self.sigma_S_d = self.sigma_S_d_B * self.K_1S_d_eff
        _print(f"σ_S(d) = {self.sigma_S_d}")

        if self.K_Fsigma == None:
            self.K_Fsigma = self.kerbe.K_Fsigma(Rz=self.Rz, sigma_B_d_eff=self.sigma_B_d_eff)
        if self.K_Ftau == None:
            self.K_Ftau = self.kerbe.K_Ftau(Rz=self.Rz, sigma_B_d_eff=self.sigma_B_d_eff)
        if zda or ba:
            _print(f"K_Fσ = {self.K_Fsigma}")
        if ta:
            _print(f"K_Fτ = {self.K_Ftau}")

        if self.beta_sigmazd == None:
            self.beta_sigmazd = kerbe.beta_sigmazd(sigma_B_d=self.sigma_B_d, sigma_S_d=self.sigma_S_d, harte_randschicht=self.harte_randschicht)
        if self.beta_sigmab == None:
            self.beta_sigmab = kerbe.beta_sigmab(sigma_B_d=self.sigma_B_d, sigma_S_d=self.sigma_S_d, harte_randschicht=self.harte_randschicht)
        if self.beta_tau == None:
            self.beta_tau = kerbe.beta_tau(sigma_B_d=self.sigma_B_d, sigma_S_d=self.sigma_S_d, harte_randschicht=self.harte_randschicht)
        if zda:
            _print(kerbe.msg_zd, end="")
            _print(f"β_σzd = {self.beta_sigmazd}")
        if ba:
            _print(kerbe.msg_b, end="")
            _print(f"β_σb = {self.beta_sigmab}")
        if ta:
            _print(kerbe.msg_t, end="")
            _print(f"β_τ = {self.beta_tau}")
  
        if self.K_sigmazd == None:
            self.K_sigmazd = K_sigmazd(self.beta_sigmazd, self.K_2zd_d, self.K_Fsigma, self.K_V)
        if self.K_sigmab == None:
            self.K_sigmab = K_sigmab(self.beta_sigmab, self.K_2b_d, self.K_Fsigma, self.K_V)
        if self.K_tau == None:
            self.K_tau = K_tau(self.beta_tau, self.K_2t_d, self.K_Ftau, self.K_V)
        if zda:
            _print(f"K_σzd = {self.K_sigmazd}")
        if ba:
            _print(f"K_σb = {self.K_sigmab}")
        if ta:
            _print(f"K_τ = {self.K_tau}")

        self.sigma_zdWK = sigma_zdWK(self.sigma_zdW_d_B, self.K_1B_d_eff, self. K_sigmazd)
        self.sigma_bWK = sigma_bWK(self.sigma_bW_d_B, self.K_1B_d_eff, self.K_sigmab)
        self.tau_tWK = tau_tWK(self.tau_tW_d_B, self.K_1B_d_eff, self.K_tau)
        if zda:
            _print(f"σ_zdWK = {self.sigma_zdWK}")
        if ba:
            _print(f"σ_bWK = {self.sigma_bWK}")
        if ta:
            _print(f"τ_tWK = {self.tau_tWK}")

        self.psi_zdsigmaK = psi_zdsigmaK(self.sigma_zdWK, self.K_1B_d_eff, self.sigma_B_d_B)
        self.psi_bsigmaK = psi_bsigmaK(self.sigma_bWK, self.K_1B_d_eff, self.sigma_B_d_B)
        self.psi_tauK = psi_tauK(self.tau_tWK, self.K_1B_d_eff, self.sigma_B_d_B)
        if zda:
            _print(f"ψ_zdσK = {self.psi_zdsigmaK}")
        if ba:
            _print(f"ψ_bσK = {self.psi_bsigmaK}")
        if ta:
            _print(f"ψ_τK = {self.psi_tauK}")
            
        self.gamma_Fzd = gamma_Fzd(self.kerbe.alpha_sigmazd if hasattr(self.kerbe, "alpha_sigmazd") else self.beta_sigmazd, self.kerbe.umdrehungskerbe, self.harte_randschicht)
        self.gamma_Fb = gamma_Fb(self.kerbe.alpha_sigmab if hasattr(self.kerbe, "alpha_sigmab") else self.beta_sigmab, self.kerbe.umdrehungskerbe, self.harte_randschicht)
        self.gamma_Ft = gamma_Ft()
        if zda or zdmax:
            _print("γ_Fzd =", self.gamma_Fzd)
        if ba or bmax:
            _print("γ_Fb =", self.gamma_Fb)
        if ta or tmax:
            _print("γ_Ft =", self.gamma_Ft)

        self.K_2Fzd = K_2Fzd()
        self.K_2Fb = K_2Fb(self.harte_randschicht, self.hohlwelle)
        self.K_2Ft = K_2Ft(self.harte_randschicht, self.hohlwelle)
        if zda or zdmax:
            _print("K_2Fzd =", self.K_2Fzd)
        if ba or bmax:
            _print("K_2Fb =", self.K_2Fb)
        if ta or tmax:
            _print("K_2Ft =", self.K_2Ft)
        
        self.sigma_zdFK = sigma_zd_bFK(self.K_1S_d_eff, self.K_2Fzd, self.gamma_Fzd, self.sigma_S_d_B)
        self.sigma_bFK = sigma_zd_bFK(self.K_1S_d_eff, self.K_2Fb, self.gamma_Fb, self.sigma_S_d_B)
        self.tau_tFK = tau_tFK(self.K_1S_d_eff, self.K_2Ft, self.gamma_Ft, self.sigma_S_d_B)
        if zda or zdmax:
            _print("σ_zdFK =", self.sigma_zdFK)
        if ba or bmax:
            _print("σ_bFK =", self.sigma_bFK)
        if ta or tmax:
            _print("τ_tFK =", self.tau_tFK)
    
        self.sigma_mv, self.tau_mv, self.sigma_zdADK, self.sigma_bADK, self.tau_tADK = ADK(self.fall, self.sigma_zdm, self.sigma_bm, self.tau_tm,
                                                                                           self.sigma_zda, self.sigma_ba, self.tau_ta,
                                                                                           self.sigma_zdFK, self.sigma_bFK, self.tau_tFK,
                                                                                           self.sigma_zdWK, self.sigma_bWK, self.tau_tWK,
                                                                                           self.psi_zdsigmaK, self.psi_bsigmaK, self.psi_tauK,
                                                                                           _print)
        if zda or ba or ta:
            _print(f"σ_mv = {self.sigma_mv}")
        if ta:
            _print(f"τ_mv = {self.tau_mv}")
        if zda:
            _print(f"σ_zdADK = {self.sigma_zdADK}")
        if ba:
            _print(f"σ_bADK = {self.sigma_bADK}")
        if ta:
            _print(f"τ_tADK = {self.tau_tADK}")

        self.S_Dauer = S_Dauer(self.sigma_zda, self.sigma_ba, self.tau_ta, self.sigma_zdADK, self.sigma_bADK, self.tau_tADK)
        self.S_Verform = S_Verform(self.sigma_zdmax, self.sigma_bmax, self.tau_tmax, self.sigma_zdFK, self.sigma_bFK, self.tau_tFK)

        assert _check_safety(self.S_Dauer, S_min, "S_Dauer", "Sicherheit gegen Dauerbruch", _print) or not _assert
        assert _check_safety(self.S_Verform, S_min, "S_Verform", "Sicherheit gegen bleibende Verformungen", _print) or not _assert

        _print()
        return
    
    pass
