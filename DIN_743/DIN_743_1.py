import math as m
from typing import Literal



def _safe_division(a : float, b : float):
    if b == 0:
        if a == 0:
            return 0
        return a * float("inf")
    return a / b


def sigma_B_d_eff(K_1B_d_eff : float, sigma_B_d_B : float):
    """743-2 Abschnitt 7"""
    return K_1B_d_eff * sigma_B_d_B
def sigma_B_d(K_1B_d_eff : float, sigma_B_d_B : float):
    """keine Ahnung"""
    return K_1B_d_eff * sigma_B_d_B
def sigma_S_d(K_1S_d_eff : float, sigma_S_d_B : float):
    """Abschnitt 5.2."""
    return K_1S_d_eff * sigma_S_d_B

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


