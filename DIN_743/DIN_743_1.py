import math as m
from typing import Literal, Optional
from . import DIN_743_2, DIN_743_3

__all__ = ["Festigkeit"]

# enable colored console output on windows
import colorama as _colorama
_colorama.just_fix_windows_console()

class Festigkeit:
    """Implementiert die Berechnungen aus DIN 743 Teil 1"""

    def _gamma_F(self, alpha):
            """Tabelle 2"""
            if not self.kerbe.umdrehungskerbe or self.harte_randschicht:
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

    def _ADK(self, mv, a, FK, WK, psi, _print):
        """Glg 10 - 19"""
        if self.fall == 1:
            r = (FK - WK) / (1 - psi)
            if mv <= r:
                _print("Fall 1,", mv, "<=", r)
                return WK - psi * mv        # Glg 10 - 12
            else:
                _print("Fall 1,", mv, ">", r)
                return FK - mv              # Glg 13 & 14
        else:
            l = mv / a
            r = (FK - WK) / (WK - FK * psi)
            if l <= r:
                _print("Fall 2,", l, "<=", r)
                return WK / (1 + psi * l)   # Glg 15 - 17
            else:
                _print("Fall 2,", l, ">", r)
                return FK / (1 + l)         # Glg 18 & 19
        pass

    def _check_value(self, temp, min_or_interv, name, full_name, _print, _assert):
        if temp == 0:
            _print("\033[33m", full_name, " nicht berechnet\033[0m", sep="")
            val = float("inf")
            res = True
        else:
            val = 1 / m.sqrt(temp)
            if isinstance(min_or_interv, tuple):
                res = min_or_interv[0] <= val <= min_or_interv[1]
                if res:
                    _print("\033[32m", name, " = ", val, " in ", min_or_interv, "\033[0m", sep="")
                else:
                    _print("\033[31m", name, " = ", val, " not in ", min_or_interv, ", ", full_name, " ist nicht erfüllt\033[0m", sep="")
            else:
                res = min_or_interv <= val
                if res:
                    _print("\033[32m", name, " = ", val, " >= ", min_or_interv, "\033[0m", sep="")
                else:
                    _print("\033[31m", name, " = ", val, " < ", min_or_interv, ", ", full_name, " ist nicht erfüllt\033[0m", sep="")

        if _assert:
            assert res

        return val
        
    def __init__(self,
                fall : Literal[1, 2],
                werkstoff : DIN_743_3.Werkstoff,
                kerbe : DIN_743_2.Kerbe,
                d_eff : float,
                F_zdm : float, F_zda : float, F_zdmax : float, M_bm : float, M_ba : float, M_bmax : float, M_tm : float, M_ta : float, M_tmax : float,
                Rz : float,
                K_V : float,
                harte_randschicht : bool,

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
                K_2Fzd : Optional[float] = None,
                K_2Fb : Optional[float] = None,
                K_2Ft : Optional[float] = None,
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
        self.K_2Fzd = K_2Fzd
        self.K_2Fb = K_2Fb
        self.K_2Ft = K_2Ft

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
  
        # Glg 8, 9
        if self.K_sigmazd == None:
            self.K_sigmazd = (self.beta_sigmazd / self.K_2zd_d + 1 / self.K_Fsigma - 1) / self.K_V
        if self.K_sigmab == None:
            self.K_sigmab = (self.beta_sigmab / self.K_2b_d + 1 / self.K_Fsigma - 1) / self.K_V
        if self.K_tau == None:
            self.K_tau = (self.beta_tau / self.K_2t_d + 1 / self.K_Ftau - 1) / self.K_V
        if zda:
            _print(f"K_σzd = {self.K_sigmazd}")
        if ba:
            _print(f"K_σb = {self.K_sigmab}")
        if ta:
            _print(f"K_τ = {self.K_tau}")

        # Glg 5-7
        self.sigma_zdWK = self.sigma_zdW_d_B * self.K_1B_d_eff /self. K_sigmazd
        self.sigma_bWK = self.sigma_bW_d_B * self.K_1B_d_eff / self.K_sigmab
        self.tau_tWK = self.tau_tW_d_B * self.K_1B_d_eff / self.K_tau
        if zda:
            _print(f"σ_zdWK = {self.sigma_zdWK}")
        if ba:
            _print(f"σ_bWK = {self.sigma_bWK}")
        if ta:
            _print(f"τ_tWK = {self.tau_tWK}")

        # Glg 20-22
        self.psi_zdsigmaK = self.sigma_zdWK / (2 * self.K_1B_d_eff * self.sigma_B_d_B - self.sigma_zdWK)
        self.psi_bsigmaK = self.sigma_bWK / (2 * self.K_1B_d_eff * self.sigma_B_d_B - self.sigma_bWK)
        self.psi_tauK = self.tau_tWK / (2 * self.K_1B_d_eff * self.sigma_B_d_B - self.tau_tWK)
        if zda:
            _print(f"ψ_zdσK = {self.psi_zdsigmaK}")
        if ba:
            _print(f"ψ_bσK = {self.psi_bsigmaK}")
        if ta:
            _print(f"ψ_τK = {self.psi_tauK}")
            
        # Glg 23, 24
        self.sigma_mv = m.sqrt((self.sigma_zdm + self.sigma_bm)**2 + 3 * self.tau_tm**2) 
        self.tau_mv = self.sigma_mv / m.sqrt(3)
        if zda or ba or ta:
            _print(f"σ_mv = {self.sigma_mv}")
        if ta:
            _print(f"τ_mv = {self.tau_mv}")

        # Tabelle 2
        self.gamma_Fzd = self._gamma_F(self.kerbe.alpha_sigmazd if hasattr(self.kerbe, "alpha_sigmazd") else self.beta_sigmazd)
        self.gamma_Fb = self._gamma_F(self.kerbe.alpha_sigmab if hasattr(self.kerbe, "alpha_sigmab") else self.beta_sigmab)
        self.gamma_Ft = 1.
        if zda or zdmax:
            _print(f"γ_Fzd = {self.gamma_Fzd}")
        if ba or bmax:
            _print(f"γ_Fb = {self.gamma_Fb}")
        if ta or tmax:
            _print(f"γ_Ft = {self.gamma_Ft}")

        # K_2F, Tabelle 3, vernachlässigen der Hohlwelle
        if self.K_2Fzd == None:
            self.K_2Fzd = 1.0
        if self.K_2Fb == None:
            self.K_2Fb = 1.0 if self.harte_randschicht else 1.2
        if self.K_2Ft == None:
            self.K_2Ft = 1.0 if self.harte_randschicht else 1.2
        if zda or zdmax:
            _print(f"K_2Fzd = {self.K_2Fzd}")
        if ba or bmax:
            _print(f"K_2Fb = {self.K_2Fb}")
        if ta or tmax:
            _print(f"K_2Ft = {self.K_2Ft}")

        # Glg 31, 32
        self.sigma_zdFK = self.K_1S_d_eff * self.K_2Fzd * self.gamma_Fzd * self.sigma_S_d_B
        self.sigma_bFK = self.K_1S_d_eff * self.K_2Fb * self.gamma_Fb * self.sigma_S_d_B
        self.tau_tFK = self.K_1S_d_eff * self.K_2Ft * self.gamma_Ft * self.sigma_S_d_B / m.sqrt(3)
        if zda or zdmax:
            _print(f"σ_zdFK = {self.sigma_zdFK}")
        if ba or bmax:
            _print(f"σ_bFK = {self.sigma_bFK}")
        if ta or tmax:
            _print(f"τ_tFK = {self.tau_tFK}")
    
        if self.sigma_zdm + self.sigma_bm < 0:
            raise NotImplementedError
        if self.sigma_mv < 0:
            raise NotImplementedError

        if self.fall == 1:
            self.sigma_zdADK = self._ADK(self.sigma_mv, self.sigma_zda, self.sigma_zdFK, self.sigma_zdWK, self.psi_zdsigmaK, _print)
            self.sigma_bADK = self._ADK(self.sigma_mv, self.sigma_ba, self.sigma_bFK, self.sigma_bWK, self.psi_bsigmaK, _print)
            self.tau_tADK = self._ADK(self.tau_mv, self.tau_ta, self.tau_tFK, self.tau_tWK, self.psi_tauK, _print)
            if zda:
                _print(f"σ_zdADK = {self.sigma_zdADK}")
            if ba:
                _print(f"σ_bADK = {self.sigma_bADK}")
            if ta:
                _print(f"τ_tADK = {self.tau_tADK}")
        else:
            if self.sigma_zda:
                self.sigma_zdADK = self._ADK(self.sigma_mv, self.sigma_zda, self.sigma_zdFK, self.sigma_zdWK, self.psi_zdsigmaK, _print)
                if zda:
                    _print(f"σ_zdADK = {self.sigma_zdADK}")
            if self.sigma_ba:
                self.sigma_bADK = self._ADK(self.sigma_mv, self.sigma_ba, self.sigma_bFK, self.sigma_bWK, self.psi_bsigmaK, _print)
                if ba:
                    _print(f"σ_bADK = {self.sigma_bADK}")
            if self.tau_ta:
                self.tau_tADK = self._ADK(self.tau_mv, self.tau_ta, self.tau_tFK, self.tau_tWK, self.psi_tauK, _print)
                if ta:
                    _print(f"τ_tADK = {self.tau_tADK}")

        temp = 0
        if self.sigma_zda:
            temp += self.sigma_zda / self.sigma_zdADK
        if self.sigma_ba:
            temp += self.sigma_ba / self.sigma_bADK
        temp = temp**2
        if self.tau_ta:
            temp += (self.tau_ta / self.tau_tADK)**2
        self.S_Dauer = self._check_value(temp, S_min, "S_Dauer", "Sicherheit gegen Dauerbruch", _print, _assert)

        temp = (self.sigma_zdmax / self.sigma_zdFK + self.sigma_bmax / self.sigma_bFK)**2 + (self.tau_tmax / self.tau_tFK)**2
        self.S_Verform = self._check_value(temp, S_min, "S_Verform", "Sicherheit gegen bleibende Verformungen", _print, _assert)

        _print()
        return
    
    pass