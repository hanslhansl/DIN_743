import math as m, sys, os
from typing import Literal, Optional
from DIN_743_2 import *
from DIN_743_3 import *


def enable_virtual_terminal_processing():
    import ctypes
    from ctypes import wintypes

    try:
        # Konstanten aus der Windows-API
        STD_OUTPUT_HANDLE = -11
        ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x0004
        INVALID_HANDLE_VALUE = ctypes.c_void_p(-1).value

        # Funktionen aus der Windows-API laden
        kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
        GetStdHandle = kernel32.GetStdHandle
        GetStdHandle.argtypes = [wintypes.DWORD]
        GetStdHandle.restype = wintypes.HANDLE

        GetConsoleMode = kernel32.GetConsoleMode
        GetConsoleMode.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
        GetConsoleMode.restype = wintypes.BOOL

        SetConsoleMode = kernel32.SetConsoleMode
        SetConsoleMode.argtypes = [wintypes.HANDLE, wintypes.DWORD]
        SetConsoleMode.restype = wintypes.BOOL

        # Handle für die Standardausgabe abrufen
        hOut = GetStdHandle(STD_OUTPUT_HANDLE)
        if hOut == INVALID_HANDLE_VALUE:
            raise ctypes.WinError(ctypes.get_last_error())

        # Aktuellen Konsolenmodus abrufen
        dwMode = wintypes.DWORD()
        if not GetConsoleMode(hOut, ctypes.byref(dwMode)):
            raise ctypes.WinError(ctypes.get_last_error())

        # Virtual Terminal Processing aktivieren
        dwMode.value |= ENABLE_VIRTUAL_TERMINAL_PROCESSING
        if not SetConsoleMode(hOut, dwMode):
            raise ctypes.WinError(ctypes.get_last_error())

        print("Virtual Terminal Processing wurde aktiviert.")
    except WindowsError as e:
        print(f"Fehler: {e}")


class Festigkeit:
    def __init__(self,
                fall : Literal[1] | Literal[2],
                werkstoff : Werkstoff,
                kerbe : Kerbe,
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
                S_min : float = 1.2,

                _print = print,
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
            self.K_1B_d_eff = K_1(werkstoff=self.werkstoff, d_eff=self.d_eff, zugfestigkeit=True)
        if self.K_1S_d_eff == None:
            self.K_1S_d_eff = K_1(werkstoff=self.werkstoff, d_eff=self.d_eff, zugfestigkeit=False)
        _print(f"K_1B(d_eff) = {self.K_1B_d_eff}")
        _print(f"K_1S(d_eff) = {self.K_1S_d_eff}")

        if self.K_2zd_d == None:
            self.K_2zd_d = K_2_zd(d=self.kerbe.d)
        if self.K_2b_d == None:
            self.K_2b_d = K_2_b(d=self.kerbe.d)
        if self.K_2t_d == None:
            self.K_2t_d = K_2_t(d=self.kerbe.d)
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
            _print(f"β_σzd = {self.beta_sigmazd}")
        if ba:
            _print(f"β_σb = {self.beta_sigmab}")
        if ta:
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
        
        def gamma_F(alpha):
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
        self.gamma_Fzd = gamma_F(self.kerbe.alpha_sigmazd if hasattr(self.kerbe, "alpha_sigmazd") else self.beta_sigmazd)
        self.gamma_Fb = gamma_F(self.kerbe.alpha_sigmab if hasattr(self.kerbe, "alpha_sigmab") else self.beta_sigmab)
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
    
        assert not (self.sigma_zdm + self.sigma_bm < 0)
        assert not (self.sigma_mv < 0)
        if self.fall == 1:
            def ADK(mv, FK, WK, psi):
                if mv <= (FK - WK) / (1 - psi):
                    _print(f"{mv} <= {(FK - WK) / (1 - psi)}")
                    return WK - psi * mv
                else:
                    _print(f"{mv} > {(FK - WK) / (1 - psi)}")
                    return FK - mv
            self.sigma_zdADK = ADK(self.sigma_mv, self.sigma_zdFK, self.sigma_zdWK, self.psi_zdsigmaK)
            self.sigma_bADK = ADK(self.sigma_mv, self.sigma_bFK, self.sigma_bWK, self.psi_bsigmaK)
            self.tau_tADK = ADK(self.tau_mv, self.tau_tFK, self.tau_tWK, self.psi_tauK)
            if zda:
                _print(f"σ_zdADK = {self.sigma_zdADK}")
            if ba:
                _print(f"σ_bADK = {self.sigma_bADK}")
            if ta:
                _print(f"τ_tADK = {self.tau_tADK}")
        else:
            def ADK(mv, a, FK, WK, psi):
                if mv / a <= (FK - WK) / (WK - FK * psi):
                    _print(f"{mv / a} <= {(FK - WK) / (WK - FK * psi)}")
                    return WK / (1 + psi * mv / a)
                else:
                    _print(f"{mv / a} > {(FK - WK) / (WK - FK * psi)}")
                    return FK / (1 + mv / a)
            
            if self.sigma_zda:
                self.sigma_zdADK = ADK(self.sigma_mv, self.sigma_zda, self.sigma_zdFK, self.sigma_zdWK, self.psi_zdsigmaK)
                if zda:
                    _print(f"σ_zdADK = {self.sigma_zdADK}")
            if self.sigma_ba:
                self.sigma_bADK = ADK(self.sigma_mv, self.sigma_ba, self.sigma_bFK, self.sigma_bWK, self.psi_bsigmaK)
                if ba:
                    _print(f"σ_bADK = {self.sigma_bADK}")
            if self.tau_ta:
                self.tau_tADK = ADK(self.tau_mv, self.tau_ta, self.tau_tFK, self.tau_tWK, self.psi_tauK)
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
        if temp == 0:
            _print("\033[93mSicherheit gegen Dauerbruch nicht berechnet\033[0m")
        else:
            self.S_Dauer = 1 / m.sqrt(temp)
            if self.S_Dauer >= S_min:
                _print(f"S_Dauer = {self.S_Dauer} >= {S_min}")
            else:
                _print(f"S_Dauer = {self.S_Dauer} < {S_min}")
                _print("\033[91mSicherheit gegen Dauerbruch ist nicht erfuellt\033[0m")

        temp = (self.sigma_zdmax / self.sigma_zdFK + self.sigma_bmax / self.sigma_bFK)**2 + (self.tau_tmax / self.tau_tFK)**2
        if temp == 0:
            _print("\033[93mSicherheit gegen bleibende Verformungen nicht berechnet\033[0m")
        else:
            self.S_Verform = 1 / m.sqrt(temp)
            if self.S_Verform >= S_min:
                _print(f"S_Verform = {self.S_Verform} >= {S_min}")
            else:
                _print(f"S_Verform = {self.S_Verform} < {S_min}")
                _print("\033[91mSicherheit gegen bleibende Verformungen ist nicht erfuellt\033[0m")

        _print()
        return
    
    pass



def test():
    _print = lambda *_, **__: None
    print_all = False

    MEL1_2024W_1 = Festigkeit(fall = 1,
        werkstoff = Werkstoff._50CrMo4,
        kerbe = Absatz(84, 4, 8),
        d_eff = 100,
        F_zdm = 443341.5553,
        F_zda = 0,
        F_zdmax = 443341.5553,
        M_bm = 23275.43165,
        M_ba = 2909.428956,
        M_bmax = 26184.86061,
        M_tm = 13965.25899,
        M_ta = 0,
        M_tmax = 13965.25899,
        Rz = 10,
        K_V = 1,
        harte_randschicht = False,
        _print=_print,
        print_all=print_all)
    assert 2.3725 <= MEL1_2024W_1.S_Dauer <= 2.3735
    assert 1.4635 <= MEL1_2024W_1.S_Verform <= 1.4645

    MEL1_2024W_2 = Festigkeit(fall = 2,
        werkstoff = Werkstoff.C50,
        kerbe = Passfeder(60, 1),
        d_eff = 80,
        F_zdm = 0,
        F_zda = 0,
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 1402.294,
        M_bmax = 1402.294 * 1.8,
        M_tm = 1336.902,
        M_ta = 0,
        M_tmax = 1336.902 * 1.8,
        Rz = 25,
        K_V = 1,
        harte_randschicht = True,
        _print=_print,
        print_all=print_all)
    assert 1.4715 <= MEL1_2024W_2.S_Dauer <= 1.4725
    assert 2.5675 <= MEL1_2024W_2.S_Verform <= 2.5685
    
    MEL1_2024W_4 = Festigkeit(fall = 2,
        werkstoff = Werkstoff.S235,
        kerbe = Spitzkerbe(60, 10),
        d_eff = 80,
        F_zdm = 40000,
        F_zda = 0,
        F_zdmax = 40000 * 1.5,
        M_bm = 12666.667,
        M_ba = 0,
        M_bmax = 12666.66667 * 1.5,
        M_tm = 0,
        M_ta = 0,
        M_tmax = 0,
        Rz = 25,
        K_V = 1,
        harte_randschicht = False,
        _print=_print,
        print_all=print_all)
    assert 0.2735 <= MEL1_2024W_4.S_Verform <= 0.2745

    MEL1_2022W = Festigkeit(fall = 2,
        werkstoff = Werkstoff.E335,
        kerbe = Querbohrung(50, 2),
        d_eff = 80,
        F_zdm = 0,
        F_zda = 0,
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 1088.802 * 1,
        M_bmax = 1088.802 * 2,
        M_tm = 596.831,
        M_ta = 0,
        M_tmax = 596.831 * 2,
        Rz = 25,
        K_V = 1,
        harte_randschicht = False,
        _print=_print,
        print_all=print_all)
    assert 1.1045 <= MEL1_2022W.S_Dauer <= 1.1055
    assert 1.6065 <= MEL1_2022W.S_Verform <= 1.6075


    return

if __name__ == "__main__":
    enable_virtual_terminal_processing()


    test()

    sys.exit()

    werkstoff = Werkstoff.S500
    K_A = 1.75
    K_S = 2.5


    print("Passfeder Lamellenkupplung")
    lamellenkupplung = Festigkeit(fall = 2,
        werkstoff = werkstoff,
        kerbe = Passfeder(d = 30, i = 2),
        d_eff = 56,
        F_zdm = 0,
        F_zda = 0,
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 0,
        M_bmax = 0,
        M_tm = 535.93,
        M_ta = 0,
        M_tmax = 535.93 * K_S,
        Rz = 16,
        K_V = 1,
        harte_randschicht = False)
    l = lamellenkupplung.werkstoff.sigma_S_d_B * 0.9 * (7 - 4) * (35) * lamellenkupplung.kerbe.d / 2 * 2 * 0.75 / 1000
    r = 535.93 * 1.75
    print(f"{l} >= {r}")
    print(l/r)
    print()
    assert l >= r

    print("Passfeder Ritzel")
    ritzel = Festigkeit(fall = 2,
        werkstoff = lamellenkupplung.werkstoff,
        kerbe = Passfeder(d = 50, i = 1),
        d_eff = 56,
        F_zdm = 0,
        F_zda = 0,
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 450.256 * K_A,
        M_bmax = 450.256 * K_S,
        M_tm = 535.93,
        M_ta = 0,
        M_tmax = 535.93 * K_S,
        Rz = 16,
        K_V = 1,
        harte_randschicht = False)
    l = ritzel.werkstoff.sigma_S_d_B * 0.9 * (9 - 5.5) * (40) * ritzel.kerbe.d / 2 * 1 * 1 / 1000
    r = 535.93 * 1.75
    print(f"{l} >= {r}")
    print(l/r)
    print()
    assert l >= r

    print("Passfeder Rad")
    rad = Festigkeit(fall = 2,
        werkstoff = werkstoff,
        kerbe = Passfeder(d = 70, i = 2),
        d_eff = 78,
        F_zdm = 0,
        F_zda = 0,
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 313.678 * K_A, 
        M_bmax = 313.678 * K_S,
        M_tm = 2933.511,
        M_ta = 0, 
        M_tmax = 2933.511 * K_S,
        Rz = 16,
        K_V = 1,
        harte_randschicht = False)
    l = rad.werkstoff.sigma_S_d_B * 0.9 * (12 - 7.5) * (48.3) * rad.kerbe.d / 2 * 2 * 0.75 / 1000
    r = 2933.511 * 1.75
    print(f"{l} >= {r}")
    print(l/r)
    print()
    assert l >= r

    print("Passfeder drehstarre Kupplung")
    drehstarr = Festigkeit(fall = 2,
        werkstoff = rad.werkstoff,
        kerbe = Passfeder(d = 55, i = 2),
        d_eff = 78,
        F_zdm = 0, 
        F_zda = 0, 
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 0, 
        M_bmax = 0,
        M_tm = 2933.511,
        M_ta = 0, 
        M_tmax = 2933.511 * K_S,
        Rz = 16,
        K_V = 1,
        harte_randschicht = False)
    l = drehstarr.werkstoff.sigma_S_d_B * 0.9 * (10 - 6) * 70 * drehstarr.kerbe.d / 2 * 2 * 0.75 / 1000
    r = 2933.511 * 1.75
    print(f"{l} >= {r}")
    print(l/r)
    print()
    assert l >= r

    print("Absatz 1")
    absatz1 = Festigkeit(fall = 2,
        werkstoff = werkstoff,
        kerbe = Absatz(d = 60, r = 1, t = 5),
        d_eff = 78,
        F_zdm = 0, 
        F_zda = 0,
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 135.977 * K_A,
        M_bmax = 135.977 * K_S,
        M_tm = 2933.511,
        M_ta = 0,
        M_tmax = 2933.511 * K_S,
        Rz = 16,
        K_V = 1,
        harte_randschicht = False)

    print("Absatz 2")
    absatz2 = Festigkeit(fall = 2,
        werkstoff = werkstoff,
        kerbe = Absatz(d = 70, r = 1, t = 4),
        d_eff = 78,
        F_zdm = 0,
        F_zda = 0,
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 131.175 * K_A,
        M_bmax = 131.175 * K_S,
        M_tm = 2933.511,
        M_ta = 0, 
        M_tmax = 2933.511 * K_S,
        Rz = 16,
        K_V = 1,
        harte_randschicht = False)

    print("Absatz 3")
    absatz3 = Festigkeit(fall = 2,
        werkstoff = werkstoff,
        kerbe = Absatz(d = 60, r = 1, t = 9),
        d_eff = 78,
        F_zdm = 0,
        F_zda = 0,
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 6.003 * K_A,
        M_bmax = 6.003 * K_S,
        M_tm = 2933.511,
        M_ta = 0,
        M_tmax = 2933.511 * K_S,
        Rz = 16,
        K_V = 1,
        harte_randschicht = False)

    print("Freistrich 1")
    freistrich1 = Festigkeit(fall = 2,
        werkstoff = werkstoff,
        kerbe = Freistrich(d = 55.4, r = 1, t = 2.3),
        d_eff = 78,
        F_zdm = 0,
        F_zda = 0,
        F_zdmax = 0,
        M_bm = 0,
        M_ba = 0,
        M_bmax = 0,
        M_tm = 2933.511,
        M_ta = 0, 
        M_tmax = 2933.511 * K_S,
        Rz = 16,
        K_V = 1,
        harte_randschicht = False)

