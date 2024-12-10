import math as m, sys, os
from DIN_743_2 import *
from DIN_743_3 import *



S_min = 1.2
print(f"S_min = {S_min}")
print()

def enable_virtual_terminal_processing():
    import ctypes
    from ctypes import wintypes

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


class HiddenPrints:
    def __init__(self, switch : bool):
        self.switch = switch

    def __enter__(self):
        if self.switch:
            self._original_stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.switch:
            sys.stdout.close()
            sys.stdout = self._original_stdout



class Festigkeit:
    def __init__(self, fall, werkstoff : Werkstoff, kerbe : Kerbe, d_eff,
    F_zdm, F_zda, M_bm, M_ba, M_tm, M_ta,
    Rz, K_A, K_S, K_V, harte_randschicht, output = True, **kwargs):
        """
        Im Druckbereich sind sigma_zdm und sigma_bm negativ
        - d: Bauteildurchmesser im Kerbquerschnitt 
        - d_eff: für die Wärmebehandlung maßgebender Durchmesser

        Werte die von einem if not hasattr umgeben sind können mit **kwargs überschrieben werden
        """
        
        assert(fall in (1, 2))

        self.fall = fall
        self.werkstoff = werkstoff
        self.kerbe = kerbe
        self.d_eff = d_eff
        self.Rz = Rz
        self.K_A = K_A
        self.K_S = K_S
        self.K_V = K_V
        self.harte_randschicht = harte_randschicht


        with HiddenPrints(not output):
            [print(f"{key} = {value}") for key, value in vars(self).items()]
            self.__dict__.update(kwargs)

            def _max(m, a):
                return m + a if a >= 0 else m - a
            self.sigma_zdm = self.kerbe.sigma_zd(F_zdm)
            self.sigma_zda = self.K_A * self.kerbe.sigma_zd(F_zda)
            self.sigma_zdmax = self.K_S * _max(self.sigma_zdm, self.kerbe.sigma_zd(F_zda))
            print(f"σ_zdm = {self.sigma_zdm}")
            print(f"σ_zda = {self.sigma_zda}")
            print(f"σ_zdmax = {self.sigma_zdmax}")

            self.sigma_bm = self.kerbe.sigma_b(M_bm)
            self.sigma_ba = self.K_A * self.kerbe.sigma_b(M_ba)
            self.sigma_bmax = self.K_S * _max(self.sigma_bm, self.kerbe.sigma_b(M_ba))
            print(f"σ_bm = {self.sigma_bm}")
            print(f"σ_ba = {self.sigma_ba}")
            print(f"σ_bmax = {self.sigma_bmax}")

            self.tau_tm = self.kerbe.tau_t(M_tm)
            self.tau_ta = self.K_A * self.kerbe.tau_t(M_ta)
            self.tau_tmax = self.K_S * _max(self.tau_tm, self.kerbe.tau_t(M_ta))
            print(f"τ_tm = {self.tau_tm}")
            print(f"τ_ta = {self.tau_ta}")
            print(f"τ_tmax = {self.tau_tmax}")

            self.sigma_B_d_B = self.werkstoff.sigma_B_d_B
            self.sigma_S_d_B = self.werkstoff.sigma_S_d_B
            print(f"σ_B(d_B) = {self.sigma_B_d_B}")
            print(f"σ_S(d_B) = {self.sigma_S_d_B}")

            self.sigma_zdW_d_B = self.werkstoff.sigma_zdW_d_B
            self.sigma_bW_d_B = self.werkstoff.sigma_bW_d_B
            self.tau_tW_d_B = self.werkstoff.tau_tW_d_B
            print(f"σ_zdW(d_B) = {self.sigma_zdW_d_B}")
            print(f"σ_bW(d_B) = {self.sigma_bW_d_B}")
            print(f"τ_tW(d_B) = {self.tau_tW_d_B}")

            if not hasattr(self, "K_1B_d_eff"):
                self.K_1B_d_eff = K_1(werkstoff=self.werkstoff, d_eff=self.d_eff, zugfestigkeit=True)
            if not hasattr(self, "K_1S_d_eff"):
                self.K_1S_d_eff = K_1(werkstoff=self.werkstoff, d_eff=self.d_eff, zugfestigkeit=False)
            print(f"K_1B(d_eff) = {self.K_1B_d_eff}")
            print(f"K_1S(d_eff) = {self.K_1S_d_eff}")

            if not hasattr(self, "K_2_d_zd"):
                self.K_2zd_d = K_2_zd(d=self.kerbe.d)
            if not hasattr(self, "K_2_d_b"):
                self.K_2b_d = K_2_b(d=self.kerbe.d)
            if not hasattr(self, "K_2_d_t"):
                self.K_2t_d = K_2_t(d=self.kerbe.d)
            print(f"K_2zd(d) = {self.K_2zd_d}")
            print(f"K_2b(d) = {self.K_2b_d}")
            print(f"K_2t(d) = {self.K_2t_d}")

            # 743-2 Abschnitt 7
            self.sigma_B_d_eff = self.K_1B_d_eff * self.sigma_B_d_B
            print(f"σ_B(d_eff) = {self.sigma_B_d_eff}")
            # keine Ahnung
            self.sigma_B_d = self.sigma_B_d_B * self.K_1B_d_eff
            print(f"σ_B(d) = {self.sigma_B_d}")
            # Abschnitt 5.2
            self.sigma_S_d = self.sigma_S_d_B * self.K_1S_d_eff
            print(f"σ_S(d) = {self.sigma_S_d}")

            if not hasattr(self, "K_Fsigma"):
                self.K_Fsigma = self.kerbe.K_Fsigma(Rz=self.Rz, sigma_B_d_eff=self.sigma_B_d_eff)
            if not hasattr(self, "K_Ftau"):
                self.K_Ftau = self.kerbe.K_Ftau(Rz=self.Rz, sigma_B_d_eff=self.sigma_B_d_eff)
            print(f"K_Fσ = {self.K_Fsigma}")
            print(f"K_Fτ = {self.K_Ftau}")

            if not hasattr(self, "beta_sigma_zd"):
                self.beta_sigma_zd = kerbe.beta_sigma_zd( sigma_B_d=self.sigma_B_d, sigma_S_d=self.sigma_S_d, harte_randschicht=self.harte_randschicht)
            if not hasattr(self, "beta_sigma_b"):
                self.beta_sigma_b = kerbe.beta_sigma_b(sigma_B_d=self.sigma_B_d, sigma_S_d=self.sigma_S_d, harte_randschicht=self.harte_randschicht)
            if not hasattr(self, "beta_tau"):
                self.beta_tau = kerbe.beta_tau(sigma_B_d=self.sigma_B_d, sigma_S_d=self.sigma_S_d, harte_randschicht=self.harte_randschicht)
            print(f"β_σzd = {self.beta_sigma_zd}")
            print(f"β_σb = {self.beta_sigma_b}")
            print(f"β_τ = {self.beta_tau}")
  
            # Glg 8, 9
            if not hasattr(self, "K_sigma_zd"):
                self.K_sigma_zd = (self.beta_sigma_zd / self.K_2zd_d + 1 / self.K_Fsigma - 1) / self.K_V
            if not hasattr(self, "K_sigma_b"):
                self.K_sigma_b = (self.beta_sigma_b / self.K_2b_d + 1 / self.K_Fsigma - 1) / self.K_V
            if not hasattr(self, "K_tau"):
                self.K_tau = (self.beta_tau / self.K_2t_d + 1 / self.K_Ftau - 1) / self.K_V
            print(f"K_σzd = {self.K_sigma_zd}")
            print(f"K_σb = {self.K_sigma_b}")
            print(f"K_τ = {self.K_tau}")

            # Glg 5-7
            self.sigma_zdWK = self.sigma_zdW_d_B * self.K_1B_d_eff /self. K_sigma_zd
            self.sigma_bWK = self.sigma_bW_d_B * self.K_1B_d_eff / self.K_sigma_b
            self.tau_tWK = self.tau_tW_d_B * self.K_1B_d_eff / self.K_tau
            print(f"σ_zdWK = {self.sigma_zdWK}")
            print(f"σ_bWK = {self.sigma_bWK}")
            print(f"τ_tWK = {self.tau_tWK}")

            # Glg 20-22
            self.psi_zdsigmaK = self.sigma_zdWK / (2 * self.K_1B_d_eff * self.sigma_B_d_B - self.sigma_zdWK)
            self.psi_bsigmaK = self.sigma_bWK / (2 * self.K_1B_d_eff * self.sigma_B_d_B - self.sigma_bWK)
            self.psi_tauK = self.tau_tWK / (2 * self.K_1B_d_eff * self.sigma_B_d_B - self.tau_tWK)
            print(f"ψ_zdσK = {self.psi_zdsigmaK}")
            print(f"ψ_bσK = {self.psi_bsigmaK}")
            print(f"ψ_τK = {self.psi_tauK}")

            # Glg 23, 24
            self.sigma_mv = m.sqrt((self.sigma_zdm + self.sigma_bm)**2 + 3 * self.tau_tm**2) 
            self.tau_mv = self.sigma_mv / m.sqrt(3)
            print(f"σ_mv = {self.sigma_mv}")
            print(f"τ_mv = {self.tau_mv}")

            # Tabelle 2
            def gamma_F(beta_sigma):
                if not self.kerbe.umdrehungskerbe or self.harte_randschicht:
                    return 1.
                else:
                    if beta_sigma < 1.5:
                        return 1.
                    elif beta_sigma < 2:
                        return 1.05
                    elif beta_sigma < 3:
                        return 1.1
                    else:
                        return 1.15
            self.gamma_Fzd = gamma_F(self.beta_sigma_zd)
            self.gamma_Fb = gamma_F(self.beta_sigma_b)
            self.gamma_Ft = 1.
            print(f"γ_Fzd = {self.gamma_Fzd}")
            print(f"γ_Fb = {self.gamma_Fb}")
            print(f"γ_Ft = {self.gamma_Ft}")

            # K_2F, Tabelle 3, vernachlässigen der Hohlwelle
            if not hasattr(self, "K_2Fzd"):
                self.K_2Fzd = 1.0
            if not hasattr(self, "K_2Fb"):
                self.K_2Fb = 1.0 if self.harte_randschicht else 1.2
            if not hasattr(self, "K_2Ft"):
                self.K_2Ft = 1.0 if self.harte_randschicht else 1.2
            print(f"K_2Fzd = {self.K_2Fzd}")
            print(f"K_2Fb = {self.K_2Fb}")
            print(f"K_2Ft = {self.K_2Ft}")

            # Glg 31, 32
            self.sigma_zdFK = self.K_1S_d_eff * self.K_2Fzd * self.gamma_Fzd * self.sigma_S_d_B
            self.sigma_bFK = self.K_1S_d_eff * self.K_2Fb * self.gamma_Fb * self.sigma_S_d_B
            self.tau_tFK = self.K_1S_d_eff * self.K_2Ft * self.gamma_Ft * self.sigma_S_d_B / m.sqrt(3)
            print(f"σ_zdFK = {self.sigma_zdFK}")
            print(f"σ_bFK = {self.sigma_bFK}")
            print(f"τ_tFK = {self.tau_tFK}")
    
            assert not (self.sigma_zdm + self.sigma_bm < 0)
            assert not (self.sigma_mv < 0)
            if self.fall == 1:
                def ADK(mv, FK, WK, psi):
                    if mv <= (FK - WK) / (1 - psi):
                        return WK - psi * mv
                    else:
                        return FK - mv
                self.sigma_zdADK = ADK(self.sigma_mv, self.sigma_zdFK, self.sigma_zdWK, self.psi_zdsigmaK)
                self.sigma_bADK = ADK(self.sigma_mv, self.sigma_bFK, self.sigma_bWK, self.psi_bsigmaK)
                self.tau_tADK = ADK(self.tau_mv, self.tau_tFK, self.tau_tWK, self.psi_tauK)
                print(f"σ_zdADK = {self.sigma_zdADK}")
                print(f"σ_bADK = {self.sigma_bADK}")
                print(f"τ_tADK = {self.tau_tADK}")
            else:
                def ADK(mv, a, FK, WK, psi):
                    if mv / a <= (FK - WK) / (WK - FK * psi):
                        print(f"{mv / a} <= {(FK - WK) / (WK - FK * psi)}")
                        return WK / (1 + psi * mv / a)
                    else:
                        print(f"{mv / a} > {(FK - WK) / (WK - FK * psi)}")
                        return FK / (1 + mv / a)
            
                if self.sigma_zda != 0:
                    self.sigma_zdADK = ADK(self.sigma_mv, self.sigma_zda, self.sigma_zdFK, self.sigma_zdWK, self.psi_zdsigmaK)
                    print(f"σ_zdADK = {self.sigma_zdADK}")
                if self.sigma_ba != 0:
                    self.sigma_bADK = ADK(self.sigma_mv, self.sigma_ba, self.sigma_bFK, self.sigma_bWK, self.psi_bsigmaK)
                    print(f"σ_bADK = {self.sigma_bADK}")
                if self.tau_ta != 0:
                    self.tau_tADK = ADK(self.tau_mv, self.tau_ta, self.tau_tFK, self.tau_tWK, self.psi_tauK)
                    print(f"τ_tADK = {self.tau_tADK}")

            temp = 0
            if self.sigma_zda != 0:
                temp += self.sigma_zda / self.sigma_zdADK
            if self.sigma_ba != 0:
                temp += self.sigma_ba / self.sigma_bADK
            temp = temp**2
            if self.tau_ta != 0:
                temp += (self.tau_ta / self.tau_tADK)**2
            if temp == 0:
                print("\033[93mSicherheit gegen Dauerbruch nicht berechnet\033[0m")
            else:
                self.S_Dauer = 1 / m.sqrt(temp)
                print(f"S_Dauer = {self.S_Dauer}")
                if not self.S_Dauer >= S_min:
                    print("\033[91mSicherheit gegen Dauerbruch ist nicht erfuellt\033[0m")

            temp = 0
            if self.sigma_zdmax:
                temp += self.sigma_zdmax / self.sigma_zdFK
            if self.sigma_bmax != 0:
                temp += self.sigma_bmax / self.sigma_bFK
            temp = temp**2
            if self.tau_tmax != 0:
                temp += (self.tau_tmax / self.tau_tFK)**2
            if temp == 0:
                print("\033[93mSicherheit gegen bleibende Verformungen nicht berechnet\033[0m")
            else:
                self.S_Verform = 1 / m.sqrt(temp)
                print(f"S_Verform = {self.S_Verform}")
                if not self.S_Verform >= S_min:
                    print("\033[91mSicherheit gegen bleibende Verformungen ist nicht erfuellt\033[0m")

            print()
        return
    
    pass


if __name__ == "__main__":
    enable_virtual_terminal_processing()
    werkstoff = Werkstoff.S500


    print("Absatz 4")
    drehstarr = Festigkeit(fall = 2,
        werkstoff = werkstoff,
        kerbe = Absatz(d = 56, r = 1, t = 2),
        d_eff = 78,
        F_zdm = 0, 
        F_zda = 0, 
        M_bm = 0,
        M_ba = 0, 
        M_tm = 2933.511,
        M_ta = 0, 
        Rz = 16,
        K_A = 1.75,
        K_S = 2.5,
        K_V = 1,
        harte_randschicht = False)
    print()






    #sys.exit()
    print("Passfeder Lamellenkupplung")
    lamellenkupplung = Festigkeit(fall = 2,
        werkstoff = werkstoff,
        kerbe = Passfeder(d = 30, i = 2),
        d_eff = 56,
        F_zdm = 0,
        F_zda = 0,
        M_bm = 0,
        M_ba = 0,
        M_tm = 535.93,
        M_ta = 0,
        Rz = 16,
        K_A = 1.75,
        K_S = 2.5,
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
        M_bm = 0,
        M_ba = 450.256,
        M_tm = 535.93,
        M_ta = 0,
        Rz = 16,
        K_A = 1.75,
        K_S = 2.5,
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
        M_bm = 0,
        M_ba = 313.678, 
        M_tm = 2933.511,
        M_ta = 0, 
        Rz = 16,
        K_A = 1.75,
        K_S = 2.5,
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
        M_bm = 0,
        M_ba = 0, 
        M_tm = 2933.511,
        M_ta = 0, 
        Rz = 16,
        K_A = 1.75,
        K_S = 2.5,
        K_V = 1,
        harte_randschicht = False)
    l = drehstarr.werkstoff.sigma_S_d_B * 0.9 * (10 - 6) * 70 * drehstarr.kerbe.d / 2 * 2 * 0.75 / 1000
    r = 2933.511 * 1.75
    print(f"{l} >= {r}")
    print(l/r)
    print()
    assert l >= r




    



    sys.exit()
    # MEL 1 2024W #1
    Festigkeit(fall = 1,
        werkstoff = Werkstoff._50CrMo4,
        kerbe = Absatz(4, 8),
        d = 84,
        d_eff = 100,
        F_zdm = 443341,
        F_zda = 0,
        M_bm = 23280,
        M_ba = 2910,
        M_tm = 13970,
        M_ta = 0,
        Rz = 10,
        K_A = 1.,
        K_S = 1.,
        K_V = 1,
        harte_randschicht = False)
    # MEL 1 2024W #2
    Festigkeit(fall = 2,
        werkstoff = Werkstoff.C50,
        kerbe = Passfeder(1),
        d = 60,
        d_eff = 80,
        F_zdm = 0,
        F_zda = 0,
        M_bm = 0,
        M_ba = 1402.294,
        M_tm = 1336.902,
        M_ta = 0,
        Rz = 25,
        K_A = 1.,
        K_S = 1.8,
        K_V = 1,
        harte_randschicht = True)
    # MEL 1 2024W #4
    Festigkeit(fall = 2,
        werkstoff = Werkstoff.S235,
        kerbe = Spitzkerbe(10),
        d = 80,
        d_eff = 80,
        F_zdm = 40000,
        F_zda = 0,
        M_bm = 10666.667,
        M_ba = 0,
        M_tm = 0,
        M_ta = 0,
        Rz = 25,
        K_A = 1.,
        K_S = 1.5,
        K_V = 1,
        harte_randschicht = False)

    # MEL 1 2022W
    Festigkeit(fall = 2,
        werkstoff = Werkstoff.E335,
        kerbe = Querbohrung(2),
        d = 50,
        d_eff = 80,
        F_zdm = 0,
        F_zda = 0,
        M_bm = 0,
        M_ba = 1088.802,
        M_tm = 596.831,
        M_ta = 0,
        Rz = 25,
        K_A = 1.,
        K_S = 2.,
        K_V = 1,
        harte_randschicht = False)