import math as m, sys
from DIN_743_2 import *
from DIN_743_3 import *


S_min = 1.2
print(f"S_min = {S_min}")
print()



def Festigkeit(fall, werkstoff : Werkstoff, kerbe : Kerbe, d, d_eff,
    F_zdm, F_zda, M_bm, M_ba, M_tm, M_ta,
    Rz, K_A, K_S, K_V, harte_randschicht, **kwargs):
    """Im Druckbereich sind sigma_zdm und sigma_bm negativ"""

    assert(fall in (1, 2))
    
    print(f"d = {d}")
    print(f"d_eff = {d_eff}")
    
    print(f"Rz = {Rz}")
    print(f"K_A = {K_A}")
    print(f"K_S = {K_S}")
    print(f"K_V = {K_V}")

    def _max(m, a):
        return m + a if a >= 0 else m - a
    sigma_zdm = kerbe.sigma_zd(F_zdm, d)
    sigma_zda = kerbe.sigma_zd(F_zda, d)
    sigma_zdmax = K_S * _max(sigma_zdm, sigma_zda)
    sigma_zda *= K_A
    print(f"sigma_zdm = {sigma_zdm}")
    print(f"sigma_zda = {sigma_zda}")
    print(f"sigma_zdmax = {sigma_zdmax}")

    sigma_bm = kerbe.sigma_b(M_bm, d)
    sigma_ba = kerbe.sigma_b(M_ba, d)
    sigma_bmax = K_S * _max(sigma_bm, sigma_ba)
    sigma_ba *= K_A
    print(f"sigma_bm = {sigma_bm}")
    print(f"sigma_ba = {sigma_ba}")
    print(f"sigma_bmax = {sigma_bmax}")

    tau_tm = kerbe.tau_t(M_tm, d)
    tau_ta = kerbe.tau_t(M_ta, d)
    tau_tmax = K_S * _max(tau_tm, tau_ta)
    tau_ta *= K_A
    print(f"tau_tm = {tau_tm}")
    print(f"tau_ta = {tau_ta}")
    print(f"tau_tmax = {tau_tmax}")

    sigma_B_d_B = werkstoff.sigma_B_d_B
    sigma_S_d_B = werkstoff.sigma_S_d_B
    print(f"sigma_B_d_B = {sigma_B_d_B}")
    print(f"sigma_S_d_B = {sigma_S_d_B}")

    sigma_zdW_d_B = werkstoff.sigma_zdW_d_B
    sigma_bW_d_B = werkstoff.sigma_bW_d_B
    tau_tW_d_B = werkstoff.tau_tW_d_B
    print(f"sigma_zdW_d_B = {sigma_zdW_d_B}")
    print(f"sigma_bW_d_B = {sigma_bW_d_B}")
    print(f"tau_tW_d_B = {tau_tW_d_B}")

    K_1B_d_eff = K_1(werkstoff=werkstoff, d_eff=d_eff, zugfestigkeit=True)
    K_1S_d_eff = K_1(werkstoff=werkstoff, d_eff=d_eff, zugfestigkeit=False)
    print(f"K_1B_d_eff = {K_1B_d_eff}")
    print(f"K_1S_d_eff = {K_1S_d_eff}")

    # 743-2 Abschnitt 7
    sigma_B_d_eff = K_1B_d_eff * sigma_B_d_B
    print(f"sigma_B_d_eff = {sigma_B_d_eff}")
    # keine Ahnung
    sigma_B_d = sigma_B_d_B * K_1B_d_eff
    print(f"sigma_B_d = {sigma_B_d}")
    # Abschnitt 5.2
    sigma_S_d = sigma_S_d_B * K_1S_d_eff
    print(f"sigma_S_d = {sigma_S_d}")

    K_2_d_zd = K_2_zd(d=d)
    K_2_d_b = K_2_b(d=d)
    K_2_d_t = K_2_t(d=d)
    print(f"K_2_d_zd = {K_2_d_zd}")
    print(f"K_2_d_b = {K_2_d_b}")
    print(f"K_2_d_t = {K_2_d_t}")

    K_Fsigma = kerbe.K_Fsigma(Rz=Rz, sigma_B_d_eff=sigma_B_d_eff)
    K_Ftau = kerbe.K_Ftau(Rz=Rz, sigma_B_d_eff=sigma_B_d_eff)
    print(f"K_Fsigma = {K_Fsigma}")
    print(f"K_Ftau = {K_Ftau}")

    beta_sigma_zd = kerbe.beta_sigma_zd(d=d, sigma_B_d=sigma_B_d, sigma_S_d=sigma_S_d, harte_randschicht=harte_randschicht)
    beta_sigma_b = kerbe.beta_sigma_b(d, sigma_B_d=sigma_B_d, sigma_S_d=sigma_S_d, harte_randschicht=harte_randschicht)
    beta_tau = kerbe.beta_tau(d, sigma_B_d=sigma_B_d, sigma_S_d=sigma_S_d, harte_randschicht=harte_randschicht)
    print(f"beta_sigma_zd = {beta_sigma_zd}")
    print(f"beta_sigma_b = {beta_sigma_b}")
    print(f"beta_tau = {beta_tau}")
  
    # Glg 8, 9
    K_sigma_zd = (beta_sigma_zd / K_2_d_zd + 1 / K_Fsigma - 1) / K_V
    K_sigma_b = (beta_sigma_b / K_2_d_b + 1 / K_Fsigma - 1) / K_V
    K_tau = (beta_tau / K_2_d_t + 1 / K_Ftau - 1) / K_V
    print(f"K_sigma_zd = {K_sigma_zd}")
    print(f"K_sigma_b = {K_sigma_b}")
    print(f"K_tau = {K_tau}")

    # Glg 5-7
    sigma_zdWK = sigma_zdW_d_B * K_1B_d_eff / K_sigma_zd
    sigma_bWK = sigma_bW_d_B * K_1B_d_eff / K_sigma_b
    tau_tWK = tau_tW_d_B * K_1B_d_eff / K_tau
    print(f"sigma_zdWK = {sigma_zdWK}")
    print(f"sigma_bWK = {sigma_bWK}")
    print(f"tau_tWK = {tau_tWK}")

    # Glg 20-22
    psi_zdsigmaK = sigma_zdWK / (2 * K_1B_d_eff * sigma_B_d_B - sigma_zdWK)
    psi_bsigmaK = sigma_bWK / (2 * K_1B_d_eff * sigma_B_d_B - sigma_bWK)
    psi_tauK = tau_tWK / (2 * K_1B_d_eff * sigma_B_d_B - tau_tWK)
    print(f"psi_zdsigmaK = {psi_zdsigmaK}")
    print(f"psi_bsigmaK = {psi_bsigmaK}")
    print(f"psi_tauK = {psi_tauK}")

    # Glg 23, 24
    sigma_mv = m.sqrt((sigma_zdm + sigma_bm)**2 + 3 * tau_tm**2) 
    tau_mv = sigma_mv / m.sqrt(3)
    print(f"sigma_mv = {sigma_mv}")
    print(f"tau_mv = {tau_mv}")

    # Tabelle 2
    def gamma_F(beta_sigma):
        if not kerbe.umdrehungskerbe or harte_randschicht:
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
    gamma_Fzd = gamma_F(beta_sigma_zd)
    gamma_Fb = gamma_F(beta_sigma_b)
    gamma_Ft = 1.
    print(f"gamma_Fzd = {gamma_Fzd}")
    print(f"gamma_Fb = {gamma_Fb}")
    print(f"gamma_Ft = {gamma_Ft}")

    # K_2F, Tabelle 3, vernachlässigen der Hohlwelle
    if harte_randschicht:
        K_2Fzd = 1.0
        K_2Fb = 1.0
        K_2Ft = 1.0
    else:
        K_2Fzd = 1.0
        K_2Fb = 1.2
        K_2Ft = 1.2
    print(f"K_2Fzd = {K_2Fzd}")
    print(f"K_2Fb = {K_2Fb}")
    print(f"K_2Ft = {K_2Ft}")

    # Glg 31, 32
    sigma_zdFK = K_1S_d_eff * K_2Fzd * gamma_Fzd * sigma_S_d_B
    sigma_bFK = K_1S_d_eff * K_2Fb * gamma_Fb * sigma_S_d_B
    tau_tFK = K_1S_d_eff * K_2Ft * gamma_Ft * sigma_S_d_B / m.sqrt(3)
    print(f"sigma_zdFK = {sigma_zdFK}")
    print(f"sigma_bFK = {sigma_bFK}")
    print(f"tau_tFK = {tau_tFK}")
    
    assert not (sigma_zdm + sigma_bm < 0)
    assert not (sigma_mv < 0)
    if fall == 1:
        def ADK(mv, FK, WK, psi):
            if mv <= (FK - WK) / (1 - psi):
                return WK - psi * mv
            else:
                return FK - mv
        sigma_zdADK = ADK(sigma_mv, sigma_zdFK, sigma_zdWK, psi_zdsigmaK)
        sigma_bADK = ADK(sigma_mv, sigma_bFK, sigma_bWK, psi_bsigmaK)
        tau_tADK = ADK(tau_mv, tau_tFK, tau_tWK, psi_tauK)
        print(f"sigma_zdADK = {sigma_zdADK}")
        print(f"sigma_bADK = {sigma_bADK}")
        print(f"tau_tADK = {tau_tADK}")
    else:
        def ADK(mv, a, FK, WK, psi):
            if mv / a <= (FK - WK) / (WK - FK * psi):
                return WK / (1 + psi * mv / a)
            else:
                return FK / (1 + mv / a)
            
        if sigma_zda != 0:
            sigma_zdADK = ADK(sigma_mv, sigma_zda, sigma_zdFK, sigma_zdWK, psi_zdsigmaK)
            print(f"sigma_zdADK = {sigma_zdADK}")
        if sigma_ba != 0:
            sigma_bADK = ADK(sigma_mv, sigma_ba, sigma_bFK, sigma_bWK, psi_bsigmaK)
            print(f"sigma_bADK = {sigma_bADK}")
        if tau_ta != 0:
            tau_tADK = ADK(tau_mv, tau_ta, tau_tFK, tau_tWK, psi_tauK)
            print(f"tau_tADK = {tau_tADK}")

    temp = 0
    if sigma_zda != 0:
        temp += sigma_zda / sigma_zdADK
    if sigma_ba != 0:
        temp += sigma_ba / sigma_bADK
    temp = temp**2
    if tau_ta != 0:
        temp += (tau_ta / tau_tADK)**2
    if temp == 0:
        print("\033[93mSicherheit gegen Dauerbruch nicht berechnet\033[0m")
    else:
        S_Dauer = 1 / m.sqrt(temp)
        print(f"S_Dauer = {S_Dauer}")
        if not S_Dauer >= S_min:
            print("\033[91mSicherheit gegen Dauerbruch ist nicht erfuellt\033[0m")

    temp = 0
    if sigma_zdmax:
        temp += sigma_zdmax / sigma_zdFK
    if sigma_bmax != 0:
        temp += sigma_bmax / sigma_bFK
    temp = temp**2
    if tau_tmax != 0:
        temp += (tau_tmax / tau_tFK)**2
    if temp == 0:
        print("\033[93mSicherheit gegen bleibende Verformungen nicht berechnet\033[0m")
    else:
        S_Verform = 1 / m.sqrt(temp)
        print(f"S_Verform = {S_Verform}")
        if not S_Verform >= S_min:
            print("\033[91mSicherheit gegen bleibende Verformungen ist nicht erfuellt\033[0m")

    print()
    return

if __name__ == "__main__":

    # 1
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
    # 2
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
    # 3
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
    # 4
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



    print("Passfeder Ritzel")
    Festigkeit(fall = 2,
        werkstoff = Werkstoff._34CrNiMo6,
        kerbe = Passfeder(1),
        d = 44,
        d_eff = 44,
        F_zdm = 0,
        F_zda = 0,
        M_bm = 0,
        M_ba = 513.681,
        M_tm = 535.93,
        M_ta = 0,
        Rz = 16,
        K_A = 1.75,
        K_S = 2.5,
        K_V = 1,
        harte_randschicht = False)

    """    F_zda = 0,
        F_zdm = 0,
        M_ba = 950.309,
        M_bm = 0,
        M_bmax = 1284.202,
        M_ta = 0,
        M_tm = 991.415,
        M_tmax = 1339.75,
        R_Z = 16,
        K_A = 1.75,
        K_S = 2.5)"""

    print("Passfeder Rad")
    Festigkeit(fall = 2,
        werkstoff = Werkstoff.C60,
        kerbe = Passfeder(1),
        d = 94,
        d_eff = 94,
        F_zdm = 0, 
        F_zda = 0, 
        M_bm = 0,
        M_ba = 513.681, 
        M_tm = 2933.511,
        M_ta = 0, 
        Rz = 16,
        K_A = 1.75,
        K_S = 2.5,
        K_V = 1,
        harte_randschicht = False)

    """F_zda = 0, 
        F_zdm = 0, 
        F_zdmax = 0,
        M_ba = 950.309, 
        M_bm = 0,
        M_bmax = 1284.202,
        M_ta = 0, 
        M_tm = 5390.812,
        M_tmax = 7284.88,"""