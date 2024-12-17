from DIN_743 import *


werkstoff = S500
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