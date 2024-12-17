


def test_DIN_743():
    import DIN_743

    _print = print #lambda *_, **__: None
    _assert = False
    print_all = True


    _print("MEL1 2024W #1")
    MEL1_2024W_1 = DIN_743.Festigkeit(fall = 1,
        werkstoff = DIN_743._50CrMo4,
        kerbe = DIN_743.Absatz(84, 4, 8),
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
    
    _print("MEL1 2024W #2")
    MEL1_2024W_2 = DIN_743.Festigkeit(fall = 2,
        werkstoff = DIN_743.C50,
        kerbe = DIN_743.Passfeder(60, 1),
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
    
    _print("MEL1 2024W #4")
    MEL1_2024W_4 = DIN_743.Festigkeit(fall = 2,
        werkstoff = DIN_743.S235,
        kerbe = DIN_743.Spitzkerbe(60, 10),
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
    
    _print("MEL1 2022W")
    MEL1_2022W = DIN_743.Festigkeit(fall = 2,
        werkstoff = DIN_743.E335,
        kerbe = DIN_743.Querbohrung(50, 2),
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
        _assert=_assert,
        print_all=print_all)
    assert 1.1045 <= MEL1_2022W.S_Dauer <= 1.1055
    assert 1.6065 <= MEL1_2022W.S_Verform <= 1.6075



test_DIN_743()