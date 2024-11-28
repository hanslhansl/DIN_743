from dataclasses import dataclass
from enum import Enum, IntEnum


@dataclass
class Werkstoff:
    class Art(IntEnum):
        Baustahl = 0 # Tabelle A.1
        schweißgeeigneterFeinkornbaustahl = 1   # Tabelle A.2
        CrNiMoEinsatzstahl = 2    # Tabelle A.3
        andererEinsatzstahl = 3    # Tabelle A.3
        vergüteterStahl = 4    # Tabelle A.4
        Nitrierstahl = 5    # Tabelle A.5

    art: Art
    sigma_B_d_B: float
    sigma_S_d_B: float
    d_B_B : float
    d_B_S : float

    @property
    def sigma_bW_d_B(self):
        if self.art in (Werkstoff.Art.Baustahl, Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, Werkstoff.Art.CrNiMoEinsatzstahl,
                        Werkstoff.Art.andererEinsatzstahl, Werkstoff.Art.vergüteterStahl, Werkstoff.Art.Nitrierstahl):
            return 0.5 * self.sigma_B_d_B
        raise NotImplementedError
    @property
    def sigma_zdW_d_B(self):
        if self.art in (Werkstoff.Art.Baustahl, Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, Werkstoff.Art.CrNiMoEinsatzstahl,
                        Werkstoff.Art.andererEinsatzstahl, Werkstoff.Art.vergüteterStahl, Werkstoff.Art.Nitrierstahl):
            return 0.4 * self.sigma_B_d_B
        raise NotImplementedError
    @property
    def tau_tW_d_B(self):
        if self.art in (Werkstoff.Art.Baustahl, Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, Werkstoff.Art.CrNiMoEinsatzstahl,
                        Werkstoff.Art.andererEinsatzstahl, Werkstoff.Art.vergüteterStahl, Werkstoff.Art.Nitrierstahl):
            return 0.3 * self.sigma_B_d_B
        raise NotImplementedError
    
    
Werkstoff.S185 = Werkstoff(Werkstoff.Art.Baustahl, 290, 185, 100, 16)
Werkstoff.S235 = Werkstoff(Werkstoff.Art.Baustahl, 360, 235, 100, 16)
Werkstoff.S275 = Werkstoff(Werkstoff.Art.Baustahl, 410, 275, 100, 16)
Werkstoff.S355 = Werkstoff(Werkstoff.Art.Baustahl, 470, 355, 100, 16)
Werkstoff.S450 = Werkstoff(Werkstoff.Art.Baustahl, 550, 450, 100, 16)
Werkstoff.S500 = Werkstoff(Werkstoff.Art.Baustahl, 580, 500, 100, 16)
Werkstoff.E295 = Werkstoff(Werkstoff.Art.Baustahl, 470, 295, 100, 16)
Werkstoff.E335 = Werkstoff(Werkstoff.Art.Baustahl, 570, 335, 100, 16)
Werkstoff.E360 = Werkstoff(Werkstoff.Art.Baustahl, 670, 360, 100, 16)

Werkstoff.S275N = Werkstoff(Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, 370, 275, 100, 16)

Werkstoff.C10E = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 500, 310, 16, 16)
Werkstoff._17Cr3 = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 800, 545, 16, 16)
Werkstoff._18CrMoS4 = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 1100, 775, 16, 16)
Werkstoff._18CrNiMo7_6 = Werkstoff(Werkstoff.Art.CrNiMoEinsatzstahl, 1200, 850, 16, 16)
Werkstoff._16MnCr5 = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 1000, 695, 16, 16)
Werkstoff._20MnCr5 = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 1200, 850, 16, 16)

Werkstoff.C25 = Werkstoff(Werkstoff.Art.vergüteterStahl, 550, 370, 16, 16)
Werkstoff.C30 = Werkstoff(Werkstoff.Art.vergüteterStahl, 600, 400, 16, 16)
Werkstoff.C35 = Werkstoff(Werkstoff.Art.vergüteterStahl, 630, 430, 16, 16)
Werkstoff.C40 = Werkstoff(Werkstoff.Art.vergüteterStahl, 650, 460, 16, 16)
Werkstoff.C45 = Werkstoff(Werkstoff.Art.vergüteterStahl, 700, 490, 16, 16)
Werkstoff.C50 = Werkstoff(Werkstoff.Art.vergüteterStahl, 750, 520, 16, 16)
Werkstoff.C55 = Werkstoff(Werkstoff.Art.vergüteterStahl, 800, 550, 16, 16)
Werkstoff.C60 = Werkstoff(Werkstoff.Art.vergüteterStahl, 850, 580, 16, 16)
Werkstoff._41Cr4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1000, 800, 16, 16)
Werkstoff._34CrMo4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1000, 800, 16, 16)
Werkstoff._42CrMo4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1100, 900, 16, 16)
Werkstoff._50CrMo4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1100, 900, 16, 16)
Werkstoff._36CrNiMo4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1100, 900, 16, 16)
Werkstoff._30CrNiMo8 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1030, 850, 16, 16)
Werkstoff._34CrNiMo6 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1200, 1000, 16, 16)

Werkstoff._32CrAlMo7_10 = Werkstoff(Werkstoff.Art.Nitrierstahl, 800, 600, 100, 100)
Werkstoff._34CrAlMo5_10 = Werkstoff(Werkstoff.Art.Nitrierstahl, 800, 700, 100, 100)