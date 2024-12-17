from dataclasses import dataclass
from enum import IntEnum

__all__ = ["Werkstoff", 
              "S185", "S235", "S275", "S355", "S450", "S500", "E295", "E335", "E360",
              "S275N", "S275NL", "S355N", "S355NL", "S420N", "S420NL", "S460N", "S460NL",
              "C10E", "_17Cr3", "_18CrMoS4", "_18CrNiMo7_6", "_16MnCr5", "_20MnCr5",
              "C25", "C30", "C35", "C40", "C45", "C50", "C55", "C60", "_41Cr4", "_34CrMo4", "_42CrMo4", "_50CrMo4", "_36CrNiMo4", "_30CrNiMo8", "_34CrNiMo6",
              "_32CrAlMo7_10", "_34CrAlMo5_10", "_41CrAlMo7_10", "_34CrAlNi7_10", "_8CrMo16_5", "_24CrMo13_6", "_31CrMo12", "_20CrMoV5_7", "_31CrMoV9", "_33CrMoV12_9", "_40CrMoV13_9"
           ]

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
    

# Tabelle A.1
S185 = Werkstoff(Werkstoff.Art.Baustahl, 290, 185, 100, 16)
S235 = Werkstoff(Werkstoff.Art.Baustahl, 360, 235, 100, 16)
S275 = Werkstoff(Werkstoff.Art.Baustahl, 410, 275, 100, 16)
S355 = Werkstoff(Werkstoff.Art.Baustahl, 470, 355, 100, 16)
S450 = Werkstoff(Werkstoff.Art.Baustahl, 550, 450, 100, 16)
S500 = Werkstoff(Werkstoff.Art.Baustahl, 580, 500, 100, 16)
E295 = Werkstoff(Werkstoff.Art.Baustahl, 470, 295, 100, 16)
E335 = Werkstoff(Werkstoff.Art.Baustahl, 570, 335, 100, 16)
E360 = Werkstoff(Werkstoff.Art.Baustahl, 670, 360, 100, 16)

# Tabelle A.2
S275N = Werkstoff(Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, 370, 275, 100, 16)
S275NL = Werkstoff(Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, 370, 275, 100, 16)
S355N = Werkstoff(Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, 470, 355, 100, 16)
S355NL = Werkstoff(Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, 470, 355, 100, 16)
S420N = Werkstoff(Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, 520, 420, 100, 16)
S420NL = Werkstoff(Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, 520, 420, 100, 16)
S460N = Werkstoff(Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, 540, 460, 100, 16)
S460NL = Werkstoff(Werkstoff.Art.schweißgeeigneterFeinkornbaustahl, 540, 460, 100, 16)

# Tabelle A.3
C10E = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 500, 310, 16, 16)
_17Cr3 = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 800, 545, 16, 16)
_18CrMoS4 = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 1100, 775, 16, 16)
_18CrNiMo7_6 = Werkstoff(Werkstoff.Art.CrNiMoEinsatzstahl, 1200, 850, 16, 16)
_16MnCr5 = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 1000, 695, 16, 16)
_20MnCr5 = Werkstoff(Werkstoff.Art.andererEinsatzstahl, 1200, 850, 16, 16)

# Tabelle A.4
C25 = Werkstoff(Werkstoff.Art.vergüteterStahl, 550, 370, 16, 16)
C30 = Werkstoff(Werkstoff.Art.vergüteterStahl, 600, 400, 16, 16)
C35 = Werkstoff(Werkstoff.Art.vergüteterStahl, 630, 430, 16, 16)
C40 = Werkstoff(Werkstoff.Art.vergüteterStahl, 650, 460, 16, 16)
C45 = Werkstoff(Werkstoff.Art.vergüteterStahl, 700, 490, 16, 16)
C50 = Werkstoff(Werkstoff.Art.vergüteterStahl, 750, 520, 16, 16)
C55 = Werkstoff(Werkstoff.Art.vergüteterStahl, 800, 550, 16, 16)
C60 = Werkstoff(Werkstoff.Art.vergüteterStahl, 850, 580, 16, 16)
_41Cr4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1000, 800, 16, 16)
_34CrMo4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1000, 800, 16, 16)
_42CrMo4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1100, 900, 16, 16)
_50CrMo4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1100, 900, 16, 16)
_36CrNiMo4 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1100, 900, 16, 16)
_30CrNiMo8 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1030, 850, 16, 16)
_34CrNiMo6 = Werkstoff(Werkstoff.Art.vergüteterStahl, 1200, 1000, 16, 16)

# Tabelle A.5
_32CrAlMo7_10 = Werkstoff(Werkstoff.Art.Nitrierstahl, 800, 600, 100, 100)
_34CrAlMo5_10 = Werkstoff(Werkstoff.Art.Nitrierstahl, 800, 700, 100, 100)
_41CrAlMo7_10 = Werkstoff(Werkstoff.Art.Nitrierstahl, 850, 650, 100, 100)
_34CrAlNi7_10 = Werkstoff(Werkstoff.Art.Nitrierstahl, 900, 800, 100, 100)
_8CrMo16_5 = Werkstoff(Werkstoff.Art.Nitrierstahl, 900, 720, 100, 100)
_24CrMo13_6 = Werkstoff(Werkstoff.Art.Nitrierstahl, 900, 720, 100, 100)
_31CrMo12 = Werkstoff(Werkstoff.Art.Nitrierstahl, 980, 785, 100, 100)
_20CrMoV5_7 = Werkstoff(Werkstoff.Art.Nitrierstahl, 900, 800, 100, 100)
_31CrMoV9 = Werkstoff(Werkstoff.Art.Nitrierstahl, 1000, 800, 100, 100)
_33CrMoV12_9 = Werkstoff(Werkstoff.Art.Nitrierstahl, 1050, 850, 100, 100)
_40CrMoV13_9 = Werkstoff(Werkstoff.Art.Nitrierstahl, 900, 720, 100, 100)