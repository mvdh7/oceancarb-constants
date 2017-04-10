import numpy as np

# suggested usage: import oceancarb as pyco2

# === INPUTS ===================================================================
#  temp: temperature / K
#  -  sal: practical salinity
# All Python scripts written by Matthew P. Humphreys
# Last updated: 2017-04-10

# === REFERENCES ===============================================================
#  BE86: Bates & Erickson 1986 JSC
# DSC07: Dickson et al.   2007 PICES Special Publication 3
#  D90a: Dickson          1990 DSRA 10.1016/0198-0149(90)90004-F
# LDK00: Lueker et al.    2000 MC   10.1016/S0304-4203(00)00022-0
# LKB10: Lee et al.       2010 GCA  10.1016/j.gca.2009.12.027
#  lM11: Le Menn          2011 OS   10.5194/os-7-651-2011
# WLD69: Wooster et al.   1969 L&O
#   W71: Warner           1971 DSR  10.1016/0011-7471(71)90030-1
# ==============================================================================

# === FUNCTIONS ================================================================
# === Dissociation/equilibrium constants =======================================
# k2amp: 2-aminopyridine    Seawater pH [BE86]   (temp, sal)
#  k1k2: carbonic acid         Total pH [LDK00]  (temp, sal)
#  kh2o: water                 Total pH []       (temp, sal)
#    kb: boric acid            Total pH []       (temp, sal)
# khso4: bisulfate ion         Total pH []       (temp, sal)
#
# === Concentrations ===========================================================
#    tb: boric acid         micromol/kg [LDK10]  (sal)
#    tf: fluoride                mol/kg [W71]    (sal)
#
# === Miscellaneous ============================================================
#  istr: Ionic strength                 []       (sal)
#
# ==============================================================================


# === Dissociation/equilibrium constants =======================================

def k2amp(temp, sal):
# The 2-aminopyridine stoichiometric dissociation constant
# Empirically determined in synthetic seawater [BE86]

    pKstarSWS_AMP = 2498.31 / temp - 15.3274 + 2.4050 * np.log(temp) \
        + sal * (0.012929 - 2.9417e-5 * temp)

    return pKstarSWS_AMP


def k1k2(temp, sal):
# Carbonic acid stoichiometric equilibrium constants, Total pH [LDK00]
    pKstarT_C1 = 3633.86 / temp - 61.2172 + 9.6777 * np.log(temp) \
        - 0.011555 * sal + 0.0001152 * sal ** 2

    pKstarT_C2 = 471.78 / temp + 25.929 - 3.16967 * np.log(temp) \
        - 0.01781 * sal + 0.0001122 * sal ** 2

    return pKstarT_C1, pKstarT_C2


def kb(temp, sal):
# Boric acid equilibrium constant, Total pH scale
# Equation 23 [D90a]

    lnKB = (-8966.90 - 2890.53 * sal ** 0.5 - 77.942 * sal \
        + 1.728 * sal ** 1.5 - 0.0996 * sal ** 2) / temp \
        + 148.0248 + 137.1942 * sal ** 0.5 + 1.62142 * sal \
        - (24.4344 + 25.085 * sal ** 0.5 + 0.2474 * sal) \
        * np.log(temp) + 0.053105 * sal ** 0.5 * temp

    pKstarT_B = -np.log10(np.exp(lnKB))

    return pKstarT_B


def kh2o(temp, sal):
    # WATER
    lnKW = 148.96502 - 13847.26 / temp - 23.6521 * np.log(temp) \
        + (118.67 / temp - 5.977 + 1.0495 * np.log(temp)) \
        * sal ** 0.5 - 0.01615 * sal

    pKstarT_W = -np.log10(np.exp(lnKW))

    return pKstarT_W


def khso4(temp, sal):
# The ion product of water on the Total pH scale, and
#  the HSO4- stoichiometric dissociation constant on the Free pH scale

    # HSO4
    lnKS =  -4276.1 / temp + 141.328 -  23.093 * np.log(temp) \
        + (-13856.0 / temp + 324.57 - 47.986 * np.log(temp)) \
        * istr(sal) ** 0.5 + (35474.0 / temp - 771.54 + 114.723 \
        * np.log(temp)) * istr(sal)  - (2698.0 / temp) * istr(sal) ** 1.5 \
        + (1776.0 / temp) * istr(sal) ** 2 + np.log(1 - 0.001005 * sal)

    pKstarF_HSO4 = -np.log10(np.exp(lnKS))

    return pKstarF_HSO4



# === Concentrations ===========================================================

def tb(sal, *args):
# Estimate total boron from practical salinity in micromol/kg
# Optional options: LKB10 [default], U74

    # Set up and select appropriate set of constants
    B_sal_dict     = {'LKB10': 0.1336, 'U74': 0.232 / 1.80655}
    B_sal_unc_dict = {'LKB10': 0.0005, 'U74': 0.005 / 1.80655}

    if len(args) == 1:
        ref = args[0]
    else:
        ref = 'LKB10'

    # TB calculation
    B_sal = B_sal_dict[ref] # boron/salinity / mg/kg
    B     = B_sal * sal     # boron / mg/kg
    B_RAM = 10.811e3        # relative atomic mass of boron / mg/mol [DSC07]
    TB    = 1e6 * B/B_RAM   # total boron / micromol/kg

    # Uncertainty in TB, assuming perfect sal
    B_sal_unc = B_sal_unc_dict[ref]
    B_unc     = B_sal_unc * sal
    B_RAM_unc = 0.007e3
    TB_unc    = TB * ((B_unc / B) ** 2 + (B_RAM_unc / B_RAM) ** 2) ** 0.5

    # Alternatively, assuming 0.0034 uncertainty in sal [lM11]
    sal_unc = 0.0034
    B_unc_sal = B * ((B_sal_unc / B_sal) ** 2 + (sal_unc / sal) ** 2) ** 0.5
    TB_unc_sal = TB * ((B_unc_sal / B) ** 2 + (B_RAM_unc / B_RAM) ** 2) ** 0.5

    return TB, TB_unc, TB_unc_sal


def tf(sal):
# Estimate total fluoride from practical salinity [W71] in mol/kg

    F_Cl = 6.75e-5
    F_Cl_unc = 0.03e-5

    F_mass = 18.998 # g/mol [CO2 guide]
    F_mass_unc = 0.001 # g/mol [CO2 guide]

    # Calculate tF
    TF = (F_Cl/F_mass) * sal / 1.80655; # mol/kg, factor 1.80655 from WLD69

    # Uncertainty, assuming perfect salinity
    TF_unc = TF * ((F_Cl_unc / F_Cl) ** 2 + (F_mass_unc / F_mass) ** 2) ** 0.5

    return TF, TF_unc


# === Miscellaneous ============================================================

def istr(sal):
# Ionic strength
    istr = 19.924 * sal / (1000.0 - 1.005 * sal)
    return istr
