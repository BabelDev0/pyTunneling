# CODE FOR CALCULATING THE TUNNELING PROBABILITY OF AN ELECTRON THROUGH 
# AN ARBITRARY 1D POTENTIAL BARRIER WITH THE TRANSFER MATRIX METHOD.

# declare constants 
HBAR = 6.582119514e-16;  # eV*s
MASS_E = 5.6856300620910916e-30; # eV*s^2/nm^2
Q_E = 1.6021766208e-19; # C
EPS0 = 1.4185972717563562e-39; # C^2/(eV*nm)
PI = 3.14159265358979323846;
CONST_C = HBAR*HBAR/(2*MASS_E); # eV*nm^2

DB = 2; # barrier width
DX = 2; # distance between barriers
MOD = DB+DX; # period of the potential
NB = 1; # number of barriers
LM = (DB*NB)+((NB-1)*DX); # last barrier position

def rectangular_barier(x):
    if x >= LM:
        return 0
    return (x%MOD) <= DB

def main():


if __name__ == '__main__':
    main()