#include <fstream>
#include <iostream>
#include <string>

// v0.1.0 Direct translation

/* Original Header
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   ASROT  - ASYMMETRIC ROTOR ENERGY LEVELS, TRANSITIONS, LINE STRENGTHS AND
C            INTENSITIES
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C        - J UP TO 300 AND UP TO 92104 LINES IN PREDICTION
C        - WATSON'S REDUCED HAMILTONIAN (repr Ir, reductions 'A' and 'S')
C          UP TO AND INCLUDING DECADIC TERMS
C        - ALL TYPES OF TRANSITIONS AND BROAD RANGE OF SELECTION RULES ARE
C          SUPPORTED, CLEAR INDICATION OF TRANSITION TYPE IS GIVEN IN THE
C          PRINTOUT
C        - LINE INTENSITIES (MAXIMUM ABSORPTION COEFFICIENTS) ARE CALCULATED
C          FOR A SPECIFIED TEMPERATURE, VALUES ASSUMED FOR OTHER KEY
C          PARAMETERS ARE PRINTED AT THE END OF THE FREQUENCY LISTING
C
C
C       This program, though much developed, is still recognisably descended
C       from ASROT of P.J.Mjoberg et al (ca 1972) - handling of Wang
C       factorisation is common with that used in ASFIT.
C
C
C       ver. 2.06.2015                             ----- Zbigniew KISIEL -----
C                          __________________________________________________
C                         | Institute of Physics, Polish Academy of Sciences |
C                         | Al.Lotnikow 32/46, Warszawa, POLAND              |
C                         |                             kisiel@ifpan.edu.pl  |
C                         |     http://info.ifpan.edu.pl/~kisiel/prospe.htm  |
C_________________________/--------------------------------------------------
C
C  Modification history:
C
c  14.05.95: Warning on reaching the maximum number of lines
c  13.04.97: Overhaul of line packing/unpacking
c  14.04.97: New line sorting code
c  15.04.97: New exit for too many lines
c  29.05.97: More screen information+some mods.
c   4.09.03: Increase to MAXDIM=91204
C  21.03.14: Code tweaks for compilation with contemporary FORTRAN compilers
C   2.06.15: Small mods for Linux compilation
C
C    This program is standard F77 and this version requires about 5 Mb of
C    memory.
C
C_____________________________________________________________________________
C
C    ASROT can be used completely interactively, questions being asked by the
C    program being (hopefully) self explanatory.
C
C    In cases where many constants are to be put in, online typing-in
C    may not be preferred and it is possible to read these from an ASFIT
C    produced or manually edited data file of constants with the following
C    format:
C
C     Line.1:   NCON =    8
C     Line.2:   A    =  48553.1122
C           .   B    =  5010.04024
C           .   C    =  4529.89707
C                DJ  =  .00230008
C                DJK =  -.0691118
C                DK  =  1.74704
C                dJ  =  .000477929
C                dK  =  .029516
C
C      This file is read in under a fixed (A7,F30) FORMAT, the alphanumeric
C      field being a dummy field for constant identification.  NCON can be set
C      as required, but no constants can be skipped i.e. they have to adhere to
C      the following sequence:
C
C      'A' reduction:
C
C      A, B, C
C      DJ, DJK, DK, dJ, dK
C      HJ, HJK, HKJ, HK, hJ, hJK, hK
C      LJ, LJJK, LJK, LKKJ, LK, lJ, lJK, lKJ, lK
C      PJ, PJJK, PJK, PKJ, PKKJ, PK, pJ, pJJK, pJK, pKKJ, pK
C
C
C      'S' reduction:
C
C      A, B, C
C      DJ, DJK, DK, d1, d2
C      HJ, HJK, HKJ, HK, h1, h2, h3
C      LJ, LJJK, LJK, LKKJ, LK, l1, l2, l3, l4
C      PJ, PJJK, PJK, PKJ, PKKJ, PK, p1, p2, p3, p4, p5
C
C
C    - On interactive input the cutoff in constants is specified by means
C      of the order of the Hamiltonian eg. 4 for up to and including quartic,
C      6 - sextic etc.  Any other values are truncated to their nearest lower
C      multiple of 2.
C
C    - Various printout combinations are possible (see appropriate question
C      in the interactive mode), options 1 (frequencies only)+sorting+discard,
C      and option 3 (energy levels only) being recommended for routine use
C
C    - ASROT has built in detection of close K-1, or K+1 doublets on sorting;
C      when these are closer than 0.001 in frequency (FRESOL in SORT0) an
C      entry for a doublet is printed - this has a D after the frequency and
C      only the common K quantum number is printed.
C
C    - ASROT writes identified transitions to a temporary file, each line
C      taking 24 bytes.  For J=300, i.e. 22801 lines this implies
C      a maximum of ca 550Kbytes for this file.  The ASCII output
C      file may be considerably larger, up to over 1.8 Mb
C
C    - if ASROT is run under options: 'frequencies only, lines to be sorted,
C      unsorted lines to be discarded' then the results can be viewed
C      graphically with ASCP_L (or ASCP).
C
C    - On repetitive use of ASROT for similar predictions interactive use
C      becomes tedious and many systems allow a pipeline mechanism such as
C      ASROT<A.IN where A.IN contains all responses required by the program
C_____________________________________________________________________________
C
C     NCON - the number of the highest non-zero constant
C   NCONST - maximum number of usable constants
C
C
C-----------------------------------------------------------------------------
C    Compilation:
C-----------------------------------------------------------------------------
C
C   There are only cosmetic issues with some FORMAT statements on output to
C   screen.  In particular, you might want to uncomment the desired version of
C   2400 format
C
c   Intel/Windows:   ifort -nopdbfile -nodebug -traceback -arch:IA32 -O3 -Qsave
c                    -ccdefault:fortran -fpscomp:filesfromcmd piform.for
C
C   gfortran/Linux:  gfortran -fno-automatic -O2 asrot.for -o asrot
C_____________________________________________________________________________
*/

void pause()
{
    char x;
    std::cin >> x;
}

// Any variable from now on in the range [A-H,O-Z] will be a double unless otherwise stated
// Any variable from now on in the range [I-N] will be a (signed) short unless otherwise stated

// FORTRAN style matrices and C style matrices are transposed relative to each other

constexpr int NCONST = 35;
constexpr int MAXJ = 35;
constexpr int MAXDIM = MAXJ / 2 + 1;
constexpr int LENIST = NCONST + 3;

struct ASRO
{
    double FMIN, FMAX;
    double A[NCONST];
    double EA, EB, EE, EPS;
    double D[3];
    double E[MAXDIM][8];
    double H[MAXDIM][MAXDIM];
    double V[MAXDIM];
    double R;
    double ID[MAXDIM];
    short I, J, L, NT, IP;
};

struct DEF
{
    short NREC;
};

struct BIG
{
    double T;
};

// Increase ISTP size from 6 to 7 for nuill terminated C-style strings
char ISTP[LENIST][7] = 
{
    "A    =", "B    =", "C    =", " DJ  =", " DJK =",
    " DK  =", " dJ  =", " dK  =",
    "HJ   =", "HJK  =", "HKJ  =", "HK   =", "hJ   =", "hJK  =",
    "hK   =",
    " LJ  =", " LJJK=", " LJK =", " LKKJ=", " LK  =", " lJ  =",
    " lJK =", " lKJ =", " lK  =",
    "PJ   =", "PJJK =", "PJK  =", "PKJ  =", "PKKJ =", "PK   =",
    "pJ   =", "pJJK =", "pJK  =", "pKKJ =", "pK   =",
    " MUA =", " MUB =", " MUC ="
};

char BLUFF[6];

// Increase ISTPS size from 6 to 7 for nuill terminated C-style strings
char ISTPS[NCONST][7] = 
{
    "A    =", "B    =", "C    =", " DJ  =", " DJK =",
    " DK  =", " d1  =", " d2  =",
    "HJ   =", "HJK  =", "HKJ  =", "HK   =", "h1   =", "h2   =",
    "h3   =",
    " LJ  =", " LJJK=", " LJK =", " LKKJ=", " LK  =", " l1  =",
    " l2  =", " l3  =", " l4  =",
    "PJ   =", "PJJK =", "PJK  =", "PKJ  =", "PKKJ =", "PK   =",
    "p1   =", "p2   =", "p3   =", "p4   =", "p5   ="
};

char FILNAM[30];
char FILROT[30];

double DIDUTY[18] = 
{
    402657856, 421828864, 402953920,
    422122624, 151267392, 136652417, 270574336,
    285779137, 134493776, 153434257, 287880465,
    268481232, 688395328, 673486593, 537142864,
    556083345, 553924672, 539309697
};

float T[MAXDIM][MAXDIM][8];
float R[MAXDIM][MAXDIM];
float FV = .95f;
float SIGMA = 1;
float GI = 1;
float HALFW = 20;
float ABUND = 1;
float TEMP;
char COMENT[72];
char LINE[79];
char NAMRED[2] = {'A', 'S'};
int NREC, NT;

int main(int argc, char** argv)
{
    // Starting from 181: there seems to be like... 25 lines of FORTRAN code that do nothing

    // I can't make heads or tails of FORTRAN's weird formatting rules  
    // so I'm winging it here based on output from the executable

    std::cout
        << "                                                                           " << '\n'
        << "                                                                           " << '\n'
        << "                                                                           " << '\n'
        << " _________________________________________________________________________ " << '\n'
        << "|  ASROT - PREDICTION OF ASYMMETRIC TOP ROTATIONAL ENERGIES AND SPECTRUM  |" << '\n'
        << "|          (reduction A or S, all terms up to decadic)                    |" << '\n'
        << "|_________________________________________________________________________|" << '\n'
        << " version C++ v0.1.0                                        Zbigniew KISIEL " << '\n'
        << std::endl;

    std::cout << "FILE NAME TO BE USED FOR OUTPUT: ";

    std::cin >> FILNAM;

    std::ofstream output_file(FILNAM, std::ios_base::out);

    std::cout
        << '\n'
        << "PRINTOUT CONTROL:       +-1 = Frequencies only" << '\n'
        << "                        +-2 = Frequencies and energy levels" << '\n'
        << "                        +-3 = Energy levels only" << '\n'
        << '\n'
        << "Positive values=reduction 'A' negative=reduction 'S' ...." << '\n'
        << std::endl;
    
    pause();
}