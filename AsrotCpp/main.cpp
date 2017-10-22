// v0.1.xxx Direct translation - Attempting to keep the code mostly the same at this point for 
//          ease of testing/general comparison. Some control structures may change, but variable 
//          names, etc. should be fairly similar.

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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

void ASROHA(short NCON);
void ASROHS();
void DIPTRA(double DIDUTY[], double RJS, short ISON, short KMMAX, bool& iexit);

// really awful pause function just so the console doesn't die on me while testing
void pause()
{
    char x;
    std::cout << "PAUSED - ENTER ANY CHARACTER TO CONTINUE: ";
    std::cin >> x;
}

constexpr int NCONST = 35;
constexpr int MAXJ = 35;
constexpr int MAXDIM = MAXJ / 2 + 1;
constexpr int LENIST = NCONST + 3;

struct asro_t
{
    double FMIN, FMAX;
    double A[NCONST];
    double EA, EB, EE, EPS;
    double D[3];
    double E[MAXDIM][8];
    // saves on init later. Will probably make this a matrix class at some point.
    double H[MAXDIM][MAXDIM] = {0};
    double V[MAXDIM];
    double R;
    double ID[MAXDIM];
    short I, J, L, NT, IP;
} ASRO;

struct def_t
{
    short NREC;
} DEF;

struct big_t
{
    double T;
} BIG;

// Increase ISTP size from 6 to 7 for null terminated C-style strings
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

// Increase ISTPS size from 6 to 7 for null terminated C-style strings
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

double DIDUTY[] = 
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
// string for getline
std::string COMENT;
char LINE[79];
char NAMRED[2] = {'A', 'S'};
int NREC, NT;

const std::string SEPARATOR = std::string(78, '-');
const std::string PREDICTION_OF = "  PREDICTION OF ASYMMETRIC TOP TRANSITION FREQUENCIES AND/OR ENERGY LEVELS";
std::string watsons_reduced = "  (Watson's reduced Hamiltonian, reduction ?, representation I)";
const std::string PRETTY_STARS = "* * * * *";

constexpr short REDUCTION_REPLACEMENT_INDEX = 43;

constexpr char VERSION[] = "0.1.005";

int main(int argc, char** argv)
{
    using namespace std;

    // Starting from 181: error and other formatting code for about 25 lines

    // I can't make heads or tails of FORTRAN's weird formatting rules  
    // so I'm winging it here based on output from the executable

    string header_version = " C++ v";
    header_version += VERSION;
    header_version += " based on 2.VI.2015 by";

    cout
        << "                                                                           \n"
        << "                                                                           \n"
        << "                                                                           \n"
        << " _________________________________________________________________________ \n"
        << "|  ASROT - PREDICTION OF ASYMMETRIC TOP ROTATIONAL ENERGIES AND SPECTRUM  |\n"
        << "|          (reduction A or S, all terms up to decadic)                    |\n"
        << "|_________________________________________________________________________|\n"
        << left << setw(59) << header_version << right << "Zbigniew KISIEL \n"
        << endl;

    cout << "FILE NAME TO BE USED FOR OUTPUT: ";

    cin >> FILNAM;

    ofstream output_file(FILNAM, ios_base::out);

    // todo [3] create error handling function FORTRAN code seems to handle bad input data well
    //          get the value, flush the input buffer, validate... all that good stuff

    do
    {
        cout
            << '\n'
            << "PRINTOUT CONTROL:       +-1 = Frequencies only\n"
            << "                        +-2 = Frequencies and energy levels\n"
            << "                        +-3 = Energy levels only\n"
            << '\n'
            << "Positive values=reduction 'A' negative=reduction 'S' .... ";

        cin >> ASRO.IP;
    } while (ASRO.IP < -3 || ASRO.IP > 3 || ASRO.IP == 0);
    cout << endl;

    short IRED = 0;
    if (ASRO.IP < 0)
    {
        IRED = 1;
        ASRO.IP = ::abs(ASRO.IP);
        for (int i = 0; i < NCONST; ++i)
            strncpy(ISTP[i], ISTPS[i], 7);
    }

    // line 245

    cout << "   JMIN =  ";
    short JMIN;
    cin >> JMIN;

    cout << "  JRMAX =  ";
    short JRMAX;
    cin >> JRMAX;

    cout << "  JQMAX =  ";
    short JQMAX;
    cin >> JQMAX;

    cout << "K-1.MAX =  ";
    short KMMAX;
    cin >> KMMAX;

    short IS;
    do
    {
        cout << '\n' << "LINES TO BE SORTED (1=YES, 0=NO) ?  ";
        cin >> IS;
    } while(IS != 0 && IS != 1);

    short ISON = 0;
    if (ASRO.IP == 1 && IS == 1)
    {
        do
        {
            cout << '\n' << "ARE UNSORTED FREQUENCIES TO BE DISCARDED ?  ";
            cin >> ISON;
        } while (ISON != 0 && ISON != 1);
    }

    cout << '\n' << "COMMENT :" << endl;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    getline(cin, COMENT);
    if (COMENT.size() >= 72)
        COMENT.erase(72);

    // replace the character with the correct reduction
    watsons_reduced[REDUCTION_REPLACEMENT_INDEX] = NAMRED[IRED];

    output_file
        << ' ' << SEPARATOR << '\n'
        << ' ' << COMENT << '\n'
        << ' ' << SEPARATOR << '\n'
        << PREDICTION_OF << '\n'
        << watsons_reduced << '\n'
        << '\n'
        << setw(47) << PRETTY_STARS << '\n'
        << '\n'
        << '\n'
        << endl;

    // line 278

    cout << "\nARE ROTATIONAL CONSTANTS AVAILABLE ON DISK ?  ";
    cin >> ASRO.I;

    short NCON;
    if (ASRO.I == 1)
    {
        // not sure how to read in constants, save for later.
        cout << '\n' << "Reading constants from files is not yet supported." << endl;

        //cout << '\n' << "NAME OF FILE CONTAINING ROTATIONAL CONSTANTS :  ";
        //cin >> FILROT;
    }

    cout << '\n' << "Order of Hamiltonian (eg. 2, 4, 6 etc.) =  ";
    int IORDER;
    cin >> IORDER;
    // ensure IORDER is within bounds (assuming 0 is not an option). 2 <= IORDER <= 10
    IORDER = max(2, min(IORDER, 10));
    // fudge IORDER to make the constant calculation easier. Division by 2 truncates (i.e. 6 and 7 both yield 3), which means no need to test for odds.
    IORDER = IORDER / 2 + 1;

    // Hamiltonian order - number of constants
    // 2 - 3
    // 4 - 8
    // 6 - 15
    // 8 - 24
    // 10 - 35

    NCON = IORDER * IORDER - 1;

    cout << '\n' << "ROTATIONAL AND CENTRIFUGAL DISTORTION CONSTANTS /MHz :" << '\n' << endl;

    // Going rogue again. Curse you FORMAT

    for (int i = 0; i < NCON; ++i)
    {
        cout << setw(25) << ISTP[i];
        // todo [2] when input as an int or in e notation, divide by 10^12. Why though? Find more edge cases...
        cin >> ASRO.A[i];
    }

    cout << '\n' << "CALCULATION WILL BE CARRIED OUT WITH FOLLOWING CONSTANTS:" << '\n' << endl;

    cout << fixed << setprecision(19);
    for (int i = 0; i < NCON; ++i)
    {
        // todo [3] line up decimal points in output, always show 19 decimals places of precision (like 12345.00...00 <- 19 zeros)
        cout << setw(36) << ISTP[i] << setw(36) << ASRO.A[i] << endl;
    }

    if (ASRO.IP != 3)
    {
        cout << '\n' << "DIPOLE MOMENT COMPONENTS /D :" << endl;
        for (int i = 0; i < 3; ++i)
        {
            cout << '\n' << ISTP[NCONST + i];
            cin >> ASRO.D[i];
        }
        cout << '\n' << "     CUTOFF LINE STRENGTH =  ";
        cin >> ASRO.EPS;
        cout << '\n' << "ROTATIONAL TEMPERATURE /K =  ";
        cin >> TEMP;
        cout << '\n' << " LOW FREQUENCY LIMIT /MHz =  ";
        cin >> ASRO.FMIN;
        cout << '\n' << "HIGH FREQUENCY LIMIT /MHz =  ";
        cin >> ASRO.FMAX;
    }

    output_file
        << "\n\n"
        << fixed << setprecision(19)
        << setw(36) << ISTP[0] << setw(30) << ASRO.A[0] << '\n'
        << setw(36) << ISTP[1] << setw(30) << ASRO.A[1] << '\n'
        << setw(36) << ISTP[2] << setw(30) << ASRO.A[2] << '\n'
        << '\n'
        << setw(36) << ISTP[35] << setw(30) << ASRO.D[0] << '\n'
        << setw(36) << ISTP[36] << setw(30) << ASRO.D[1] << '\n'
        << setw(36) << ISTP[37] << setw(30) << ASRO.D[2] << '\n'
        << endl;

    // line 367

    JMIN = ::abs(JMIN);
    JRMAX = min(::abs(JRMAX), MAXJ);
    JQMAX = min(::abs(JQMAX), MAXJ);
    JMIN = ::abs(JMIN);
    short JMAX = max(JRMAX, JQMAX);
    if (KMMAX == 0)
        KMMAX = JMAX;

    short K = 1;
    if (JMIN > JMAX)
    { /* todo [1] error */}

    K = 2;
    if (ASRO.A[0] < ASRO.A[1] || ASRO.A[1] < ASRO.A[2])
    { /* todo [1] error */ }

    double RJR = 0;
    double RJQ = 0;

    if (ASRO.IP == 3)
    {
        cout << NAMRED[IRED];
    }
    else
    {
        ASRO.EPS = ::abs(ASRO.EPS);
        ASRO.FMIN = ::abs(ASRO.FMIN);
        ASRO.FMAX = ::abs(ASRO.FMAX);

        K = 4;

        if (ASRO.FMAX < ASRO.FMIN)
        { /* todo [1] error */ }

        for (int i = 1; i <= 3; ++i)
        {
            RJR *= 16;
            RJQ *= 4;
            if (ASRO.D[3-i] == 0) continue;
            RJR += 15;
            RJQ += 12288;
        }
    }

    // line 400

    cout << "***** CALCULATION IN PROGRESS, NOW WORKING ON J =  " << endl;

    short JQMIN = JMIN;
    if (JQMAX <= JMIN)
    {
        JQMAX = 0;
        JQMIN = 0;
    }

    output_file
        << fixed << setprecision(1)
        << "J limits for P&R-type: " << setw(3) << JMIN << " - " << setw(3) << JMAX << "     "
        << setprecision(6)
        << "Frequency limits:  " << setw(9) << ASRO.FMIN << " - " << setw(9) << ASRO.FMAX << '\n'
        << "J limits for Q-type:   " << setw(3) << JQMIN << " - " << setw(3) << JQMAX << "     "
        << "Line strength cutoff: " << setw(11) << ASRO.EPS << '\n'
        << "K-1 limits:              0 - " << setw(3) << KMMAX << "     "
        << "Alpha cutoff:            1x10-12" << endl;

    for (short J = JMIN; J <= JMAX; ++J)
    {
        cout << setw(4) << J;
        ASRO.L = 15;
        if (IRED != 1)
            ASROHA(NCON);
        else
            ASROHS();
        if (ASRO.IP >= 1)
        {
        }
        double RJS = 1e-6;

        if (J != JMIN && J <= JRMAX)
            RJS += RJR;

        if (J <= JQMAX)
            RJS += RJQ;

        bool iexit = false;
        DIPTRA(DIDUTY, RJS, ISON, KMMAX, iexit);
    }

    output_file.close();

    pause();
    return EXIT_SUCCESS;
}

void ASROHA(short NCON)
{

}

void ASROHS()
{

}

void DIPTRA(double DIDUTY[], double RJS, short ISON, short KMMAX, bool& iexit)
{

}