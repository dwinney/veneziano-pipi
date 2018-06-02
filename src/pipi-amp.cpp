// Contains functions related specifically to pi-pi scattering.
//
// Dependencies: venez-amp.cpp, cgamma.cpp
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@giu.edu
// ---------------------------------------------------------------------------

//TODO: Uncouple from big header.
#include  "veneziano.h"

//TODO: Figure out how to get rid of global variables.

//Isospin projected Veneziano amplitudes for pi-pi scattering.
// iso - isospin (0, 1, or 2)
// coup - couplings matrix a_n,i
// alph - Regge trajectory parameters
// s, t, u - Mandelstam Variables
cd isospin_amp(int iso, double coup[][maxN+1], double alph[], double s, double t, double u)
{
        double couplings[6];
        cd amp = 0.;
        double n_coup[maxN+1];
        for (int n = 1; n < maxN + 1; n++)
        {
                for (int i = 1; i < maxN+1; i++)
                {
                        n_coup[i] = coup[n][i];
                }
                switch (iso) {
                case 0: amp += n_amp(n, alph, n_coup, t, u) - 3. *(n_amp(n, alph, n_coup, s, u) + n_amp(n, alph, n_coup, s, t));
                        break;
                case 1: amp += 2.* (n_amp(n, alph, n_coup, s, u) - n_amp(n, alph, n_coup, s, t));
                        break;
                case 2: amp += -2.*n_amp(n, alph, n_coup, t, u);
                        break;
                default: cout << "Entered isospin not allowed (only 0, 1, or 2)." << endl;
                        exit(1);
                }
        }

        return amp;
}
