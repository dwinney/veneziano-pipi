//
// Definitions of Definite Isospin amplitudes for pipi scattering
// using Veneziano model.
// ---------------------------------------------------------------------------

#include "veneziano.h"

//TODO Check indexing on couplings to make sure they line up

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
                }
        }
        return amp;
}
