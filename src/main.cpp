#include "veneziano.h"

int main()
{
        //couplings from file
        double a[maxN+1][maxN+1]; //couplings matrix
        // int num = 0;
        //
        // for (int c = 0; c < maxN + 1; c++)
        // {
        //         num += c;
        // }
        // int n[num], i[num];
        // double ani[num];
        //
        // ifstream couplings;
        // // string filename = "couplings";
        // string line;
        // couplings.open("couplings.txt");
        // if(couplings.fail()) // checks to see if file opended
        // {
        //         cout << "Couldn't open couplings.txt" << endl;
        //         return 1; // no point continuing if the file didn't open...
        // }
        // for (int ii = 0; ii < num; ii++)
        // {
        //         couplings >> n[ii];
        //         couplings >> i[ii];
        //         couplings >> ani[ii];
        //         a[n[ii]][i[ii]] = ani[ii];
        // }
        // couplings.close();
        for (int c = 0; c < 4; c++)
        {
                for(int b = 0; b < 4; b++)
                {
                        a[c][b] = 0.;
                }
        }
        a[1][1] = 11.;
        a[2][1] = 21.;
        a[2][2] = 22.;
        a[3][1] = 31.;
        a[3][2] = 32.;
        a[3][3] = 33.;

        double alph[2] = {.3, .9};

        double s;
        cd amp;

        ofstream pwave;
        pwave.open("./output/rho.dat");
        pwave << left << setw(30) << "#s" << setw(30) << "Re[Amp]" << setw(30) << "Im[Amp]" << endl;
        for (int i = 0; i < 201; i++)
        {
                s = sthPi + .01*double(i);

                amp = VENEZ_isospin_amp(1, a, alph, s, 1. );

                pwave << left << setw(30) << s << setw(30) << real(amp)  << setw(30) << imag(amp) << endl;
        }
        pwave.close();
        system("gnuplot ./src/gnuplot/graph.gnu");
        system("okular ./output/rho.pdf");

        // for (int i = 0; i < 11; i++)
        // {
        // s = 0. + .2*double(i);
        // t = t_man(s, 1.);
        // u = u_man(s, 1.);
        // cout << n_amp(1, alph, a[1], s, t) << endl;
        // }

        return 0;
}
