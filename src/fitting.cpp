#include "veneziano.h"
#include <TFitter.h>
#include <TROOT.h>
#include <TMath.h>

int main(int argc, char* argv[])
{
        setprecision(10);

        int np = 100;
        int make = -1;
        int iso = -1;
        for (int ii = 0; ii < argc; ii++)
        {
                if (strcmp(argv[ii], "-i")==0) iso = atoi(argv[ii+1]);
                if (strcmp(argv[ii], "-make")==0) make = 1;
                if (strcmp(argv[ii], "-n")==0) np = atoi(argv[ii+1]);
        }

        if (iso < 0)
        { cout << "Input Isospin to fit [-i isospin]. Quitting..." << endl;}

        //If need to make data points
        if (make >=0)
        {
                double s;
                double smax, step;
                complex<double> amp;

                smax = pow(1.42, 2.);
                step = (smax - sthPi) / double(np); //get stepsize

                cout << "Printing Isospin-" << iso << " amplitude to .dat" << endl;
                cout << "Using " << np << " data points. Step size = " << step << " GeV^2" << endl;
                ofstream data;
                data.open("./output/data.dat");
                for (int i = 0; i < np; i++)
                {
                        s = sthPi + 1e-8 + step*double(i);
                        amp = GKPRY_iso_amp(iso, s, 1.);
                        data << left << setw(30) << s << setw(30) << real(amp) << setw(30) << imag(amp) << endl;
                }
                data.close();
        }
        return 0;
}
