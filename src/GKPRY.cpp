#include "veneziano.h"

int main(int argc, char ** argv)
{
        setprecision(17);
        double s, delta, eta, sigma;
        double smax, step;
        int nstep;
        complex<double> amp;

        double l = -600.;
        double iso = -600.;
        int plot = -1;
        int open = -1;

        for (int ii=0; ii<argc; ii++)
        {
                if (strcmp(argv[ii],"-l")==0) l = atof(argv[ii+1]);
                if (strcmp(argv[ii],"-i")==0) iso = atof(argv[ii+1]);
                if (strcmp(argv[ii],"-plot")==0) plot = 1;
                if (strcmp(argv[ii],"-open")==0) open = 1;
        }

        smax = pow(1.42, 2.);
        step = .0001;

        nstep = (int) ((smax - sthPi) / step);

        ofstream gkpry;
        gkpry.open("./output/GKPRY.dat");
        gkpry << left << setw(30) << "#sqrt(s)" << setw(30) << "delta"
              << setw(30) << "eta" << setw(30) << "A" << setw(30)
              << "sigma" << endl;
        for (int i = 0; i < nstep; i++)
        {
                s = sthPi + 1e-5 + step * double(i); //offset threshold

                sigma = GKPRY_cross_section(s);
                gkpry << left << setw(30) << sqrt(s);
                cout << sigma << endl;

                if (l >= 0)
                {
                        delta = phase_shift(l, iso, s);
                        eta = inelasticity(l, iso, s);
                        amp = GKPRY_partial_wave(l, iso, s);
                        gkpry << setw(30) << delta << setw(30) <<  eta
                              << setw(30) << real(amp) << setw(30) << imag(amp);
                }
                gkpry << setw(30) << sigma << endl;
        }

        //plot if wanted
        if (plot != -1)
        {
                if(iso >= 0) system("gnuplot ./src/gnuplot/phase.gnu");
                system("gnuplot ./src/gnuplot/crosssection.gnu");
        }
        if (open != -1)
        {
                if(iso >= 0) system("gnuplot ./src/gnuplot/phase.gnu");
                system("okular ./output/crosssection.pdf");
        }
        return 0;
}
