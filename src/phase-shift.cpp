#include "veneziano.h"

int main(int argc, char ** argv)
{
        setprecision(17);
        double s, delta, eta;
        double smax, step;
        int nstep;

        double l = -600.;
        double iso = -600.;
        int plot = -1;

        for (int ii=0; ii<argc; ii++)
        {
                if (strcmp(argv[ii],"-l")==0) l = atof(argv[ii+1]);
                if (strcmp(argv[ii],"-i")==0) iso = atof(argv[ii+1]);
                if (strcmp(argv[ii],"-plot")==0) plot = 1;
        }

        if ((l < 0) || (iso < 0))
        {
                cout << "Enter valid Angular momentum and Isospin. Quitting..." << endl;
                cout << "[USAGE] -l (ang_mom) -i (isospin) " << endl;
                exit(1);
        }

        smax = pow(1.42, 2.);
        step = .0001;

        nstep = (int) ((smax - sthPi) / step);

        ofstream pshift;
        pshift.open("./output/phase_shift.dat");
        pshift << left << setw(30) << "#sqrt(s)" << setw(30) << "Phase_shift" << endl;
        for (int i = 0; i < nstep; i++)
        {
                s = sthPi + .000001 + step * double(i); //offset threshold
                delta = phase_shift(l, iso, s);
                eta = inelasticity(l, iso, s);
                pshift << left << setw(30) << sqrt(s) << setw(30) << delta << setw(30) <<  eta <<  endl;
        }

        //plot if wanted
        if (plot != -1)
        {
                system("gnuplot ./src/gnuplot/phase.gnu");
        }
        return 0;
}
