#include "veneziano.h"

int main()
{
        setprecision(17);
        double s, delta, eta;
        double smax, step;
        int nstep;

        smax = pow(1.42, 2.);
        step = .0001;

        nstep = (int) ((smax - sthPi) / step);

        ofstream pshift;
        pshift.open("./output/phase_shift.dat");
        pshift << left << setw(30) << "#sqrt(s)" << setw(30) << "Phase_shift" << endl;
        for (int i = 0; i < nstep - 1; i++)
        {
                s = sthPi + .000001 + step * double(i); //offset threshold
                delta = phase_shift(0, 2, s);
                eta = inelasticity(0, 2, s);
                pshift << left << setw(30) << sqrt(s) << setw(30) << delta << setw(30) <<  eta <<  endl;
        }
        system("gnuplot ./src/gnuplot/phase.gnu");
        return 0;
}
