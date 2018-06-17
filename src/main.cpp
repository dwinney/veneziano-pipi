
// Main program, parses user inputs to fit and/or plot functions.
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "veneziano.h"

int isochoice;

int WaveTranslate(string OPTION)
{
        int wave;
        if (OPTION == "S0") {wave = 0;}
        else if (OPTION == "S2") {wave = 02;}
        else if (OPTION == "P1") {wave = 11;}
        else if (OPTION == "D0") {wave = 20;}
        else if (OPTION == "D2") {wave = 22;}
        else if (OPTION == "F1") {wave = 31;}
        else {cout << "Invalid Partial Wave " << OPTION << ". Quitting..." << endl; cout << endl; exit(1);}
        return wave;
}

int main(int argc, char* argv[])
{
        cout << " " << endl; // Line break to make terminal look pretty :)
        int print = -1;
        int plot = -1;
        string MODEL, AMP, OPTION, OPTION2;

        //Before anything check for plot
        for (int i = 0; i < argc; i++)
        {
                if (strcmp(argv[i], "-plot")==0) plot = 1;
        }

        for (int ii = 0; ii < argc; ii++)
        {
                if (strcmp(argv[ii], "print")==0)
                {
                        if (argc <  4)
                        {
                                cout << "Not enough arguments for plot, check README. Quitting..." << endl;
                                cout << endl;
                                exit(1);
                        }
                        print = 1;
                        MODEL = argv[ii+1]; //GKPY or VENEZ
                        AMP = argv[ii+2]; //isospin, partial, total
                        if (argc > 4)
                        {
                                OPTION = argv[ii+3]; //wave
                        }
                        if (argc > 5)
                        {
                                OPTION2 = argv[ii+4]; //coupling inputfilename
                        }
                        break;
                }
        }

        if (print > 0)
        {
                string OUTPUT;
                OUTPUT = "./output/" + MODEL + "-" + AMP + "-" + OPTION;
                int CASE;
                if (MODEL == "GKPY")
                {

                        if (AMP == "isospin")
                        {
                                CASE = 1;
                        }
                        else if (AMP == "total")
                        {
                                CASE = 2;
                        }
                        else if (AMP == "phaseshift")
                        {
                                CASE = 3;
                        }
                        else if (AMP == "partial")
                        {
                                CASE = 4;
                        }
                        else
                        {
                                cout << " Invalid option for GKPY parameterization. Quitting..." << endl;
                                cout << endl;
                                exit(1);
                        }
                        plotGKPY(CASE, plot, OPTION, OUTPUT);
                }
                else if (MODEL == "VENEZ")
                {
                        if (AMP == "isospin")
                        {
                                CASE = 1;
                        }
                        if (AMP == "total")
                        {
                                CASE = 2;
                                OPTION2 = OPTION;
                        }
                        else if (AMP == "partial")
                        {
                                CASE = 3;
                        }
                        else
                        {
                                cout << " Invalid option for VENEZ parameterization. Quitting..." << endl;
                                cout << endl;
                                exit(1);
                        }
                        plotVENEZ(CASE, plot, OPTION, OPTION2, OUTPUT);
                }
                else
                {
                        cout << "Invalid MODEL chosen. Quiting..." << endl;
                        cout << endl;
                        exit(1);
                }
                //TODO: Use ROOT to plot stuff...
                // if (plot > 0)
                // {
                //         string temp = "gnuplot -p -e \"plot " + OUTPUT +"\" u 1:2 w l";
                //         system(temp.c_str());
                // }
        }
        else
        {
                cout << "Invalid commands entered. See README. Quitting..." << endl;
                cout << endl;
                exit(1);
        }
        return 6;
}
