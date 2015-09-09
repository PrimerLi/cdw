#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define PI acos(-1.0)

double epsilon(double kx, double ky, double kz, double t)
{
    return -2*t*(cos(kx) + cos(ky) + cos(kz));
}

double square(double a)
{
    return a*a;
}

double energy(double epsilon, double Delta)
{
    return sqrt(square(epsilon) + square(Delta));
}

struct Parameters
{
    double lambda;
    double K;
    double t;
    int Nx, Ny, Nz;
    double beta;
    Parameters(double lambda, double K, double t, int Nx, int Ny, int Nz, double beta)
    {
        this->lambda = lambda;
        this->K = K;
        this->t = t;
        this->Nx = Nx;
        this->Ny = Ny;
        this->Nz = Nz;
        this->beta = beta;
    }
};

double f(double Delta, const Parameters &p)
{
    double s = 0;
    double kx, ky, kz;
    int nx, ny, nz;
    for (nx = -p.Nx/4; nx < p.Nx/4; ++nx)
    {
        for (ny = -p.Ny/2; ny < p.Ny/2; ++ny)
        {
            for (nz = -p.Nz/2; nz < p.Nz/2; ++nz)
            {
                kx = 2*nx*PI/p.Nx;
                ky = 2*ny*PI/p.Ny;
                kz = 2*nz*PI/p.Nz;
                s = s + (2*square(p.lambda)/p.K)*(1/energy(epsilon(kx, ky, kz, p.t), Delta))*tanh(p.beta*energy(epsilon(kx, ky, kz, p.t), Delta));
            }
        }
    }
    s = s/(p.Nx/2*p.Ny*p.Nz);
    return s-1;
}

double fprime(double Delta, const Parameters &p)
{
    double delta = 0.0001;
    return (f(Delta + delta, p) - f(Delta, p))/delta;
}

double newton(double (*f)(double, const Parameters &), const Parameters &p)
{
    double x0 = 0.3;
    double x1 = x0 - f(x0, p)/fprime(x0, p);
    int count = 0;
    int iterationMax = 6;
    while(fabs(x0 - x1) > 0.01)
    {
        count++;
        x0 = x1;
        x1 = x0 - f(x0, p)/fprime(x0, p);
        if (count > iterationMax) exit(-1);
    }
    return x1;
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cout << "T = argv[1]. " << endl;
        return -1;
    }
    int Nx, Ny, Nz;
    int nx, ny, nz;
    double t;
    double K;
    double lambda;
    double T;
    double beta;
    int N;

    T = atof(argv[1]);
    beta = 1.0/T;
    t = 1;
    K = 1;
    lambda = 0.65;
    Nx = 100;
    Ny = 100;
    Nz = 100;

    double Delta;
    Parameters p(lambda, K, t, Nx, Ny, Nz, beta);

    /*double delta = 0.05;
    int grids = 40;
    ofstream ofile;
    ofile.open("function.txt");
    for (int i = 0; i < grids; ++i)
    {
        double temp = (i+0.5)*delta;
        ofile << temp << "    " << f(temp, p) << endl;
    }
    ofile.close();
    ofile.open("derived.txt");
    for (int i = 0; i < grids; ++i)
    {
        double temp = (i+0.5)*delta;
        ofile << temp << "    " << fprime(temp, p) << endl;
    }
    ofile.close();*/

    Delta = newton(f, p);
    cout << T << "    " << Delta << endl;

    /*double Delta;
    double delta;
    int count;
    int iterationMax = 5000;
    delta = 0.001;
    count = 0;
    Delta = 0;
    while(true)
    {
        count++;
        double s = 0;
        double kx, ky;
        for (nx = -Nx/4; nx < Nx/4; ++nx)
        {
            for (ny = -Ny/2; ny < Ny/2; ++ny)
            {
                kx = 2*PI*nx/Nx;
                ky = 2*PI*ny/Ny;
                s = s + (2*square(lambda)/K)*(1.0/energy(epsilon(kx, ky, t), Delta))*tanh(beta*energy(epsilon(kx, ky, t), Delta));
            }
        }
        s = 2*s/N;
        cout << s << endl;
        if (count > iterationMax) break;
        if (fabs(s - 1) < delta) 
        {
            cout << T << "    " << Delta << endl;
            break;
        }
        else Delta = Delta + delta;
    }*/
    return 0;
}
