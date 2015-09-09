#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define PI acos(-1.0)

double square(double a)
{
    return a*a;
}

double epsilon(double k, double t)
{
    return -2*t*cos(k);
}

double energy(double epsilon, double Delta)
{
    return sqrt(epsilon*epsilon + Delta*Delta);
}

struct Parameters
{
    double lambda, K;
    double t;
    int Nx;
    double beta;
    Parameters(double lambda, double K, double t, int Nx, double beta)
    {
        this->lambda = lambda;
        this->K = K;
        this->t = t;
        this->Nx = Nx;
        this->beta = beta;
    }
};

double f(double Delta, const Parameters &p)
{
    double s = 0;
    double kx;
    int nx;
    for (nx = -p.Nx/4; nx < p.Nx/4; ++nx)
    {
        kx = 2*nx*PI/p.Nx;
        s = s + (2*square(p.lambda)/p.K)*(1/energy(epsilon(kx, p.t), Delta))*tanh(p.beta*energy(epsilon(kx, p.t), Delta));
    }
    s = s/(p.Nx/2);
    return s-1;
}

double fprime(double Delta, const Parameters &p)
{
    double delta = 0.0001;
    return (f(Delta + delta, p) - f(Delta, p))/delta;
}

double newton(double (*f)(double, const Parameters &),  const Parameters &p)
{
    double x0 = 0.8;
    double x1 = x0 - f(x0, p)/fprime(x0, p);
    int count = 0;
    int iterationMax = 20;
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
    int nx;
    int N;
    double lambda;
    double K;
    double t;

    double Delta;
    double delta;
    double beta;
    double T;

    N = 100;
    lambda = 0.75;
    K = 1;
    t = 1;
    T = atof(argv[1]);
    beta = 1.0/T;

    Parameters p(lambda, K, t, N, beta);

    delta = 0.05;
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
    ofile.close();

    Delta = newton(f, p);
    cout << T << "    " << Delta << endl;
    
    /*Delta = 0;
    int count = 0;
    int iterationMax = 500;
    while(true)
    {
        count++;
        double s = 0;
        for (int i = -N/4; i < N/4; ++i)
        {
            s = s + (2*lambda*lambda/K)*(1/energy(epsilon(2*PI*i/N, t), Delta))*tanh(beta*energy(epsilon(2*PI*i/N, t), Delta));
        }
        s = 2*s/N;
        if (count > iterationMax)
        {
            break;
        }
        if (fabs(s - 1) < delta)
        {
            cout << T << "    " << Delta << endl;
            break;
        }
        else 
        {
            Delta = Delta + delta;
            continue;
        }
    }*/
    return 0;
}
