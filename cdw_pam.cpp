#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#define PI acos(-1.0)
using namespace std;

struct Parameters
{
    double g, K, V, t;
    Parameters(double g, double K, double V, double t)
    {
        this->g = g;
        this->K = K;
        this->V = V;
        this->t = t;
    }
};

double square(double a)
{
    return a*a;
}

double epsilon(double kx, double ky, double t)
{
    return -2*t*(cos(kx) + cos(ky));
}

double energy(double kx, double ky, double t, double Delta)
{
    return sqrt(square(epsilon(kx,ky,t)) + square(Delta));
}

double summand(double Ek, double V, double beta)
{
    double one, two;
    one = 1.0/Ek*(1.0 - Ek/sqrt(square(Ek) + 4*square(V)))*tanh(0.25*beta*(Ek - sqrt(square(Ek) + 4*square(V))));
    two = 1.0/Ek*(1.0 + Ek/sqrt(square(Ek) + 4*square(V)))*tanh(0.25*beta*(Ek + sqrt(square(Ek) + 4*square(V))));
    return one + two;
}

double f(const Parameters &p, double beta, double Delta)
{
    double g, K, V, t;
    g = p.g;
    K = p.K;
    V = p.V;
    t = p.t;

    double s = 0;
    int Nx = 100;
    int Ny = 100;
    
    double Ek;
    double kx, ky;
    double increment;
    for (int nx = -Nx/2; nx <= Nx/2; ++nx)
    {
        for (int ny = -Ny/2; ny <= Ny/2; ++ny)
        {
            kx = 2*nx*PI/Nx;
            ky = 2*ny*PI/Ny;
            Ek = energy(kx, ky, t, Delta);
            s = s + summand(Ek, V, beta);
        }
    }
    s = s*2*square(g)/(K*Nx*Ny);
    s = s - 1.0;
    return s;
}

double fprime(const Parameters &p, double beta, double Delta)
{
    double dDelta = 0.00001;
    return (f(p,beta,Delta+dDelta) - f(p,beta,Delta))/dDelta;
}

double newton(double (*f)(const Parameters &p, double beta, double Delta), const Parameters &p, double beta)
{
    double x0 = 3.0;
    double x1 = x0 - f(p, beta, x0)/fprime(p,beta,x0);
    int count = 0;
    int iterationMax = 8;
    while(fabs(x0 - x1) > 0.001)
    {
        count = count + 1;
        x0 = x1;
        x1 = x0 - f(p,beta,x0)/fprime(p,beta,x0);
        if (count > iterationMax) 
        {
            cout << "No solution found. " << endl;
            break;
            //exit(-1);
        }
    }
    return x1;
}

int main(int argc, char **argv)
{
    double g, K, V, t;
    double beta;
    double Delta;
    double increment;

    g = 1.0;
    K = 1.0;
    V = 0;
    t = 1.0;
    Parameters p(g,K,V,t);

    if (false)
    {
        increment = 0.1;
        beta = 1000.0;
        ofstream ofile;
        ofile.open("f.txt");
        for (int i = 0; i < 100; ++i)
        {
            Delta = increment*(i+0.5);
            ofile << Delta << "    " << f(p, beta, Delta) << endl;
        }
        ofile.close();
        return 0;
    }

    if (false)
    {
        double V0 = 0;
        double V1 = 1.8269075;
        int NV = 100;
        double dV = (V1 - V0)/NV;
        beta = 10000;
        for (int i = 0; i <= NV; ++i)
        {
            V = V0 + i*dV;
            Parameters parameters(g, K, V, t);
            cout << V << "    " << newton(f,parameters,beta) << endl;
        }
        return 0;
    }

    vector<double> T;
    double T0;
    double T1;
    T0 = 1.88;
    T1 = 1.8846;
    int NT = 30;
    double dT = (T1 - T0)/NT;
    for (int i = 0; i <= NT; ++i)
    {
        T.push_back(T0+i*dT);
    }
    for (int i = 0; i < T.size(); ++i)
    {
        beta = 1.0/T[i];
        cout << T[i] << "    " << newton(f, p, beta) << endl;
    }
    return 0;
}
