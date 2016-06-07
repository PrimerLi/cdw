#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#define PI acos(-1.0)
using namespace std;

struct Parameters
{
    double g, K, V;
    Parameters(double g, double K, double V)
    {
        this->g = g;
        this->K = K;
        this->V = V;
    }
};

double square(double a)
{
    return a*a;
}

double Gaussian(double epsilon)
{
    return exp(-epsilon*epsilon)/sqrt(PI);
}

double integrand(const Parameters &p, double beta, double epsilon, double Delta)
{
    double g, K, V;
    g = p.g;
    K = p.K;
    V = p.V;
    double one, two;
    double result;
    one = (1/sqrt(epsilon*epsilon + Delta*Delta) - 1/sqrt(epsilon*epsilon + Delta*Delta + 4*V*V))*tanh(0.25*beta*(sqrt(epsilon*epsilon + Delta*Delta) - sqrt(epsilon*epsilon + Delta*Delta + 4*V*V)));
    two = (1/sqrt(epsilon*epsilon + Delta*Delta) + 1/sqrt(epsilon*epsilon + Delta*Delta + 4*V*V))*tanh(0.25*beta*(sqrt(epsilon*epsilon + Delta*Delta) + sqrt(epsilon*epsilon + Delta*Delta + 4*V*V)));
    result = 2*g*g/K*(one + two)*Gaussian(epsilon);
    return result;
}

double f(const Parameters &p, double beta, double Delta)
{
    double s;
    double upper;
    double lower;
    double denergy;
    int Nenergy, ienergy;
    double energy;

    upper = 8;
    lower = -upper;
    Nenergy = 800;
    denergy = (upper - lower)/Nenergy;

    s = 0;
    for (ienergy = 0; ienergy <= Nenergy; ++ienergy)
    {
        energy = lower + ienergy*denergy;
        s = s + denergy*integrand(p, beta, energy, Delta);
    }
    s = s - 1;
    return s;
}

double fprime(const Parameters &p, double beta, double Delta)
{
    double dDelta = 0.00001;
    return (f(p, beta, Delta + dDelta) - f(p, beta, Delta))/dDelta;
}

double newton(double (*f)(const Parameters &p, double beta, double Delta), const Parameters &p, double beta)
{
    double x0 = 3.0;
    double x1 = x0 - f(p, beta, x0)/fprime(p, beta, x0);
    int count = 0;
    int iterationMax = 8;
    while(fabs(x0 - x1) > 0.001)
    {
        count++;
        x0 = x1;
        x1 = x0 - f(p, beta, x0)/fprime(p, beta, x0);
        if (count > iterationMax)
        {
            cout << "No solution found. " << endl;
            return -1;
        }
    }
    return x1;
}

int main(int argc, char **argv)
{
    double g, K, V;
    double beta;
    g = 1;
    K = 1;
    V = 1;
    
    if (false)
    {
        double Delta;
        double x0, x1;
        int N = 100;
        x0 = 0;
        x1 = 6;
        double dx = (x1 - x0)/N;
        Parameters p(g,K,V);
        beta = 1000;
        ofstream ofile;
        ofile.open("f.txt");
        for (int i = 0; i <= N; ++i)
        {
            Delta = x0 + (i+0.5)*dx;
            ofile << Delta << "    " << f(p, beta, Delta) << endl;
        }
        ofile.close();
        return 0;
    }

    if (true)
    {
        V = 2.366;
        double T0, T1;
        double dT;
        int NT;
        vector<double> T;
        T0 = 0.951;
        T1 = 1.03;
        NT = 40;
        dT = (T1 - T0)/NT;
        for (int i = 0; i <= NT; ++i)
        {
            T.push_back(T0 + (i+0.5)*dT);
        }

        Parameters p(g, K, V);
        for (int i = 0; i <= NT; ++i)
        {
            beta = 1.0/T[i];
            cout << T[i] << "    " << newton(f, p, beta) << endl;
        }
    }
    return 0;
}
