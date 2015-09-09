#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
using namespace std;

class Vector
{
    private:
        int dimension;
        double *vector;
        void checkBound(int i) const
        {
            if (i < 0 || i >= dimension)
            {
                cout << "Index out of range. " << endl;
                exit(-1);
            }
        }
    public:
        Vector()
        {
            this->dimension = 2;
            vector = new double [dimension];
            for (int i = 0; i < dimension; ++i)
            {
                vector[i] = 0;
            }
        }
        explicit Vector(int dimension)
        {
            this->dimension = dimension;
            vector = new double [dimension];
            for (int i = 0; i < dimension; ++i)
            {
                vector[i] = 0;
            }
        }
        Vector(const Vector &parameter)
        {
            this->dimension = parameter.dimension;
            vector = new double [dimension];
        }
        ~Vector()
        {
            delete []vector;
        }
        double & operator[] (int i)
        {
            this->checkBound(i);
            return vector[i];
        }
        double operator[] (int i) const
        {
            this->checkBound(i);
            return vector[i];
        }
        int getDimension() const
        {
            return this->dimension;
        }
        friend ostream & operator<< (ostream &os, const Vector &parameter)
        {
            int dimension = parameter.dimension;
            for (int i = 0; i < dimension; ++i)
            {
                os << parameter[i] << "  ";
            }
            os << endl;
            return os;
        }
};

class Matrix
{
    private:
        int dimension;
        double **matrix;
        void checkBound(int i) const
        {
            if (i < 0 || i >= dimension)
            {
                cout << "Index out of bounds exception occurred. " << endl;
                exit(-1);
            }
        }
    public:
        Matrix()
        {
            dimension = 2;
            matrix = new double *[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                matrix[i] = new double [dimension];
            }
            for (int i = 0; i < dimension; ++i)
            {
                for (int j = 0; j < dimension; ++j)
                {
                    matrix[i][j] = 0;
                }
            }
        }
        explicit Matrix(int dimension)
        {
            this->dimension = dimension;
            matrix = new double *[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                matrix[i] = new double [dimension];
            }
            for (int i = 0; i < dimension; ++i)
            {
                for (int j = 0; j < dimension; ++j)
                {
                    matrix[i][j] = 0;
                }
            }
        }
        Matrix(const Matrix &parameter)
        {
            dimension = parameter.dimension;
            matrix = new double *[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                matrix[i] = new double [dimension];
            }
            for (int i = 0; i < dimension; ++i)
            {
                for (int j = 0; j < dimension; ++j)
                {
                    matrix[i][j] = parameter.matrix[i][j];
                }
            }
        }
        ~Matrix()
        {
            for (int i = 0; i < dimension; ++i)
            {
                delete []matrix[i];
            }
            delete []matrix;
        }
        double det() const
        {
            return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
        }
        Matrix inverse() const
        {
            if (this->det() == 0) 
            {
                cout << "This matrix is singular. " << endl;
                exit(-1);
            }
            Matrix result(this->dimension);
            result.matrix[0][0] = matrix[1][1]/det();
            result.matrix[0][1] = -matrix[0][1]/det();
            result.matrix[1][0] = -matrix[1][0]/det();
            result.matrix[1][1] = matrix[0][0]/det();
            return result;
        }
        const Matrix& operator= (const Matrix &parameter) const
        {
            if (this == &parameter) return *this;
            if (this->dimension != parameter.dimension)
            {
                cout << "Matrix dimension wrong. " << endl;
                exit(-1);
            }
            for (int i = 0; i < dimension; ++i)
            {
                for (int j = 0; j < dimension; ++j)
                {
                    this->matrix[i][j] = parameter.matrix[i][j];
                }
            }
            return *this;
        }
        double & operator() (int i, int j)
        {
            this->checkBound(i);
            this->checkBound(j);
            return matrix[i][j];
        }
        double operator() (int i, int j)const
        {
            this->checkBound(i);
            this->checkBound(j);
            return matrix[i][j];
        }
        Vector operator* (const Vector &parameter) const
        {
            if (parameter.getDimension() != dimension)
            {
                cout << "Dimension wrong. " << endl;
                exit(-1);
            }
            Vector result(dimension);
            double s = 0;
            for (int i = 0; i < dimension; ++i)
            {
                s = 0;
                for (int j = 0; j < dimension; ++j)
                {
                    s = s + matrix[i][j]*parameter[j];
                }
                result[i] = s;
            }
            return result;
        }
        friend ostream & operator<< (ostream &os, const Matrix &parameter)
        {
            int dimension = parameter.dimension;
            for (int i = 0; i < dimension; ++i)
            {
                for (int j = 0; j < dimension; ++j)
                {
                    os << parameter(i, j) << "  ";
                }
                os << endl;
            }
            os << endl;
            return os;
        }
};

double sum(const vector<double> &array)
{
    double s = 0;
    for (int i = 0; i < array.size(); ++i)
    {
        s = s + array[i];
    }
    return s;
}

double square(const vector<double> &array)
{
    double s = 0;
    for (int i = 0; i < array.size(); ++i)
    {
        s = s + array[i]*array[i];
    }
    return s;
}

double inner(const vector<double> &A, const vector<double> &B)
{
    if (A.size() != B.size())
    {
        cout << "Dimension wrong. " << endl;
        exit(-1);
    }
    double s = 0;
    for (int i = 0; i < A.size(); ++i)
    {
        s = s + A[i]*B[i];
    }
    return s;
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cout << "inputFile = argv[1]. " << endl;
        exit(-1);
    }

    vector<double> x;
    vector<double> y;
    double a, b;
    
    string inputFile = string(argv[1]);
    string line;
    int count = 0;
    ifstream ifile;
    ifile.open(inputFile.c_str());
    while(!ifile.eof())
    {
        getline(ifile, line);
        if (line != "")
        {
            count++;
        }
    }
    ifile.close();

    ifile.open(inputFile);
    double lambda, Tc;
    for (int i = 0; i < count; ++i)
    {
        ifile >> lambda >> Tc;
        //cout << lambda << "    " << Tc << endl;
        x.push_back(1.0/(lambda*lambda));
        y.push_back(log(Tc));
    }
    ifile.close();

    int dimension = 2;
    Matrix matrix(dimension);
    matrix(0, 0) = count;
    matrix(0, 1) = -sum(x);
    matrix(1, 0) = sum(x);
    matrix(1, 1) = -square(x);

    Vector vector(dimension);
    vector[0] = sum(y);
    vector[1] = inner(x, y);

    Vector result = matrix.inverse()*vector;

    a = exp(result[0]);
    b = result[1];
    cout << "a = " << a << ", b = " << b << endl;
    return 0;
}
