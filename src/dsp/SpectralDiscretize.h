// SpectralDiscretize.h - Champ Darabundit 10/1/21
// Has all the spectral discretization helpers

#pragma once

#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/Polynomials>
#include "coefs.h"

// Make this any order too

template<typename T, int order>
class SpectralDiscretize
{
public:
    SpectralDiscretize()
    {
        // Initialize constant matrices
        // G2, Dp, and invDs
        for (int i = 0; i < N + 1; i ++)
        {
            G2(i) = std::pow(2, N - i);
            invG(i) = 1.0;
        }
        
        Dp = Eigen::MatrixXd::Constant(2*N+1, N+1, 0.0);
        bandstopCreate(&Dp, -1.0, -1.0);
        
        Eigen::Matrix<double, 2*N+1, N+1> tmp;
        tmp = Eigen::MatrixXd::Constant(2*N+1, N+1, 0.0);
        bandstopCreate(&tmp, 1.0, 1.0);
        invDs = tmp.completeOrthogonalDecomposition().pseudoInverse();
        std::cout << invDs << std::endl;
    }
    
    ~SpectralDiscretize()
    {}
    
    void set_fs(double samplingRate)
    {
        fs = samplingRate;
        // Set F and invF;
        for (int i = 0; i < N+1; i++)
        {
            double powFs = std::pow(M_PI*fs, N - i);
            F(i) = powFs;
            invF(i)  = 1.0/powFs;
         }
    }
    
    void discretize(coefs<T>* inCoefs)
    {
        Eigen::Matrix<double, 2*N+1, 1> bp, ap;

        // 1: Forward transform
        bp = Dp * Eigen::MatrixXd(G2.asDiagonal()) * Eigen::MatrixXd(F.asDiagonal())*inCoefs->bcoefs;
        ap = Dp * Eigen::MatrixXd(G2.asDiagonal()) * Eigen::MatrixXd(F.asDiagonal())*inCoefs->acoefs;
        // Normalize
        bp = bp/ap(2*N);
        ap = ap/ap(2*N);
        
        // 2: Stabilize
        if (N == 2)
        {
            stabilize_quartic(&bp);
            stabilize_quartic(&ap);
            
        }
        else
        {
            stabilize_sys(&bp, &ap);
        }
        // Normalize
        bp = bp/ap(2*N);
        ap = ap/ap(2*N);
        
        // 3: Inverse transform
        Eigen::Matrix<double, N+1, 1> bc = Eigen::MatrixXd(invF.asDiagonal())*Eigen::MatrixXd(invG.asDiagonal())*(invDs*bp);
        Eigen::Matrix<double, N+1, 1> ac = Eigen::MatrixXd(invF.asDiagonal())*Eigen::MatrixXd(invG.asDiagonal())*(invDs*ap);
        inCoefs->set(bc, ac);
    }
    
    void setG(double g1)
    {
        for (int i = 0; i < N+1; i++)
            invG(i) = 1.0/std::pow(g1, N - i);
    }
    
private:
    double fs = 1.0;
    static const int N = order;
    // Forward Matricies
    Eigen::Matrix<double, N+1, 1> F; // Diagonal frequency matrix represented as vector
    Eigen::Matrix<double, N+1, 1> G2;
    Eigen::Matrix<double, 2*N+1, N+1> Dp;
    
    // Backward matricies
    Eigen::Matrix<double, N+1, 1> invF; // Diagonal frequency matrix represented as vector
    Eigen::Matrix<double, N+1, 1> invG;
    Eigen::Matrix<double, N+1, 2*N+1> invDs;
    
    // Helper private functions
    double binom( int n, int k)
    {
        double res = 1.0;
        for (int i = 1; i <= k; ++i)
            res = res * (n - k + i) / i;
        return static_cast<double>( (int)(res + 0.01) ); // go to int then to double
    }
    
    // Create the bandstop matricies
    void bandstopCreate(Eigen::Matrix<double, 2*N+1, N+1>* A, double g2, double g3)
    {
        for (int j = 0; j < N + 1; j++)
        {
            for (int i = 0; i < N + 1; i++)
            {
                int r_ind = j + 2*i;
                if (r_ind < 2*N+1)
                    (*A)(2*N-r_ind,N-j) = std::pow(g2, j)*binom(N-j, i)*std::pow(g3, (N-j-i));
                else
                    break;
            }
        }
        std::cout << * A << std::endl;
    }
    
    // Iteratively compute filter coefficients
    double * vieta(std::vector<std::complex<double>> roots, unsigned long length)
    {
        std::complex<double> * coefs;
        coefs = new std::complex<double> [length];
        memset( coefs, 0, sizeof(*coefs));
        unsigned int i = 1;
        unsigned long rSize = static_cast<unsigned long>(roots.size());
        
        unsigned int rInd = 1;
        if (rSize+1 != length){
            rInd = (unsigned int) length - rSize;
            coefs[rInd-1] = 1.0;
        }
        else
        {
            coefs[0] = 1.0;
        }
        for (i = rInd-1; i < length-1; i++)
        {
            coefs[i+1] = -roots[i]*coefs[i];
            for (unsigned int j = i; j > 0; j--)
                coefs[j] = coefs[j] - coefs[j-1]*roots[i];
        }
        double * ret;
        ret = new double [length];
        for (i = 0; i < length; i++)
            ret[i] = static_cast<double>( coefs[i].real() );
        return ret;
    }
    
    
    int stabilize_quartic(Eigen::Matrix<double, 2*N+1, 1> * coefs)
    {
        // Make use of fact that resulting polynomial is a quasi-palindromic
        // quartic
        std::complex<double> yr1, yr2, r1, r2, r3, r4;
        std::complex<double> a1, b1, c1, d1, b2, c2, d2;
        if ( (*coefs)(1) == 0.0)
        {
            if ((*coefs)(0) == 0.0)
                return 0; // All roots zero, no need to stabilize
            
            a1 = (*coefs)(0); b1 = (*coefs)(2); c1 = (*coefs)(4);
            b1 = b1/a1*0.5;
            c1 = c1/a1;
            d1 = std::sqrt(b1*b1 - c1);
            yr1 = -b1 + d1;
            yr2 = -b1 - d1;
            r1 = std::sqrt(yr1);
            r2 = -r1;
            r3 = std::sqrt(yr2);
            r4 = -r3;
        }
        else
        {
            if ((*coefs)(0) == 0.0 && (*coefs)(2) == 0.0)
            {
                double * newCoefs;
                newCoefs = new double [2*N+1];
                memset(newCoefs, 0.0, sizeof(double)*(2*N+1));
                newCoefs[1] = 1.0;
                newCoefs[2] = 2.0;
                newCoefs[3] = 1.0;
                double k = std::abs((*coefs)(1));
                *coefs = k*Eigen::Map<Eigen::Matrix<double, 2*N+1, 1>>(newCoefs);
                return 0;
            }
            a1  = (*coefs)(0); b1 = (*coefs)(1); c1 = (*coefs)(2) - (2.0*(*coefs)(0)*(*coefs)(3))/(*coefs)(1);
            b1 = b1/a1*0.5;
            c1 = c1/a1;
            d1 = std::sqrt(b1*b1 - c1);
            yr1 = -b1 + d1;
            yr2 = -b1 - d1;
            // Second set of quadratics
            b1 = -yr1*0.5; c1 = (*coefs)(3)/(*coefs)(1);
            b2 = -yr2*0.5; c2 = (*coefs)(3)/(*coefs)(1);
            d1 = std::sqrt(b1*b1 - c1);
            d2 = std::sqrt(b2*b2 - c2);
            r1 = -b1 + d1;
            r2 = -b1 - d1;
            r3 = -b2 + d2;
            r4 = -b2 - d2;
        }
        std::vector<std::complex<double>> p;
        p.resize(4);

        // Stabilize
        p[0] = std::complex<double>( -1.0*std::abs(r1.real()), r1.imag());
        p[1] = std::complex<double>( -1.0*std::abs(r2.real()), r2.imag());
        p[2] = std::complex<double>( -1.0*std::abs(r3.real()), r3.imag());
        p[3] = std::complex<double>( -1.0*std::abs(r4.real()), r4.imag());
        
        // Now compute new coeffs
        double * newCoefs = new double [2*N+1];
        newCoefs[0] = (p[0]*p[1]*p[2]*p[3]).real();
        newCoefs[1] = (- p[0]*p[1]*p[2] - p[0]*p[1]*p[3] - p[0]*p[2]*p[3] - p[1]*p[2]*p[3]).real();
        newCoefs[2] = (p[0]*p[1] + p[0]*p[2] + p[0]*p[3] + p[1]*p[2] + p[1]*p[3] + p[2]*p[3]).real();
        newCoefs[3] = (- p[0] - p[1] - p[2] - p[3]).real();
        newCoefs[4] = 1.0;
        *coefs = Eigen::Map<Eigen::Matrix<double, 2*N+1, 1>>(newCoefs);
        return 1;
    }
    
    // Companion matrix method for root stabilization
    int stabilize_sys(Eigen::Matrix<double, 2*N+1, 1>*b, Eigen::Matrix<double, 2*N+1, 1>*a)
    {
        int aCnt, bCnt;
        aCnt = bCnt = -1;
        // Get true order of polynomial
        for (int i = 0; i < 2*N+1; i++)
        {
            if((*b)(i) != 0.0 && bCnt == -1)
                bCnt = i;
            if ((*a)(i) != 0.0 && aCnt == -1)
                aCnt = i;
            if( aCnt != -1 && bCnt != -1)
                break;
        }
        // Write to an Eigen vector for roots finding
        if (aCnt == -1 and bCnt == -1)
            return 0;  // Arrays are all zeros so we shouldn't procede
        // First we normalize
        (*b) = (*b)/((*a)(2*N));
        (*a) = (*a)/((*a)(2*N));
        // Then we get k
        double k = (*b)(2*N-bCnt);
        (*b) = (*b)/k;
        Eigen::VectorXd acoefs(2*N+1 - aCnt);
        Eigen::VectorXd bcoefs(2*N+1 - bCnt);
        // Solve numerator
        acoefs = ((*a).segment(aCnt, 2*N+1-aCnt)).reverse();
        Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
        solver.compute( acoefs );
        auto r1 = solver.roots();
        std::vector<std::complex<double>> p;
        if (r1.rows() > 0)
        {
            for (unsigned int i = 0; i < r1.rows(); i++)
            {
                p.push_back(r1(i));
                if (p[i].real() > 0.0)
                    p[i] = {-1.0*p[i].real(), p[i].imag()};
            }
        }
        else
            p.push_back(0);
        
        // Solve denominator
        bcoefs = ((*b).segment(bCnt, 2*N+1-bCnt)).reverse();
        solver.compute(  bcoefs );
        auto r2 = solver.roots();
        std::vector<std::complex<double>> z;
        if (r2.rows() > 0)
        {
            for (unsigned int i = 0; i < r2.rows(); i++)
            {
                z.push_back(r2(i));
                if (z[i].real() > 0.0)
                    z[i] = {-1.0*z[i].real(), z[i].imag()};
            }
        }
        else
            z.push_back(0);
        
        // Compute coefs
        double * newACoefs = vieta(p, 2*N+1);
        double * newBCoefs = vieta(z, 2*N+1);
        *a = Eigen::Map<Eigen::Matrix<double, 2*N+1, 1>>(newACoefs);
        *b = Eigen::Map<Eigen::Matrix<double, 2*N+1, 1>>(newBCoefs);
        (*b) *= k;
        
    }
};
