// Handles coefficients for sos
#ifndef coefs_h
#define coefs_h

#include <Eigen/Dense>

template <typename T>
struct coefs{
    T b0, b1, b2, a0, a1, a2;
    Eigen::Matrix<T, 3, 1> bcoefs;
    Eigen::Matrix<T, 3, 1> acoefs;
    coefs() {}
    coefs(T* c)
    {
        b0 = c[0]; b1 = c[1]; b2 = c[2];
        a0 = c[3]; a1 = c[4]; a2 = c[5];
        updateVecs();
    }

    // Normalize coefs so a0 = 1
    void norm(){
        b0 = b0/a0;
        b1 = b1/a0;
        b2 = b2/a0;
        a1 = a1/a0;
        a2 = a2/a0;
        a0 = static_cast<T>(1.0);
        updateVecs();
    }
    
    void set(T* c)
    {
        b0 = c[0]; b1 = c[1]; b2 = c[2];
        a0 = c[3]; a1 = c[4]; a2 = c[5];
        updateVecs();
    }
    
    void set(Eigen::Matrix<T, 3, 1> newB, Eigen::Matrix<T, 3, 1> newA)
    {
        b0 = newB[2]; b1 = newB[1]; b2 = newB[0];
        a0 = newA[2]; a1 = newA[1]; a2 = newA[0];
        bcoefs = newB; acoefs = newA;
    }
    
    void updateVecs()
    {
        bcoefs << b2, b1, b0;
        acoefs << a2, a1, a0;
    }
    
    void operator = (const coefs<T> &c)
    {
        b0 = c.b0; b1 = c.b1; b2 = c.b2;
        a0 = c.a0; a1 = c.a1; a2 = c.a2;
        updateVecs();
    }

};

#endif /* coefs_h */
