// SOS.h - Champ Darabundit 09/28/21
// Second order system filter handling
#ifndef sos_h
#define sos_h

#include <Eigen/Dense>

template <typename T>
struct coefs{
    T b0, b1, b2, a0, a1, a2;
    coefs() {}
    coefs(T* c)
    {
        b0 = c[0]; b1 = c[1]; b2 = c[2];
        a0 = c[3]; a1 = c[4]; a2 = c[5];
    }
    
    // Normalize coefs so a0 = 1
    void norm(){
        b0 = b0/a0;
        b1 = b1/a0;
        b2 = b2/a0;
        a1 = a1/a0;
        a2 = a2/a0;
        a0 = static_cast<T>(1.0);
    }
};

// BIQUAD TEMPLATE - Base for SOS
template <typename T>
class biquad
{
public:
    biquad()
    {}
    ~biquad()
    {}
    
    void set_fs(double newfs)
    {
        fs = newfs;
    }
    // Process DF2T
    float process(float input)
    {
        float output = dc.b0*input + state[0];
        state[0] = static_cast<T>( dc.b1*input - dc.a1*output + state[1]);
        state[1] = static_cast<T>( dc.b2*input - dc.a2*output);
        return output;
    }
    
    void reset_dsp()
    {
        state[0] = static_cast<T>(0);
        state[1] = static_cast<T>(0);
    }
    // Setters for the digital coefs
    void setDigitalCoefs(coefs<T> newCoefs)
    {
        dc = newCoefs;
        dc.norm();
    }
    
    void setAnalogCoefs(coefs<T> newCoefs)
    {
        ac = newCoefs;
    }
    
    // Overwrite the digital coefs with bilinear transform analog coefs
    void bilinear()
    {
        T b0z, b1z, b2z, a0z, a1z, a2z;
        T c = 2.0*fs;
        // b0*(z^-2 + 2z^-1 + 1)
        // c*b1*(z^-2 - 1)
        // c*c*b2*(z^-2 - 2z^-1 + 1)
        b0z = ac.b0 + c*c*ac.b2 - c*ac.b1;
        b1z = 2*(ac.b0 - c*c*ac.b2);
        b2z = ac.b0 + c*ac.b1 + c*c*ac.b2;
        a0z = ac.a0 + c*c*ac.a2 - c*ac.a1;
        a1z = 2*(ac.a0 - c*c*ac.a2);
        a2z = ac.a0 + c*ac.a1 + c*c*ac.a2;
        T tmp[6] = {b0z, b1z, b2z, a0z, a1z, a2z};
        dc = new coefs<T> (tmp);
        dc.norm();
    }
    
    // Compute analog response at frequency f
    float freqs(T f)
    {
        std::complex<T> s(0, 2.0*M_PI*f);
        auto num = ac.b0 + ac.b1*s + ac.b2*s*s;
        auto den = ac.a0 + ac.a1*s + ac.a2*s*s;
        return std::abs(num/den);
    }
    
    float freqz(T f)
    {
        std::complex<T> j(0, 1);
        std::complex<T> z = std::exp( -j * 2.0*M_PI*(f/fs));
        auto num = dc.b0 + dc.b1*z + dc.b2*z*z;
        auto den = dc.a0 + dc.a1*z + dc.a2*z*z;
        return std::abs(num/den);
    }
    
private:
    double fs = 0.0;
    T state[2];
    coefs<T> ac; // Analog coefs
    coefs<T> dc; // Digital coefs
};

template <typename T, int numSections>
class sos{
public:
    using fiter = typename std::vector<biquad<T>>::iterator;
    using citer = typename std::vector<coefs<T>>::iterator;
    sos()
    {
        N = numSections;
        T c[6] = {1.f, 0.f, 0.f, 1.f, 0.f, 0.f};
        coefs<T> tmp(c);
        for (int i = 0; i < N; ++i)
        {
            filters.push_back(biquad<T>());
            // Default to pass through
            filters[i].setAnalogCoefs(tmp);
            filters[i].setDigitalCoefs(tmp);
        }
        // Write our iterators
        start = filters.begin();
        stop = filters.end();
    }
    
    ~sos() {}
    
    void set_fs(float newfs)
    {
        fs = newfs;
        for (fiter it = start; it != stop; ++it)
        {
            it->set_fs(fs);
        }
    }
    
    float process(float input)
    {
        float filtIn = input;
        float output;
        for (fiter it = start; it!= stop; ++it)
        {
            output = it->process(filtIn);
            filtIn = output; // Output of prev section is input to new section
        }
        return output;
    }
    
    float freqs(T f)
    {
        float mag = 1.0f;
        for (fiter it = start; it != stop; ++it)
        {
            mag *= it->freqs(f);
        }
        return mag;
    }
    
    float freqz(T f)
    {
        float mag = 1.0f;
        for (fiter it = start; it != stop; ++it)
        {
            mag *= it->freqz(f);
        }
        return mag;
    }
    
    float bilinear()
    {
        for (fiter it = start; it != stop; ++it)
            it->bilinear();
    }
    
    void setAnalogCoefs(coefs<T>* newCoefs, int sec)
    {
        filters[sec].setAnalogCoefs(*newCoefs);
    }
    
    void setDigitalCoefs(coefs<T>* newCoefs, int sec)
    {
        filters[sec].setDigitalCoefs(*newCoefs);
    }
    
//    // This asks you to set shape the coefs correctly
//    void setAnalogCoefs(citer startCoefs, citer stopCoefs)
//    {
//        fiter fIt = start;
//        for (citer cIt = startCoefs; cIt != stopCoefs; cIt++, fIt++)
//        {
//            fIt->setAnalogCoefs(*cIt);
//        }
//    }
    
    // Resize number of second order sections. Doesn't reallocate
    void resize(int newNSecs)
    {
        N = newNSecs;
        filters.resize(N);
        start = filters.begin();
        stop = filters.end();
    }
    
    void set_fs(double newfs)
    {
        fs = newfs;
        for (fiter it = start; it != stop; it++)
            it->set_fs(fs);
    }
private:
    int N;
    double fs = 0.0;
    std::vector<biquad<T>> filters;
    fiter start;
    fiter stop;
};
#endif /* sos_h */
