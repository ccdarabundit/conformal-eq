// SOS.h - Champ Darabundit 09/28/21
// Second order system filter handling
#ifndef sos_h
#define sos_h

#include <limits>

#include "coefs.h"
#include "SpectralDiscretize.h"

// BIQUAD TEMPLATE - Base for SOS
template <typename T>
class biquad
{
public:
    biquad()
    {reset_dsp();}
    ~biquad()
    {}
    
    void set_fs(double newfs)
    {
        fs = newfs;
//        digitizer.prepareToPlay(fs);
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
        b0z = ac.b0 + c*c*ac.b2 + c*ac.b1;
        b1z = 2*(ac.b0 - c*c*ac.b2);
        b2z = ac.b0 - c*ac.b1 + c*c*ac.b2;
        a0z = ac.a0 + c*c*ac.a2 + c*ac.a1;
        a1z = 2*(ac.a0 - c*c*ac.a2);
        a2z = ac.a0 - c*ac.a1 + c*c*ac.a2;
        T tmp[6] = {b0z/a0z, b1z/a0z, b2z/a0z, 1.0, a1z/a0z, a2z/a0z};
        dc.set(tmp);
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
    
    // Compute digital response at f
    float freqz(T f)
    {
        if (f > 0.5*fs) // Don't plot above Nyquist
            return NAN;
        else
        {
        std::complex<T> j(0, 1);
        std::complex<T> z = std::exp( -j * 2.0*M_PI*(f/fs));
        auto num = dc.b0 + dc.b1*z + dc.b2*z*z;
        auto den = dc.a0 + dc.a1*z + dc.a2*z*z;
        return std::abs(num/den);
        }
    }
    
    coefs<T> * getAnalogPtr()
    {
        return &ac;
    }
    
    coefs<T> * getDigitalPtr()
    {
        return &dc;
    }
private:
    double fs = 0.0;
    T state[2];
    coefs<T> ac; // Analog coefs
    coefs<T> dc; // Digital coefs
};

template <typename T, int maxNumSections>
class sos{
public:
    using fiter = typename std::vector<biquad<T>>::iterator;
    using citer = typename std::vector<coefs<T>>::iterator;
    sos()
    {
        N = 0;
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
    void reset_dsp()
    {
        for (fiter it = start; it != stop; ++it)
        {
            it->reset_dsp();
        }
    }
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
    
    void sd()
    {
        for (fiter it = start; it != stop; ++it)
        {
            biquad<T> * tmp = new biquad<T>();
            *tmp = *it;
            digitizer.discretize(tmp->getAnalogPtr());
            tmp->bilinear();
            it->setDigitalCoefs(*tmp->getDigitalPtr());
        }
    }
    
    void setAnalogSectionCoefs(coefs<T>* newCoefs, int sec)
    {
        filters[sec].setAnalogCoefs(*newCoefs);
    }
    
    void setDigitalSectionCoefs(coefs<T>* newCoefs, int sec)
    {
        filters[sec].setDigitalCoefs(*newCoefs);
    }
    
    // This asks you to set shape the coefs correctly
    void setAnalogCoefs(citer startCoefs, citer stopCoefs)
    {
        int nSecs = std::distance(startCoefs, stopCoefs);
        if (filters.size() < nSecs)
        {
            filters.resize(nSecs);
            N = nSecs;
            start = filters.begin();
            stop = filters.end();
        }
        else
        {
            N = nSecs; // No resizing we just don't process
            stop = std::next(start, N);
        }
        fiter fIt = start;
        
        for (citer cIt = startCoefs; cIt != stopCoefs; cIt++, fIt++)
        {
            fIt->setAnalogCoefs(*cIt);
        }
    }

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
        digitizer.set_fs(newfs);
        for (fiter it = start; it != stop; it++)
            it->set_fs(fs);
        calc_warp();
    }
    
    void set_warp(double newf0)
    {
        f0 = newf0;
        double beta, g;
        if (newf0 >= 0.5*fs)
            g = 2.333; // Optimized g
        else
        {
            beta = newf0;
            double bnorm = 2*beta/fs;
            double bh = (2 * fs * std::tan( M_PI*beta/fs ))/(M_PI*fs);
            g = 2.0*bh*std::sqrt(1 - bnorm*bnorm) /bnorm;
        }
        digitizer.setG( g );
    }
    
    void calc_warp()
    {
        double beta, g;
        if (f0 >= 0.5*fs)
            g = 2.333; // Optimized g
        else
        {
            beta = f0;
            double bnorm = 2*beta/fs;
            double bh = (2 * fs * std::tan( M_PI*beta/fs ))/(M_PI*fs);
            g = 2.0*bh*std::sqrt(1 - bnorm*bnorm) /bnorm;
        }
        digitizer.setG( g );
    }
    
private:
    int N;
    double fs = 1.0;
    double f0 = 0.0;
    std::vector<biquad<T>> filters;
    fiter start;
    fiter stop;
    SpectralDiscretize<T, 2> digitizer;
};
#endif /* sos_h */
