#include "FilterDesign.h"

FilterDesign::FilterDesign()
{
    system.resize(1);
    // Initialize to 1 kHz lowpass
    type = Lowpass;
    Q = 1.0/std::sqrt(2);
    F0 = 1000.0f;
    
    updateCoefs(type, Q, F0);

}

void FilterDesign::setType(FilterType t)
{
    type = t;
    updateCoefs(t, Q, F0);
    system.reset_dsp();
}

void FilterDesign::setQ(float q)
{
    Q = q;
    updateCoefs(type, q, F0);
}

void FilterDesign::setFreq(float freq)
{
    F0 = freq;
    updateCoefs(type, Q, freq);
}

sos<double, 1>* FilterDesign::getSystem()
{
    return &system;
}
void FilterDesign::prepareToPlay(double sampleRate)
{
    fs = sampleRate;
    system.set_fs(fs);
    system.sd();
}
void FilterDesign::updateCoefs(FilterType t, float q, float freq)
{
    std::vector<coefs<double>> cVec;
    coefs<double> analogCoefs;
    double wc = 2.0*M_PI*freq;
    system.set_warp(freq);
    switch(t)
    {
        case Lowpass:
            analogCoefs.b0 = wc*wc;
            analogCoefs.b1 = analogCoefs.b2 = 0;
            analogCoefs.a0 = wc*wc;
            analogCoefs.a1 = wc/q;
            analogCoefs.a2 = 1.0;
            system.setAnalogSectionCoefs(&analogCoefs, 0);
            break;
        case Highpass:
            analogCoefs.b0 = analogCoefs.b1 = 0;
            analogCoefs.b2 = 1;
            analogCoefs.a0 = wc*wc;
            analogCoefs.a1 = wc/q;
            analogCoefs.a2 = 1.0;
            system.setAnalogSectionCoefs(&analogCoefs, 0);
            break;
        case Bandpass:
            analogCoefs.b0 = analogCoefs.b2 = 0;
            analogCoefs.b1 = wc/q;
            analogCoefs.a0 = wc*wc;
            analogCoefs.a1 = wc/q;
            analogCoefs.a2 = 1.0;
            system.setAnalogSectionCoefs(&analogCoefs, 0);
            break;
        case Notch:
            analogCoefs.b1 = 0;
            analogCoefs.b2 = 1;
            analogCoefs.b0 = -wc*wc;
            analogCoefs.a0 = wc*wc;
            analogCoefs.a1 = wc/q;
            analogCoefs.a2 = 1.0;
            system.setAnalogSectionCoefs(&analogCoefs, 0);
            break;
    }
    // system.bilinear(); if you want raw bilinear
    system.sd();
}
