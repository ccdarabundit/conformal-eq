// FilterDesign.h - Champ Darabundit
// Filter type and structure for building filters are here
#pragma once
#include "sos.h"
#include "SpectralDiscretize.h"

enum FilterType{
    Lowpass = 0,
    Highpass,
    Bandpass,
    Notch
};
class FilterDesign
{
public:
    FilterDesign();
    void setType(FilterType t);
    void setQ(float q);
    void setFreq(float freq);
    void updateCoefs(FilterType t, float q, float freq);
    void prepareToPlay(double sampleRate);
    sos<double, 1>* getSystem();
private:
    FilterType type;
    double fs;
    float Q;
    float F0;
    sos<double, 1> system;
};
