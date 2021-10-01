#pragma once

#include <JuceHeader.h>
#include "../dsp/sos.h"

// Custom filter plot using own sos filter structure

class FilterPlot
{
public:
    FilterPlot()
    {
        frequencies.resize (300);
            for (size_t i = 0; i < frequencies.size(); ++i)
                frequencies [i] = 20.0 * std::pow ((48000.0f/20.0f), (float) i/300.0 );

    }
    ~FilterPlot() {}
    
    template <typename T, int nSecs>
    void freqs(sos<T, nSecs> * system, int ind)
    {
        if (ind + 1 > magnitudes.size())
        {
            N += 1;
            std::vector<double>mag(300);
            magnitudes.push_back( mag );
        }
        for (int i = 0; i < static_cast<int>(frequencies.size()); i++)
            magnitudes[ind][i] = system->freqs(frequencies[i]);
    }
    
    template <typename T, int nSecs>
    void freqz(sos<T, nSecs> * system, int ind)
    {
        if (ind + 1 > magnitudes.size())
        {
            N += 1;
            std::vector<double>mag(300);
            magnitudes.push_back( mag );
        }
        for (int i = 0; i < static_cast<int>(frequencies.size()); i++)
            magnitudes[ind][i] = system->freqz(frequencies[i]);
    }
    
    std::vector<double> get_frequencies()
    {
        return frequencies;
    }
    
    void prepareToPlay (double sampleRate, int samplesPerBlockExpected)
    {
        fs = sampleRate;
    }
    
    int getN()
    {
        return N;
    }
    
    std::vector<double> get_mag(int ind)
    {
        return magnitudes[ind];
    }
private:
    int N = 0;
    juce::ReadWriteLock plotLock;
    std::vector<double> frequencies;
    std::vector<std::vector<double>> magnitudes;
    double maxDB = 30.0f;
    double fs = 0.0;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (FilterPlot)
};
