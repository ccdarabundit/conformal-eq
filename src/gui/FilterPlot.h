#include <JuceHeader.h>
#include "../dsp/sos.h"

// Custom filter plot using own sos filter structure

class FilterPlot : public foleys::MagicPlotSource
{
public:
    FilterPlot()
    {
        frequencies.resize (300);
            for (size_t i = 0; i < frequencies.size(); ++i)
                frequencies [i] = 20.0 * std::pow (2.0, i / 30.0);

        magnitudes.resize (frequencies.size());
    }
    ~FilterPlot() {}
    // This does nothing
    void pushSamples(const juce::AudioBuffer< float > &buffer)
    {}
    
    template <typename T, int nSecs>
    void freqs(sos<T, nSecs> * system)
    {
        for (int i = 0; i < static_cast<int>(frequencies.size()); i++)
            magnitudes[i] = system->freqs(frequencies[i]);
        
        resetLastDataFlag();
    }
    
    template <typename T, int nSecs>
    void freqz(sos<T, nSecs> * system)
    {
        for (int i = 0; i < static_cast<int>(frequencies.size()); i++)
            magnitudes[i] = system->freqz(frequencies[i]);
        
        resetLastDataFlag();
    }
    
    void createPlotPaths (juce::Path &path, juce::Path &filledPath, juce::Rectangle< float > bounds, foleys::MagicPlotComponent &component)
    {
        const juce::ScopedReadLock readLock (plotLock);
        
        const auto yFactor = 2.0f * bounds.getHeight() / juce::Decibels::decibelsToGain (maxDB);
        const auto xFactor = static_cast<double> (bounds.getWidth())/ frequencies.size();
        
        path.clear();
        path.startNewSubPath (bounds.getX(), float (magnitudes [0] > 0 ? bounds.getCentreY() - yFactor * std::log (magnitudes [0]) / std::log (2) : bounds.getBottom()));
            for (size_t i=1; i < frequencies.size(); ++i)
                path.lineTo (float (bounds.getX() + i * xFactor),
                             float (magnitudes [i] > 0 ? bounds.getCentreY() - yFactor * std::log (magnitudes [i]) / std::log (2) : bounds.getBottom()));

            filledPath = path;
            filledPath.lineTo (bounds.getBottomRight());
            filledPath.lineTo (bounds.getBottomLeft());
            filledPath.closeSubPath();
    }
    void prepareToPlay (double sampleRate, int samplesPerBlockExpected)
    {
        fs = sampleRate;
    }
    
private:
    juce::ReadWriteLock plotLock;
    std::vector<double> frequencies;
    std::vector<double> magnitudes;
    double maxDB = 30.0f;
    double fs = 0.0;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (FilterPlot)
};
