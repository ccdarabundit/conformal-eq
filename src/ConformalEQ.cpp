// Conformal EQ - Champ Darabundit 09/27/21

#include "ConformalEQ.h"
#include "gui/Figure.h"
namespace IDs
{
static juce::String filterQ {"filterQ"};
static juce::String filterFreq
{"filterFreq"};
static juce::String filterType{"filterType"};
static juce::Identifier  filterPlot("filterPlot");
static juce::Identifier filterPlotZ("filterPlotZ");
}

juce::AudioProcessorValueTreeState::ParameterLayout ConformalEQ::createParameterLayout()
{
    juce::AudioProcessorValueTreeState::ParameterLayout layout;
    
    // FILTER PARAMETERS
    auto filter = std::make_unique<juce::AudioProcessorParameterGroup>("filter", "Filters", "|");
    filter->addChild(std::make_unique<juce::AudioParameterChoice>(IDs::filterType, "Type", juce::StringArray(
                                                                                                            "LP", "HP", "BP", "N"), 0));
    filter->addChild(std::make_unique<juce::AudioParameterFloat>(IDs::filterFreq, "Frequency",
                                                                 juce::NormalisableRange<float>(20.0f, 48000.0f, (48000.0-20.0)/(1 << 10)), 1000.0f));
    filter->addChild(std::make_unique<juce::AudioParameterFloat>(IDs::filterQ, "Q",
                                                                 juce::NormalisableRange<float>(0.25f, 10.0f, (10.0f-0.25f)/(1 << 7)), std::sqrt(2)));
    layout.add( std::move (filter)
               );
    return layout;
}

ConformalEQ::ConformalEQ()
#ifndef JucePlugin_PreferredChannelConfigurations
: MagicProcessor (BusesProperties()
                #if ! JucePlugin_IsMidiEffect
                 #if ! JucePlugin_IsSynth
                  .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                 #endif
                  .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                #endif
                  ),
#else
:
#endif
treeState (*this, nullptr, JucePlugin_Name, createParameterLayout())
{
    // Check that everything is hooked up
    frequency = treeState.getRawParameterValue(IDs::filterFreq);
    jassert( frequency != nullptr);
    Q = treeState.getRawParameterValue(IDs::filterQ);
    jassert( Q != nullptr);
    filterType = treeState.getRawParameterValue(IDs::filterType);
    jassert( filterType != nullptr);
    
    treeState.addParameterListener(IDs::filterFreq, this);
    treeState.addParameterListener(IDs::filterQ, this);
    treeState.addParameterListener(IDs::filterType, this);
    magicState.setGuiValueTree( BinaryData::gui_xml, BinaryData::gui_xmlSize);
}

void ConformalEQ::parameterChanged(const String &parameterID, float newValue)
{
    if (parameterID == IDs::filterQ)
        filter.setQ(newValue);
    else if(parameterID == IDs::filterFreq)
        filter.setFreq(newValue);
    else if(parameterID == IDs::filterType)
    {
        // Convert range to filter type
        FilterType newType = static_cast<FilterType>(static_cast<int>(newValue));
        filter.setType(newType);
    }
    
    // Replot
    filterPlot->freqs(filter.getSystem(), 0);
    filterPlot->freqz(filter.getSystem(), 1);
}

void ConformalEQ::prepareToPlay(double sampleRate, int samplesPerBlock)
{
    filterPlot->prepareToPlay(sampleRate,samplesPerBlock);
    
    magicState.prepareToPlay(sampleRate, samplesPerBlock);
    
    filter.prepareToPlay(sampleRate);
    // attach the filter
    filterPlot->freqs(filter.getSystem(), 0);
    filterPlot->freqz(filter.getSystem(), 1);
}

void ConformalEQ::releaseResources()
{}

void ConformalEQ::processBlock(AudioBuffer<float> &buffer, MidiBuffer &midiMessages)
{ processAudioBlock(buffer);}

void ConformalEQ::processAudioBlock(AudioBuffer<float>& buffer)
{
    const auto numChannels = buffer.getNumChannels();
    const auto numSamples = buffer.getNumSamples();
    // sum input to mono
    monoBuffer.setSize (1, numSamples, false, false, true);
    monoBuffer.clear();

    for (int samp = 0; samp < numSamples; ++samp)
    {
        double wnoise = random.nextFloat()*.25 - .125;
        float output = filter.getSystem()->process(wnoise);
        monoBuffer.setSample(0, samp, output);
    }
    // copy back to stereo (or however many channels)
        for (int ch = 0; ch < numChannels; ++ch)
            FloatVectorOperations::copy (buffer.getWritePointer (ch),
                                         monoBuffer.getReadPointer (0),
                                         numSamples);
}

AudioProcessorEditor *     ConformalEQ::createEditor ()
{
    auto builder = std::make_unique<foleys::MagicGUIBuilder>(magicState);
    builder->registerJUCEFactories();
    builder->registerLookAndFeel("champLAF", std::make_unique<champLAF>());
    builder->registerFactory("Figure", &FigureItem::factory);
    auto editor = new foleys::MagicPluginEditor (magicState, std::move (builder));
    return editor;
}

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new ConformalEQ();
}



