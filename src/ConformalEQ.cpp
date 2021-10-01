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
    
    treeState.addParameterListener(IDs::filterFreq, this);
    treeState.addParameterListener(IDs::filterQ, this);
    
    magicState.setGuiValueTree( BinaryData::gui_xml, BinaryData::gui_xmlSize);
}

void ConformalEQ::parameterChanged(const String &parameterID, float newValue)
{}

void ConformalEQ::prepareToPlay(double sampleRate, int samplesPerBlock)
{
    filterPlot->prepareToPlay(sampleRate,samplesPerBlock);
    
    magicState.prepareToPlay(sampleRate, samplesPerBlock);
    double tmp[6] = {3.947841760435743e7, 0.0, 0.0, 3.947841760435743e7, 8.885765876316731e3, 1.0};
    auto c = new coefs<double>(tmp);
    filter.setAnalogCoefs( c, 0);
    filter.set_fs(sampleRate);
    filter.bilinear();
    // attach the filter
    filterPlot->freqs(&filter, 0);
    filterPlot->freqz(&filter, 1);
}

void ConformalEQ::releaseResources()
{}

void ConformalEQ::processBlock(AudioBuffer<float> &buffer, MidiBuffer &midiMessages)
{}

void ConformalEQ::processAudioBlock(AudioBuffer<float>& buffer)
{}

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



