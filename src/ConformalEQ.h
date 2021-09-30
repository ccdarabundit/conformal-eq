// Conformal EQ - Champ Darabundit 09/27/21

#pragma once

#include <JuceHeader.h>
#include "gui/champLAF.h"
#include "dsp/sos.h"
#include "gui/FilterPlot.h"

class ConformalEQ : public foleys::MagicProcessor,
                    private juce::AudioProcessorValueTreeState::Listener
{
public:
    ConformalEQ(); 
	~ConformalEQ() override {}
    
    // Process functions
	void prepareToPlay (double sampleRate, int samplesPerBlock) override;
	void releaseResources() override;
    void processAudioBlock (AudioBuffer<float>& buffer);
    void processBlock ( AudioBuffer<float>&buffer, MidiBuffer &midiMessages);
    
    // Decay control
    double getTailLengthSeconds () const override           { return 0.0;       }
    
    // MIDI functions
    bool acceptsMidi () const override                      { return 0;         }
    bool producesMidi () const override                     { return 0;         }
    
    // Audio processor program functions
    int getNumPrograms () override                          { return 1;         }
    int getCurrentProgram () override                       { return 0;         }
    void setCurrentProgram (int index) override             {                   }
    const String getProgramName (int) override              { return "None";    }
    void changeProgramName (int, const String&) override    {                   }
    
    void getStateInformation (juce::MemoryBlock &destData) override
                                                            {                   }
    void setStateInformation (const void *data, int sizeInBytes) override
                                                            {                   }
    
    // PGM/Editor Functions
    AudioProcessorEditor * 	createEditor () override {
    	auto builder = std::make_unique<foleys::MagicGUIBuilder>(magicState);
    	builder->registerJUCEFactories();
    	builder->registerLookAndFeel("champLAF", std::make_unique<champLAF>());
    	auto editor = new foleys::MagicPluginEditor (magicState, std::move (builder));
    	return editor;
    }
    bool hasEditor () const override                        { return true;       }
    void parameterChanged (const String& parameterID, float newValue) override;
    static void addParameters (juce::AudioProcessorValueTreeState::Parameter& params);
    const String getName() const override           { return "Conformal EQ"; }
private: 
	AudioBuffer<float> monoBuffer;
    // PGM/Editor variables
    juce::AudioProcessorValueTreeState treeState;
    juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();
    foleys::MagicProcessorState magicState { *this };
    FilterPlot* filterPlot = nullptr;
    FilterPlot* filterPlotZ = nullptr;
    // Filter variables
    sos<double, 1> filter;
    std::atomic<float>* frequency = nullptr;
    std::atomic<float>* Q = nullptr;
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ConformalEQ)
};
