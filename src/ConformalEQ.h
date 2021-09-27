// Conformal EQ - Champ Darabundit 09/27/21

#pragma once

#include <JuceHeader.h>

class ConformalEQ : public AudioProcessor
{
public: 
	ConformalEQ() {}

    static void addParameters (juce::AudioProcessorValueTreeState::Parameter& params) {}
	void prepareToPlay (double sampleRate, int samplesPerBlock) override;
	void releaseResources() override;
    void processAudioBlock (AudioBuffer<float>& buffer);
private: 
	AudioBuffer<float> monoBuffer;

	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ConformalEQ)
};
