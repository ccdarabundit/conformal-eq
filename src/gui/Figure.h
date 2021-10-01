// FilterFig - Champ Darabundit 09/30/21
// Custom item with inspiration from Jatin's FilterViewer
// https://github.com/Chowdhury-DSP/ChowKick/blob/main/src/gui/FilterViewer.h
#pragma once


#include "../ConformalEQ.h"
#include "FilterPlot.h"

class Figure :  public Component,
                public SettableTooltipClient,
                private Timer
{
public:
    Figure(FilterPlot* fp);
    
    enum ColourIDs
    {
        backgroundColour,
        gridColour
    };
    
    void resized() override;
    void paint (Graphics& g) override;
    
    void timerCallback() override;
private:
    FilterPlot* plotter;
    Colour plotColor[4] = {
        Colours::blue,
        Colours::orange,
        Colours::green,
        Colours::red
    };
    std::vector<double> frequencies;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Figure)
};

class FigureItem : public foleys::GuiItem
{
public:
    FOLEYS_DECLARE_GUI_FACTORY(FigureItem)

    FigureItem(foleys::MagicGUIBuilder& builder, const juce::ValueTree& node)
    : foleys::GuiItem (builder, node)
    {
        auto plugin = dynamic_cast<ConformalEQ*> (builder.getMagicState().getProcessor());
        comp = std::make_unique<Figure> (plugin->getFilterPlot());
        addAndMakeVisible(comp.get());
        setColourTranslation({
            { "background", Figure::backgroundColour },
            { "grid", Figure::gridColour},
            });
    }
    
    void update() override
    {
        comp->resized();
    }
    
    juce::Component* getWrappedComponent() override
    {
        return comp.get();
    }
private:
    std::unique_ptr<Figure> comp;
};
