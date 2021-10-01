// FilterFig - Champ Darabundit 09/30/21
// Custom item with inspiration from Jatin's FilterViewer
// https://github.com/Chowdhury-DSP/ChowKick/blob/main/src/gui/FilterViewer.h

#include "Figure.h"

namespace
{
constexpr float lowFreq = 20.0f;
constexpr float highFreq = 48000.0f;
constexpr float dbRange = 72.0f;
constexpr double fs = 48000.0;

float getFreqForInd (int i, int N)
{
    return lowFreq * std::pow ((highFreq / lowFreq), (float) i/N);
}

float getXForFreq (float freq, int N, float width, float centerX)
{
    auto ind = N*std::log(freq/lowFreq)/std::log(highFreq/lowFreq);
    auto xScale = width/N;
    auto xPos = xScale * ind;
    return xPos - 0.125*(width/2 - centerX);
}

float getYForMagnitude (float mag, float height, float center)
{
    auto magDB = jmin (Decibels::gainToDecibels (mag), 70.0f);
    auto yScale = height/dbRange;
    auto yPos =  center - yScale * magDB;
    return yPos;
}
}

Figure::Figure(FilterPlot* fp)
{
    setColour(backgroundColour, Colours::black);
    setColour(gridColour, Colours::dimgrey);
    plotter = fp;
    
    // Define our frequencies
    auto f = fp->get_frequencies();
    frequencies.resize(f.size());
    frequencies = f;
}

void Figure::resized()
{
    repaint();
}

void Figure::paint(Graphics& g)
{
    // Paint graph
    g.fillAll(findColour (backgroundColour));
    const auto width = (float) getWidth();
    const auto height = (float) getHeight();
    const auto left = (float) getX();
    const auto top = (float) getY();
    const auto centerY = (float) top + 0.5*height;
    const auto centerX = (float) left + 0.5*width;
    const int numPoints = frequencies.size();
    
    // Draw X grid lines
    int pow10 = 0;
    float lineFreq = (float) std::pow(10, pow10);
    while ( lineFreq < lowFreq)
    {
        pow10++;
        lineFreq = (float) std::pow(10, pow10);
    }
    while (lineFreq < highFreq)
    {
        auto xPos = getXForFreq(lineFreq, numPoints, width, centerX);
        Path gridline;
        gridline.startNewSubPath( xPos, top);
        gridline.lineTo( xPos, top + height);
        g.setColour( findColour (gridColour));
        juce::PathStrokeType stroke(0.5f);
        g.strokePath(gridline, stroke);
        if (lineFreq == 100.0)
            g.drawSingleLineText("100", xPos + 2, top+10);
        if (lineFreq == 1000.0)
            g.drawSingleLineText("1k", xPos + 2, top+10);
        if (lineFreq == 10000.0)
            g.drawSingleLineText("10k", xPos + 2, top+10);
        for (int i = 2; i < 9; i++)
        {
            Path minorgrid;
            xPos = getXForFreq(i*lineFreq, numPoints, width, centerX);
            minorgrid.startNewSubPath( xPos, top);
            minorgrid.lineTo( xPos, top + height);
            g.setColour( findColour (gridColour));
            const float dash[2] = {1, 4};
            juce::PathStrokeType stroke(0.5f);
            stroke.createDashedStroke(minorgrid, minorgrid, dash, 2);
            g.strokePath(minorgrid, stroke);
        }
        pow10++;
        lineFreq = (float) std::pow(10, pow10);
    }
    
    for (int i = 1; i < 4; i++)
    {
        Path gridline;
        double dBVal = -0.5*dbRange + i*dbRange/4;
        float magVal = Decibels::decibelsToGain(dBVal);
        float yPos = getYForMagnitude(magVal, height, centerY);
        gridline.startNewSubPath( left, top+yPos);
        gridline.lineTo( left+width, top+yPos);
        g.setColour( findColour (gridColour));
        juce::PathStrokeType stroke(0.5f);
        g.strokePath(gridline, stroke);
        g.drawSingleLineText(std::to_string((int) dBVal), left + 2, yPos+top-2);
    }
    
    // Use the plotter to iterate and plot all graphs
    for (int i = 0; i < plotter->getN(); i++)
    {
        Path line;
        auto magnitude = plotter->get_mag(i);
        float yPos = getYForMagnitude((float) magnitude[0], height, centerY);
        float xPos = getXForFreq(frequencies[0], numPoints, width, centerX);
        line.startNewSubPath(xPos, top + yPos);
        for (int j = 1; j < static_cast<int>(magnitude.size()) ; j++)
        {
            yPos = getYForMagnitude((float) magnitude[j], height, centerY);
            xPos = getXForFreq(frequencies[j], numPoints, width, centerX);
            line.lineTo( xPos , top + yPos );
        }
        g.setColour( plotColor[i & 3]);
        g.strokePath(line, PathStrokeType(2.0f, PathStrokeType::JointStyle::curved));
    }

}

void Figure::timerCallback()
{
    repaint();
}
