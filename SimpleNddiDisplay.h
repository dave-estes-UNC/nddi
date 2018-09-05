#ifndef SIMPLE_NDDI_DISPLAY_H
#define SIMPLE_NDDI_DISPLAY_H

#include "BaseNddiDisplay.h"

using namespace nddi;

/**
 * This class adds a public method that returns a reference to the frame buffer
 * so that a GLUT-based application can render it to the window.
 */
class SimpleNddiDisplay : public BaseNddiDisplay {

public:
    SimpleNddiDisplay() {}
    SimpleNddiDisplay(std::vector<unsigned int> &frameVolumeDimensionalSizes,
                  unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                  bool headless = false, unsigned char logcosts = NO_CHARGES,
                  bool fixed8x8Macroblocks = false, bool useSingleCoeffcientPlane = false);
    SimpleNddiDisplay(std::vector<unsigned int> &frameVolumeDimensionalSizes,
                  unsigned int displayWidth, unsigned int displayHeight,
                  unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                  bool headless = false, unsigned char logcosts = NO_CHARGES,
                  bool fixed8x8Macroblocks = false, bool useSingleCoeffcientPlane = false);
    ~SimpleNddiDisplay();

    /**
     * Use to register in bulk all of the rendering cost. This will updates the counts correctly,
     * but any more detailed information is lost.
     */
    Pixel* GetFrameBuffer();
    Pixel* GetFrameBuffer(unsigned int sub_x, unsigned int sub_y, unsigned int sub_w, unsigned int sub_h);

protected:
    void Render();
    void Render(unsigned int sub_x, unsigned int sub_y, unsigned int sub_w, unsigned int sub_h);
    void ComputePixels(unsigned int x, unsigned int y, unsigned int length, bool doCostCalculation);


protected:
    Pixel    *frameBuffer_;
    std::vector<std::vector<int> > emptyCoefficientMatrix;
};

#endif // SIMPLE_NDDI_DISPLAY_H
