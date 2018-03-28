#ifndef GL_NDDI_DISPLAY_H
#define GL_NDDI_DISPLAY_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#include "BaseNddiDisplay.h"

using namespace nddi;

/**
 * This class adds a public method that returns a reference to the frame buffer
 * so that a GLUT-based application can render it to the window.
 */
class GlNddiDisplay : public nddi::BaseNddiDisplay {

public:
    GlNddiDisplay() {}
    GlNddiDisplay(std::vector<unsigned int> &frameVolumeDimensionalSizes,
                  unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                  bool headless = false,
                  bool fixed8x8Macroblocks = false, bool useSingleCoeffcientPlane = false);
    GlNddiDisplay(std::vector<unsigned int> &frameVolumeDimensionalSizes,
                  unsigned int displayWidth, unsigned int displayHeight,
                  unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                  bool headless = false,
                  bool fixed8x8Macroblocks = false, bool useSingleCoeffcientPlane = false);
    ~GlNddiDisplay();

    /**
     * Used by GLUT-based application to get access to the frame buffer in order to display it.
     *
     * @return Texture holding a rendered frame.
     */
    GLuint GetFrameBufferTex();
    GLuint GetFrameBufferTex(unsigned int sub_x, unsigned int sub_y, unsigned int sub_w, unsigned int sub_h);

    /**
     * Use to register in bulk all of the rendering cost. This will updates the counts correctly,
     * but any more detailed information is lost.
     */
    Pixel* GetFrameBuffer();

    /**
     * Triggers a simulated render that only records the cost estimated cost involved.
     */
    void SimulateRender();

private:
    void Render();
    void Render(unsigned int sub_x, unsigned int sub_y, unsigned int sub_w, unsigned int sub_h);
    void ComputePixels(unsigned int x, unsigned int y, unsigned int length, bool doCostCalculation);
    void RegisterBulkRenderCost();


protected:
    GLuint    texture_;
    Pixel    *frameBuffer_;
};

#endif // GL_NDDI_DISPLAY_H
