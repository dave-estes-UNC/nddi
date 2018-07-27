#ifndef GL_NDDI_DISPLAY_H
#define GL_NDDI_DISPLAY_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#include "SimpleNddiDisplay.h"

using namespace nddi;

/**
 * This class adds a public method that returns a reference to the frame buffer
 * so that a GLUT-based application can render it to the window.
 */
class GlNddiDisplay : public SimpleNddiDisplay {

public:
    GlNddiDisplay() {}
    GlNddiDisplay(std::vector<unsigned int> &frameVolumeDimensionalSizes,
                  unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                  bool headless = false, bool logcosts = false,
                  bool fixed8x8Macroblocks = false, bool useSingleCoeffcientPlane = false);
    GlNddiDisplay(std::vector<unsigned int> &frameVolumeDimensionalSizes,
                  unsigned int displayWidth, unsigned int displayHeight,
                  unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                  bool headless = false, bool logcosts = false,
                  bool fixed8x8Macroblocks = false, bool useSingleCoeffcientPlane = false);
    ~GlNddiDisplay();

    /**
     * Used by GLUT-based application to get access to the frame buffer in order to display it.
     *
     * @return Texture holding a rendered frame.
     */
    GLuint GetFrameBufferTex();
    GLuint GetFrameBufferTex(unsigned int sub_x, unsigned int sub_y, unsigned int sub_w, unsigned int sub_h);

protected:
    GLuint    texture_;
};

#endif // GL_NDDI_DISPLAY_H
