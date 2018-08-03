#include <iostream>
#include <stdio.h>
#include <sys/time.h>

#include "Features.h"
#include "GlNddiDisplay.h"

// public

GlNddiDisplay::GlNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                             unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                             bool headless, unsigned char logcosts, bool fixed8x8Macroblocks, bool useSingleCoeffcientPlane) {
    texture_ = 0;
    GlNddiDisplay(frameVolumeDimensionalSizes, 320, 240, numCoefficientPlanes, inputVectorSize);
}

GlNddiDisplay::GlNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                             unsigned int displayWidth, unsigned int displayHeight,
                             unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                             bool headless, unsigned char logcosts, bool fixed8x8Macroblocks, bool useSingleCoeffcientPlane)
: SimpleNddiDisplay(frameVolumeDimensionalSizes, displayWidth, displayHeight, numCoefficientPlanes, inputVectorSize, headless, logcosts, fixed8x8Macroblocks, useSingleCoeffcientPlane) {

    // allocate a texture name
    glGenTextures( 1, &texture_ );
}

// TODO(CDE): Why is the destructor for GlNddiDisplay being called when we're using a ClNddiDisplay?
GlNddiDisplay::~GlNddiDisplay() {
    glDeleteTextures(1, &texture_);
}

// Private

GLuint GlNddiDisplay::GetFrameBufferTex() {
    return GetFrameBufferTex(0, 0, displayWidth_, displayHeight_);
}

GLuint GlNddiDisplay::GetFrameBufferTex(unsigned int sub_x, unsigned int sub_y, unsigned int sub_w, unsigned int sub_h) {

#ifdef SUPRESS_EXCESS_RENDERING
    if (changed_)
        Render(sub_x, sub_y, sub_w, sub_h);
#endif

// TODO(CDE): Temporarily putting this here until GlNddiDisplay and ClNddiDisplay
//            are using the exact same kind of GL textures
#ifndef USE_CL
    // select our current texture
    glBindTexture( GL_TEXTURE_2D, texture_ );

    // select modulate to mix texture with color for shading
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );

    // when texture area is small, bilinear filter the closest mipmap
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR_MIPMAP_NEAREST );
    // when texture area is large, bilinear filter the first mipmap
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    // if wrap is true, the texture wraps over at the edges (repeat)
    //       ... false, the texture ends at the edges (clamp)
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
                    GL_CLAMP );

    // build our texture mipmaps
    gluBuild2DMipmaps( GL_TEXTURE_2D, 3, displayWidth_, displayHeight_,
                      GL_RGBA, GL_UNSIGNED_BYTE, frameBuffer_ );
#endif

    return texture_;
}
