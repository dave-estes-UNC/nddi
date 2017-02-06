//
//  BlendingGlNddiDisplay.h
//  pixelbridge
//
//  Created by Dave Estes on 1/31/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_BlendingGlNddiDisplay_h
#define pixelbridge_BlendingGlNddiDisplay_h

#include "GlNddiDisplay.h"
#include "NDimensionalDisplayInterfaceExtended.h"

using namespace nddi;
using namespace std;

/**
 * Blending version of a GL NDDI Display.
 */
class BlendingGlNddiDisplay : public GlNddiDisplay, public NDimensionalDisplayInterfaceExtended {

public:
    BlendingGlNddiDisplay(unsigned int frameVolumeDimensionality,
                          unsigned int* frameVolumeDimensionalSizes,
                          unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                          bool headless = false);
    BlendingGlNddiDisplay(unsigned int frameVolumeDimensionality,
                          unsigned int* frameVolumeDimensionalSizes,
                          unsigned int displayWidth, unsigned int displayHeight,
                          unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                          bool headless = false);
    ~BlendingGlNddiDisplay();

    // To satisfy the NDimensionalDisplayInterfaceExtended interface
    void CopyFrameVolume(unsigned int* start, unsigned int* end, unsigned int* dest, bool blend);

private:
    nddi::Pixel BlendPixel(nddi::Pixel pTo, nddi::Pixel pFrom);

};

#endif
