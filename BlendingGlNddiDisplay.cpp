//
//  BlendingGlNddiDisplay.cpp
//  pixelbridge
//
//  Created by Dave Estes on 1/31/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>

#include "Features.h"
#include "BlendingGlNddiDisplay.h"

#define CLAMPTOBYTE(color)                  \
if((color) & (~255)) {                      \
color = (unsigned char)((-(color)) >> 31);  \
} else {                                    \
color = (unsigned char)(color);             \
}

using namespace nddi;

// public

BlendingGlNddiDisplay::BlendingGlNddiDisplay(unsigned int frameVolumeDimensionality,
                                             unsigned int frameVolumeDimensionalSizes[],
                                             unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                                             bool headless)
: GlNddiDisplay(frameVolumeDimensionality, frameVolumeDimensionalSizes, numCoefficientPlanes, inputVectorSize)
{
}

BlendingGlNddiDisplay::BlendingGlNddiDisplay(unsigned int frameVolumeDimensionality,
                                             unsigned int frameVolumeDimensionalSizes[],
                                             unsigned int displayWidth, unsigned int displayHeight,
                                             unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                                             bool headless)
: GlNddiDisplay(frameVolumeDimensionality, frameVolumeDimensionalSizes, displayWidth, displayHeight, numCoefficientPlanes, inputVectorSize)
{
}

// Implements the frame volume function directly using getPixel and setPixel since frame volume does not yet
// support blending.
void BlendingGlNddiDisplay::CopyFrameVolume(unsigned int* start, unsigned int* end, unsigned int* dest, bool blend) {

    unsigned int positionFrom[frameVolumeDimensionality_];
    unsigned int positionTo[frameVolumeDimensionality_];
    bool copyFinished = false;
    int pixelsCopied = 0;

    memcpy(positionFrom, start, sizeof(unsigned int) * frameVolumeDimensionality_);
    memcpy(positionTo, dest, sizeof(unsigned int) * frameVolumeDimensionality_);

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_FV_COORD_TUPLES(3), // Three Coordinate Tuples
                                          0);

    // Move from start to end, filling in each location with the provided pixel
    do {
        // Set pixel in frame volume at position
        if (!blend) {
            frameVolume_->setPixel(positionTo, frameVolume_->getPixel(positionFrom));
        } else {
            // Register the cost of the blending
            costModel->registerPixelBlendCharge(1);

            frameVolume_->setPixel(positionTo, BlendPixel(frameVolume_->getPixel(positionTo), frameVolume_->getPixel(positionFrom)));
        }
        pixelsCopied++;

        // Move to the next position
        int fvDim = 0;
        bool overflow;
        do {
            overflow = false;
            positionFrom[fvDim]++;
            positionTo[fvDim]++;
            if ( (positionFrom[fvDim] >= frameVolumeDimensionalSizes_[fvDim])
                || (positionFrom[fvDim] > end[fvDim]) ) {
                overflow = true;
                positionFrom[fvDim] = start[fvDim];
                positionTo[fvDim] = dest[fvDim];
                if (++fvDim >= frameVolumeDimensionality_)
                    copyFinished = true;
            }
        } while (overflow && !copyFinished);

    } while (!copyFinished);
#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

// private

// Don't register the blending cost here, since we want to bypass it sometimes for speed.
nddi::Pixel BlendingGlNddiDisplay::BlendPixel(nddi::Pixel pTo, nddi::Pixel pFrom) {
    int r, g, b;

    r = ((int)pFrom.r * (int)pFrom.a + (int)pTo.r * (255 - (int)pFrom.a)) / 255;
    g = ((int)pFrom.g * (int)pFrom.a + (int)pTo.g * (255 - (int)pFrom.a)) / 255;
    b = ((int)pFrom.b * (int)pFrom.a + (int)pTo.b * (255 - (int)pFrom.a)) / 255;

    CLAMPTOBYTE(r);
    CLAMPTOBYTE(g);
    CLAMPTOBYTE(b);

    pTo.r = (unsigned char)r;
    pTo.g = (unsigned char)g;
    pTo.b = (unsigned char)b;

    return pTo;
}
