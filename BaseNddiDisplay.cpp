#include <iostream>
#include <cassert>
#include <cmath>

#include "Features.h"
#include "BaseNddiDisplay.h"

using namespace nddi;

// public

BaseNddiDisplay::BaseNddiDisplay() :
        displayWidth_(0),
        displayHeight_(0),
        inputVector_(NULL),
        frameVolume_(NULL),
        coefficientPlanes_(NULL),
        costModel(NULL),
        quiet_(false),
        changed_(false)
{}

BaseNddiDisplay::BaseNddiDisplay(unsigned int frameVolumeDimensionality,
                                 unsigned int* frameVolumeDimensionalSizes,
                                 unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                                 bool headless)
: displayWidth_(0),
  displayHeight_(0),
  numPlanes_(numCoefficientPlanes),
  inputVector_(NULL),
  frameVolumeDimensionality_(frameVolumeDimensionality),
  frameVolume_(NULL),
  coefficientPlanes_(NULL),
  costModel(NULL),
  quiet_(false),
  changed_(false) {
    frameVolumeDimensionalSizes_ = (unsigned int*)malloc(sizeof(unsigned int) * frameVolumeDimensionality);
    memcpy(frameVolumeDimensionalSizes_, frameVolumeDimensionalSizes, sizeof(unsigned int) * frameVolumeDimensionality);
}

BaseNddiDisplay::BaseNddiDisplay(unsigned int frameVolumeDimensionality,
                                 unsigned int* frameVolumeDimensionalSizes,
                                 unsigned int displayWidth, unsigned int displayHeight,
                                 unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                                 bool headless)
: displayWidth_(displayWidth),
  displayHeight_(displayHeight),
  numPlanes_(numCoefficientPlanes),
  inputVector_(NULL),
  frameVolumeDimensionality_(frameVolumeDimensionality),
  frameVolume_(NULL),
  coefficientPlanes_(NULL),
  costModel(NULL),
  quiet_(false),
  changed_(false) {
    frameVolumeDimensionalSizes_ = (unsigned int*)malloc(sizeof(unsigned int) * frameVolumeDimensionality);
    memcpy(frameVolumeDimensionalSizes_, frameVolumeDimensionalSizes, sizeof(unsigned int) * frameVolumeDimensionality);
}

BaseNddiDisplay::~BaseNddiDisplay() {
    free((void*)frameVolumeDimensionalSizes_);
}

unsigned int BaseNddiDisplay::DisplayWidth() {
    return displayWidth_;
}

unsigned int BaseNddiDisplay::DisplayHeight() {
    return displayHeight_;
}

unsigned int BaseNddiDisplay::NumCoefficientPlanes() {
    return numPlanes_;
}

void BaseNddiDisplay::PutPixel(Pixel p, unsigned int* location) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(1) +         // One Pixel
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(1), // One Coordinate Tuple
                                          0);

    // Set the single pixel
    frameVolume_->PutPixel(p, location);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::CopyPixelStrip(Pixel* p, unsigned int* start, unsigned int* end) {

    int dimensionToCopyAlong;
    bool dimensionFound = false;

    // Find the dimension to copy along
    for (int i = 0; !dimensionFound && (i < frameVolumeDimensionality_); i++) {
        if (start[i] != end[i]) {
            dimensionToCopyAlong = i;
            dimensionFound = true;
        }
    }
    int pixelsToCopy = end[dimensionToCopyAlong] - start[dimensionToCopyAlong] + 1;

    // Register transmission cost now that we know the length of the strip sent
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(pixelsToCopy) +    // A strip of pixels
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(2),       // Two Coordinate Tuples
                                          0);

    // Copy the pixels
    frameVolume_->CopyPixelStrip(p, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::CopyPixels(Pixel* p, unsigned int* start, unsigned int* end) {

    // Register transmission cost first
    int pixelsToCopy = 1;
    for (int i = 0; i < frameVolumeDimensionality_; i++) {
        pixelsToCopy *= end[i] - start[i] + 1;
    }
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(pixelsToCopy) +    // Range of pixels
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(2),       // Two Coordinate Tuples
                                          0);

    // Copy pixels
    frameVolume_->CopyPixels(p, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::CopyPixelTiles(Pixel** p, unsigned int* starts, unsigned int* size, size_t count) {

    // Register transmission cost first
    int pixelsToCopy = 1;
    for (int i = 0; i < 2; i++) {
        pixelsToCopy *= size[i];
    }
    pixelsToCopy *= count;
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(pixelsToCopy) +       // t tiles of x by y tiles of pixels
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(count + 1) + // t start coordinate tuples + 1 tuple for tile size dimensions
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1),       // 1 X by Y tile dimension double
                                          0);

    // Copy pixels
    unsigned int start[frameVolumeDimensionality_];
    unsigned int end[frameVolumeDimensionality_];
    for (int i = 0; i < count; i++) {
        start[0] = starts[i * frameVolumeDimensionality_ + 0];
        start[1] = starts[i * frameVolumeDimensionality_ + 1];
        end[0] = start[0] + size[0]- 1; if (end[0] >= frameVolumeDimensionalSizes_[0]) end[0] = frameVolumeDimensionalSizes_[0] - 1;
        end[1] = start[1] + size[1] - 1; if (end[1] >= frameVolumeDimensionalSizes_[1]) end[1] = frameVolumeDimensionalSizes_[1] - 1;
        for (int j = 2; j < frameVolumeDimensionality_; j++) {
            start[j] = end[j] = starts[i * frameVolumeDimensionality_ + j];
        }
        frameVolume_->CopyPixels(p[i], start, end);
    }

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillPixel(Pixel p, unsigned int* start, unsigned int* end) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(1) +         // One Pixel
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(2), // Two Coordinate Tuples
                                          0);

    // Fill pixels
    frameVolume_->FillPixel(p, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::CopyFrameVolume(unsigned int* start, unsigned int* end, unsigned int* dest) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_FV_COORD_TUPLES(3), // Three Coordinate Tuples
                                          0);

    // Copy pixels
    frameVolume_->CopyFrameVolume(start, end, dest);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::UpdateInputVector(int* input) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_IV_UPDATE(), // Input Vector
                                          0);

    // Update the input vector
    inputVector_->UpdateInputVector(input);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::PutCoefficientMatrix(int* coefficientMatrix, unsigned int* location) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_CMS(1) +             // One coefficient matrix
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(1), // One Coefficient Plane Coordinate triple
                                          0);

    // Update the coefficient matrix
    coefficientPlanes_->PutCoefficientMatrix(coefficientMatrix, location);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillCoefficientMatrix(int* coefficientMatrix, unsigned int* start, unsigned int* end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_CMS(1) +             // One coefficient matrix
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2), // Two Coefficient Plane Coordinate triples
                                          0);

    // Fill the coefficient matrices
    coefficientPlanes_->FillCoefficientMatrix(coefficientMatrix, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillCoefficient(int coefficient, unsigned int row, unsigned int col, unsigned int* start, unsigned int* end) {
    assert(row >= 0 && row < CM_HEIGHT);
    assert(col >= 0 && col < CM_WIDTH);

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_COEFF * 1 +                // One coefficient
                                          CALC_BYTES_FOR_CM_COORD_DOUBLES(1) + // One Coefficient Matrix Coordinate double
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2),  // Two Coefficient Plane Coordinate triples
                                          0);

    // Fill the coefficient matrices
    coefficientPlanes_->FillCoefficient(coefficient, row, col, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillCoefficientTiles(int* coefficients, unsigned int* positions, unsigned int* starts, unsigned int* size, size_t count) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_COEFF * count +                 // count coefficients
                                          CALC_BYTES_FOR_CM_COORD_DOUBLES(count) +  // count Coefficient Matrix Coordinate doubles
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(count) +  // count Coefficient Plane Coordinate triples
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1),     // 1 X by Y tile dimension double
                                          0);

    // Fill the coefficient matrices
    unsigned int* start = starts;
    unsigned int end[3];
    unsigned int* position = positions;
    for (size_t i = 0; i < count; i++) {
        end[0] = start[0] + size[0] - 1; if (end[0] >= displayWidth_) end[0] = displayWidth_ - 1;
        end[1] = start[1] + size[1] - 1; if (end[1] >= displayHeight_) end[1] = displayHeight_ - 1;
        end[2] = start[2];
        coefficientPlanes_->FillCoefficient(coefficients[i], position[0], position[1], start, end);
        start += 3;
        position += 2;
    }

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillScaler(Scaler scaler, unsigned int* start, unsigned int* end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_SCALER * 1 +              // One Scaler
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2), // Two Coefficient Plane Coordinate triples
                                          0);

    // Fill the coefficient matrices
    coefficientPlanes_->FillScaler(scaler, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillScalerTiles(Scaler* scalers, unsigned int* starts, unsigned int* size, size_t count) {
    Scaler s;

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_SCALER * count +                // count scalers
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(count) +  // count Coefficient Plane Coordinate triples
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1),     // One X by Y tile dimension double
                                          0);

    unsigned int* start = starts;
    unsigned int end[3];
    for (size_t i = 0; i < count; i++) {
        end[0] = start[0] + size[0] - 1; if (end[0] >= displayWidth_) end[0] = displayWidth_ - 1;
        end[1] = start[1] + size[1] - 1; if (end[1] >= displayHeight_) end[1] = displayHeight_ - 1;
        end[2] = start[2];
        s.packed = scalers[i].packed;
        coefficientPlanes_->FillScaler(s, start, end);
        start += 3;
    }

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::FillScalerTileStack(Scaler* scalers, unsigned int* start, unsigned int* size, size_t count) {
    Scaler s;

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_SCALER * count +            // count scalers
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(1) +  // One Coefficient Plane Coordinate triples
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1), // One X by Y tile dimension double
                                          0);

    unsigned int st[] = {start[0], start[1], start[2]};
    unsigned int end[3];
    for (size_t i = 0; i < count; i++) {
        end[0] = start[0] + size[0] - 1; if (end[0] >= displayWidth_) end[0] = displayWidth_ - 1;
        end[1] = start[1] + size[1] - 1; if (end[1] >= displayHeight_) end[1] = displayHeight_ - 1;
        end[2] = start[2];
        s.packed = scalers[i].packed;
        coefficientPlanes_->FillScaler(s, st, end);
        st[2]++;
    }

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void BaseNddiDisplay::SetPixelByteSignMode(SignMode mode) {
    assert(mode == UNSIGNED_MODE || mode == SIGNED_MODE);
    pixelSignMode_ = mode;
}

void BaseNddiDisplay::SetFullScaler(uint16_t scaler) {
    if (scaler & (scaler - 1)) {
        cout << "ERROR: THE FULL_SCALER specified is not a power of two." << endl;
    } else {
        // Set the full scaler
        fullScaler_ = scaler;

        // Initialize the shifter used during accumulation
        double s = log2((double)fullScaler_);
        accumulatorShifter_ = int(s);
    }
}

CostModel* BaseNddiDisplay::GetCostModel() {
    return costModel;
}
