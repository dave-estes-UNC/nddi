//
//  CostModel.h
//  pixelbridge
//
//  Created by Dave Estes on 10/18/11.
//  Copyright 2011 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_CostModel_h
#define pixelbridge_CostModel_h

#include <cassert>
#include <iostream>

// TODO(CDE): Sort out a clean way of including this, which is normally in NDimensionalDisplayInterface.h.
#define COEFFICIENT_UNCHANGED INT16_MAX

using namespace std;

/*
 * Definitions for widths for the various
 * data in the NDDI display. Modify to perform
 * various cost experiments.
 */
#ifdef USE_NARROW_DATA_FIELDS
#ifdef USE_ALPHA_CHANNEL
#define BYTES_PER_PIXEL     4
#else
#define BYTES_PER_PIXEL     3
#endif
#define BYTES_PER_FV_COORD  4
#define BYTES_PER_CP_COORD  2
#define BYTES_PER_CM_COORD  1
#define BYTES_PER_IV_VALUE  4
#define BYTES_PER_COEFF     4
#define BYTES_PER_SCALER    6
#else
#define BYTES_PER_PIXEL     4
#define BYTES_PER_FV_COORD  4
#define BYTES_PER_CP_COORD  4
#define BYTES_PER_CM_COORD  4
#define BYTES_PER_IV_VALUE  4
#define BYTES_PER_COEFF     4
#define BYTES_PER_SCALER    8
#endif

/*
 * Helper macros to be used when registering
 * cost charges.
 */
#define CALC_BYTES_FOR_PIXELS(c)              (BYTES_PER_PIXEL * c)
#define CALC_BYTES_FOR_FV_COORD_TUPLES(c)     (BYTES_PER_FV_COORD * frameVolumeDimensionalSizes_.size() * c)
#define CALC_BYTES_FOR_TILE_COORD_DOUBLES(c)  (BYTES_PER_CP_COORD * 2 * c)
#define CALC_BYTES_FOR_IV_UPDATE()            (BYTES_PER_IV_VALUE * (inputVector_->getSize() - 2))
#define CALC_BYTES_FOR_CMS(c)                 (BYTES_PER_COEFF * inputVector_->getSize() * frameVolumeDimensionalSizes_.size() * c)
#define CALC_BYTES_FOR_CM_COORD_DOUBLES(c)    (BYTES_PER_CM_COORD * 2 * c)
#define CALC_BYTES_FOR_CP_COORD_TRIPLES(c)    (BYTES_PER_CP_COORD * 3 * c)

using namespace std;

namespace nddi {

    typedef enum {
        NDDI_LINK_COMPONENT,
        INPUT_VECTOR_COMPONENT,
        COEFFICIENT_MATRIX_COMPONENT,
        SCALER_COMPONENT,
        FRAME_VOLUME_COMPONENT
    } component_t;

    typedef enum {
        READ_ACCESS,
        WRITE_ACCESS
    } memory_access_t;
    static const char* memory_access_text[] = {"READ_ACCESS", "WRITE_ACCESS"};

    class Charge {
    public:
        unsigned int sequenceNumber;
        Charge(unsigned int sequenceNumber) : sequenceNumber(sequenceNumber) {}
        virtual void print() {}
    };

    /**
     * Represents the details of data read from or written to the Input Vector
     */
    class InputVectorCharge : public Charge {
    public:
        memory_access_t  access;
        unsigned int     start;
        unsigned int     end;
        InputVectorCharge(unsigned int sequenceNumber, memory_access_t access, unsigned int start, unsigned int end)
        : Charge(sequenceNumber), access(access), start(start), end(end) {}
        void print() {
            cout << "{" << endl;
            cout << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            cout << "  \"inputVectorCharge\" : {" << endl;
            cout << "    \"access\" : " << memory_access_text[access] << "," << endl;
            cout << "    \"start\" : " << start << "," << endl;
            cout << "    \"end\" : " << end << "," << endl;
            cout << "  }" << endl;
            cout << "}," << endl;
        }
    };

    /**
     * Represents the details of data read from or written to single Coefficients
     */
    class CoefficientCharge : public Charge {
    public:
        memory_access_t       access;
        vector<unsigned int>  start;
        vector<unsigned int>  end;
        int                   row;
        int                   col;
        CoefficientCharge(unsigned int sequenceNumber, memory_access_t access, vector<unsigned int> start, vector<unsigned int> end, int row, int col)
        : Charge(sequenceNumber), access(access), start(start), end(end), row(row), col(col) {}
        void print() {
            cout << "{" << endl;
            cout << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            cout << "  \"coefficientCharge\" : {" << endl;
            cout << "    \"access\" : " << memory_access_text[access] << "," << endl;
            cout << "    \"start\" : " << "[" << start[0];
            for (int i = 1; i < start.size(); i++) { cout << "," << start[i]; }
            cout << "]," << endl;
            cout << "    \"end\" : " << "[" << end[0];
            for (int i = 1; i < end.size(); i++) { cout << "," << end[i]; }
            cout << "]," << endl;
            cout << "    \"row\" : " << row << "," << endl;
            cout << "    \"col\" : " << col << "," << endl;
            cout << "  }" << endl;
            cout << "}," << endl;
        }
    };

    /**
     * Represents the details of data read from or written to Coefficient Matrices
     */
    class CoefficientMatrixCharge : public Charge {
    public:
        memory_access_t        access;
        vector<unsigned int>   start;
        vector<unsigned int>   end;
        vector< vector<int> >  cm;
        CoefficientMatrixCharge(unsigned int sequenceNumber, memory_access_t access, vector<unsigned int> start, vector<unsigned int> end, vector< vector<int> > cm)
        : Charge(sequenceNumber), access(access), start(start), end(end), cm(cm) {}
        void print() {
            cout << "{" << endl;
            cout << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            cout << "  \"coefficientMatrixCharge\" : {" << endl;
            cout << "    \"access\" : " << memory_access_text[access] << "," << endl;
            cout << "    \"start\" : " << "[" << start[0];
            for (int i = 1; i < start.size(); i++) { cout << "," << start[i]; }
            cout << "]," << endl;
            cout << "    \"end\" : " << "[" << end[0];
            for (int i = 1; i < end.size(); i++) { cout << "," << end[i]; }
            cout << "]," << endl;
            cout << "    \"cm\" : " << "[";
            for (int row = 0; row < cm.size(); row++) {
                if (row != 0) { cout << ","; }
                cout << "[";
                for (int col = 0; col < cm[row].size(); col++) {
                    if (col != 0) { cout << ","; }
                    cout << cm[row][col];
                }
                cout << "]";
            }
            cout << "]," << endl;
            cout << "  }" << endl;
            cout << "}," << endl;
        }
    };

    /**
     * Represents the details of data read from or written to Scalers.
     */
    class ScalerCharge : public Charge {
    public:
        memory_access_t       access;
        vector<unsigned int>  start;
        vector<unsigned int>  end;
        ScalerCharge(unsigned int sequenceNumber, memory_access_t access, vector<unsigned int> start, vector<unsigned int> end)
        : Charge(sequenceNumber), access(access), start(start), end(end) {}
        void print() {
            cout << "{" << endl;
            cout << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            cout << "  \"scalerCharge\" : {" << endl;
            cout << "    \"access\" : " << memory_access_text[access] << "," << endl;
            cout << "    \"start\" : " << "[" << start[0];
            for (int i = 1; i < start.size(); i++) { cout << "," << start[i]; }
            cout << "]," << endl;
            cout << "    \"end\" : " << "[" << end[0];
            for (int i = 1; i < end.size(); i++) { cout << "," << end[i]; }
            cout << "]," << endl;
            cout << "  }" << endl;
            cout << "}," << endl;
        }
    };

    /**
     * Represents the details of data read from or written to the Frame Volume.
     */
    class FrameVolumeCharge : public Charge {
    public:
        memory_access_t       access;
        vector<unsigned int>  start;
        vector<unsigned int>  end;
        FrameVolumeCharge(unsigned int sequenceNumber, memory_access_t access, vector<unsigned int> start, vector<unsigned int> end)
        : Charge(sequenceNumber), access(access), start(start), end(end) {}
        void print() {
            cout << "{" << endl;
            cout << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            cout << "  \"frameVolumeCharge\" : {" << endl;
            cout << "    \"access\" : " << memory_access_text[access] << "," << endl;
            cout << "    \"start\" : " << "[" << start[0];
            for (int i = 1; i < start.size(); i++) { cout << "," << start[i]; }
            cout << "]," << endl;
            cout << "    \"end\" : " << "[" << end[0];
            for (int i = 1; i < end.size(); i++) { cout << "," << end[i]; }
            cout << "]," << endl;
            cout << "  }" << endl;
            cout << "}," << endl;
        }
    };

    /**
     * Represents the number of bytes sent over the NDDI Link.
     */
    // TODO(CDE): Convert to a proper charge object.
    typedef struct {
        unsigned long    numBytes;
    } link_charge_t;

    /**
     * Represents the number of pixel blend operations where each operation alpha
     * blends two pixels.
     */
    // TODO(CDE): Convert to a proper charge object.
    typedef struct {
        unsigned long     numBlends;
    } pixel_blend_charge_t;

    /**
     * Represents the number of pixel mapping operations, where each operation is
     * a matrix multiplication.
     */
    // TODO(CDE): Convert to a proper charge object.
    typedef struct {
        unsigned long     numMappings;
    } pixel_mapping_charge_t;


    /**
     * The CostModel allows different types of charges to be made and will run reports
     * later.
     */
    class CostModel {

    private:

        unsigned long linkCommandsSent;
        unsigned long linkBytesSent;

        unsigned long pixelsBlended;

        unsigned long pixelsMapped;

        unsigned long inputVectorReads;
        unsigned long inputVectorWrites;
        unsigned long inputVectorBytesRead;
        unsigned long inputVectorBytesWritten;

        unsigned long coefficientPlaneReads;
        unsigned long coefficientPlaneWrites;
        unsigned long coefficientPlaneBytesRead;
        unsigned long coefficientPlaneBytesWritten;

        unsigned long frameVolumeReads;
        unsigned long frameVolumeWrites;
        unsigned long frameVolumeBytesRead;
        unsigned long frameVolumeBytesWritten;

        vector<Charge*> charges;

        bool headless = false;
        bool logcosts = false;

    public:

        CostModel(bool headless, bool logcosts) : headless(headless), logcosts(logcosts) {
            clearCosts();
        }

        ~CostModel() {
            for (int i = 0; i < charges.size(); i++) {
                delete(charges[i]);
            }
        }

        void clearCosts() {
            linkCommandsSent = 0;
            linkBytesSent = 0;
            pixelsBlended = 0;
            pixelsMapped = 0;
            inputVectorReads = 0;
            inputVectorWrites = 0;
            inputVectorBytesRead = 0;
            inputVectorBytesWritten = 0;
            coefficientPlaneReads = 0;
            coefficientPlaneWrites = 0;
            coefficientPlaneBytesRead = 0;
            coefficientPlaneBytesWritten = 0;
            frameVolumeReads = 0;
            frameVolumeWrites = 0;
            frameVolumeBytesRead = 0;
            frameVolumeBytesWritten = 0;
        }

        void registerInputVectorMemoryCharge(
                memory_access_t access,
                unsigned int start,
                unsigned int end) {

            unsigned int bytes = (end - start + 1) * BYTES_PER_IV_VALUE;
            if (access == READ_ACCESS) {
#pragma omp atomic
                inputVectorReads++;
#pragma omp atomic
                inputVectorBytesRead += bytes;
            } else {
#pragma omp atomic
                inputVectorWrites++;
#pragma omp atomic
                inputVectorBytesWritten += bytes;
            }

            if (logcosts) {
                InputVectorCharge* c = new InputVectorCharge(charges.size(), access, start, end);
                charges.push_back(c);
            }
        }

        void registerCoefficientMemoryCharge(
                memory_access_t access,
                vector<unsigned int> &start,
                vector<unsigned int> &end,
                int cmRow, int cmCol) {

            assert(start.size() == end.size());
            unsigned int bytes = BYTES_PER_COEFF;
            for (int i = 0; i < start.size(); i++) {
                bytes *= end[i] - start[i] + 1;
            }
            if (access == READ_ACCESS) {
#pragma omp atomic
                coefficientPlaneReads++;
#pragma omp atomic
                coefficientPlaneBytesRead += bytes;
            } else {
#pragma omp atomic
                coefficientPlaneWrites++;
#pragma omp atomic
                coefficientPlaneBytesWritten += bytes;
            }

            if (logcosts) {
                CoefficientCharge* c = new CoefficientCharge(charges.size(), access, start, end, cmRow, cmCol);
                charges.push_back(c);
            }
        }

        void registerCoefficientMatrixMemoryCharge(
                memory_access_t access,
                vector<unsigned int> &start,
                vector<unsigned int> &end,
                vector< vector<int> > &coefficientMatrix) {

            assert(start.size() == end.size());
            unsigned int bytes = 0;
            for (int row = 0; row < coefficientMatrix.size(); row++) {
                for (int col = 0; col < coefficientMatrix[row].size(); col++) {
                    if (coefficientMatrix[row][col] != COEFFICIENT_UNCHANGED) {
                        bytes += BYTES_PER_COEFF;
                    }
                }
            }
            for (int i = 0; i < start.size(); i++) {
                bytes *= end[i] - start[i] + 1;
            }
            if (access == READ_ACCESS) {
#pragma omp atomic
                coefficientPlaneReads++;
#pragma omp atomic
                coefficientPlaneBytesRead += bytes;
            } else {
#pragma omp atomic
                coefficientPlaneWrites++;
#pragma omp atomic
                coefficientPlaneBytesWritten += bytes;
            }

            if (logcosts) {
                CoefficientMatrixCharge* c = new CoefficientMatrixCharge(charges.size(), access, start, end, coefficientMatrix);
                charges.push_back(c);
            }
        }

        void registerScalerMemoryCharge(
                memory_access_t access,
                vector<unsigned int> &start,
                vector<unsigned int> &end) {

            assert(start.size() == end.size());
            unsigned int bytes = BYTES_PER_SCALER;
            for (int i = 0; i < start.size(); i++) {
                bytes *= end[i] - start[i] + 1;
            }
            // TODO(CDE): Update to use specific scaler counts. PixelBridge statistics and csv will need to be udpated.
            if (access == READ_ACCESS) {
#pragma omp atomic
                coefficientPlaneReads++;
#pragma omp atomic
                coefficientPlaneBytesRead += bytes;
            } else {
#pragma omp atomic
                coefficientPlaneWrites++;
#pragma omp atomic
                coefficientPlaneBytesWritten += bytes;
            }

            if (logcosts) {
                ScalerCharge* c = new ScalerCharge(charges.size(), access, start, end);
                charges.push_back(c);
            }
        }

        void registerFrameVolumeMemoryCharge(
                memory_access_t access,
                vector<unsigned int> &start,
                vector<unsigned int> &end) {

            assert(start.size() == end.size());
            unsigned int bytes = BYTES_PER_PIXEL;
            for (int i = 0; i < start.size(); i++) {
                bytes *= end[i] - start[i] + 1;
            }
            if (access == READ_ACCESS) {
#pragma omp atomic
                frameVolumeReads++;
#pragma omp atomic
                frameVolumeBytesRead += bytes;
            } else {
#pragma omp atomic
                frameVolumeWrites++;
#pragma omp atomic
                frameVolumeBytesWritten += bytes;
            }

            if (logcosts) {
                FrameVolumeCharge* c = new FrameVolumeCharge(charges.size(), access, start, end);
                charges.push_back(c);
            }
        }

        void registerBulkMemoryCharge(component_t component,
                                      unsigned long accessCount,
                                      memory_access_t access,
                                      unsigned long numBytes) {

            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    if (access == READ_ACCESS) {
#pragma omp atomic
                        inputVectorReads += accessCount;
#pragma omp atomic
                        inputVectorBytesRead += numBytes;
                    } else {
#pragma omp atomic
                        inputVectorWrites += accessCount;
#pragma omp atomic
                        inputVectorBytesWritten += numBytes;
                    }
                    break;
                case COEFFICIENT_MATRIX_COMPONENT:
                // TODO(CDE): Add seperate counts for SCALER_COMPONENT. Will require changes anywhere bulk memory charges are made.
                case SCALER_COMPONENT:
                    if (access == READ_ACCESS) {
#pragma omp atomic
                        coefficientPlaneReads += accessCount;
#pragma omp atomic
                        coefficientPlaneBytesRead += numBytes;
                    } else {
#pragma omp atomic
                        coefficientPlaneWrites += accessCount;
#pragma omp atomic
                        coefficientPlaneBytesWritten += numBytes;
                    }
                    break;
                case FRAME_VOLUME_COMPONENT:
                    if (access == READ_ACCESS) {
#pragma omp atomic
                        frameVolumeReads += accessCount;
#pragma omp atomic
                        frameVolumeBytesRead += numBytes;
                    } else {
#pragma omp atomic
                        frameVolumeWrites += accessCount;
#pragma omp atomic
                        frameVolumeBytesWritten += numBytes;
                    }
                    break;
                default:
                    break;
            }
        }

        void registerTransmissionCharge(unsigned long numBytes, unsigned long time) {
#pragma omp atomic
            linkCommandsSent++;

#pragma omp atomic
            linkBytesSent += numBytes;
        }

        void registerPixelBlendCharge(unsigned long numBlends) {
#pragma omp atomic
            pixelsBlended += numBlends;
        }

        void registerPixelMappingCharge(unsigned long numMappings) {
#pragma omp atomic
            pixelsMapped += numMappings;
        }

        unsigned long getLinkCommandsSent() {
            return linkCommandsSent;
        }

        unsigned long getLinkBytesTransmitted() {
            return linkBytesSent;
        }

        unsigned long getReadAccessCount(component_t component) {

            unsigned long count = 0;

            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    count = inputVectorReads;
                    break;
                case COEFFICIENT_MATRIX_COMPONENT:
                    count = coefficientPlaneReads;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    count = frameVolumeReads;
                    break;
                default:
                    break;
            }
            return count;
        }

        unsigned long getWriteAccessCount(component_t component) {

            unsigned long count = 0;

            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    count = inputVectorWrites;
                    break;
                case COEFFICIENT_MATRIX_COMPONENT:
                    count = coefficientPlaneWrites;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    count = frameVolumeWrites;
                    break;
                default:
                    break;
            }
            return count;
        }

        unsigned long getBytesRead(component_t component) {

            unsigned long count = 0;

            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    count = inputVectorBytesRead;
                    break;
                case COEFFICIENT_MATRIX_COMPONENT:
                    count = coefficientPlaneBytesRead;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    count = frameVolumeBytesRead;
                    break;
                default:
                    break;
            }
            return count;
        }

        unsigned long getBytesWritten(component_t component) {

            unsigned long count = 0;

            switch (component) {
                case INPUT_VECTOR_COMPONENT:
                    count = inputVectorBytesWritten;
                    break;
                case COEFFICIENT_MATRIX_COMPONENT:
                    count = coefficientPlaneBytesWritten;
                    break;
                case FRAME_VOLUME_COMPONENT:
                    count = frameVolumeBytesWritten;
                    break;
                default:
                    break;
            }
            return count;
        }

        unsigned long getPixelsBlended() {
            return pixelsBlended;
        }

        unsigned long getPixelsMapped() {
            return pixelsMapped;
        }

        bool isHeadless() {
            return headless;
        }

        void printCharges() {
            for (int i = 0; i < charges.size(); i++) {
                charges[i]->print();
            }
        }
    };
}
#endif
