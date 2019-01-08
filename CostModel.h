#ifndef pixelbridge_CostModel_h
#define pixelbridge_CostModel_h

/**
 * \file CostModel.h
 *
 * \brief This file holds the implementation of the CostModel.
 *
 * This file holds the implementation of the CostModel. The CostModel is used
 * by nDDI implementations to register the cost of the various operations triggered
 * by nDDI commands.
 */

#include <cassert>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

///\cond
// TODO(CDE): Sort out a clean way of including this, which is normally in NDimensionalDisplayInterface.h.
#define COEFFICIENT_UNCHANGED INT16_MAX
///\endcond

/*
 * Definitions for widths for the various
 * data in the NDDI display. Modify to perform
 * various cost experiments.
 */
#ifdef USE_NARROW_DATA_FIELDS
/// \brief Configures the number of bytes for a pixel.
///
/// Configures the number of bytes for a pixel.
#ifdef USE_ALPHA_CHANNEL
#define BYTES_PER_PIXEL     4
#else
#define BYTES_PER_PIXEL     3
#endif

/// \brief Configures the number of bytes for a frame volume coordinate.
///
/// Configures the number of bytes for a frame volume coordinate.
#define BYTES_PER_FV_COORD  4

/// \brief Configures the number of bytes for a coefficient plane coordinate.
///
/// Configures the number of bytes for a coefficient plane coordinate.
#define BYTES_PER_CP_COORD  2

/// \brief Configures the number of bytes for coefficient matrix coordinate.
///
/// Configures the number of bytes for a coefficient matrix coordinate.
#define BYTES_PER_CM_COORD  1

/// \brief Configures the number of bytes for an input vector value.
///
/// Configures the number of bytes for an input vector value.
#define BYTES_PER_IV_VALUE  4

/// \brief Configures the number of bytes for a coefficient.
///
/// Configures the number of bytes for a coefficient.
#define BYTES_PER_COEFF     4

/// \brief Configures the number of bytes for a scaler.
///
/// Configures the number of bytes for a scaler.
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
/// \brief Helper macro that calculates the bytes used to encode c pixels.
///
/// Helper macro that calculates the bytes used to encode c pixels.
#define CALC_BYTES_FOR_PIXELS(c)              (BYTES_PER_PIXEL * c)

/// \brief Helper macro that calculates the bytes used to encode c frame volume coordinate tuples.
///
/// Helper macro that calculates the bytes used to encode c frame volume coordinate tuples.
#define CALC_BYTES_FOR_FV_COORD_TUPLES(c)     (BYTES_PER_FV_COORD * frameVolumeDimensionalSizes_.size() * c)

/// \brief Helper macro that calculates the bytes used to encode c tile coordinate tuples.
///
/// Helper macro that calculates the bytes used to encode c tile coordinate tuples.
#define CALC_BYTES_FOR_TILE_COORD_DOUBLES(c)  (BYTES_PER_CP_COORD * 2 * c)

/// \brief Helper macro that calculates the bytes used to encode c input vector values.
///
/// Helper macro that calculates the bytes used to encode c input vector values.
#define CALC_BYTES_FOR_IV_UPDATE()            (BYTES_PER_IV_VALUE * (inputVector_->getSize() - 2))

/// \brief Helper macro that calculates the bytes used to encode c coefficient matrices.
///
/// Helper macro that calculates the bytes used to encode c coefficent matrices.
#define CALC_BYTES_FOR_CMS(c)                 (BYTES_PER_COEFF * inputVector_->getSize() * frameVolumeDimensionalSizes_.size() * c)

/// \brief Helper macro that calculates the bytes used to encode c coefficient matrix coordinate tuples
///
/// Helper macro that calculates the bytes used to encode c coefficient matrix coordinate tuples
#define CALC_BYTES_FOR_CM_COORD_DOUBLES(c)    (BYTES_PER_CM_COORD * 2 * c)

/// \brief Helper macro that calculates the bytes used to encode c coefficient plane coordinate tuples
///
/// Helper macro that calculates the bytes used to encode c coefficient plane coordinate tuples
#define CALC_BYTES_FOR_CP_COORD_TRIPLES(c)    (BYTES_PER_CP_COORD * 3 * c)


using namespace std;

/**
 * \brief Namespace for the entire nDDI API.
 *
 * Namespace for the entire nDDI API.
 */
namespace nddi {

    /**
     * \brief Enumerates the five different NDDI components for which charges are tracked.
     *
     * Enumerates the five different NDDI components for which charges are tracked: nDDI
     * link, input vector, coefficient matrix, scaler, frame volume. The coefficient matrix
     * and scaler components of the coefficient planes are tracked seperately.
     */
    typedef enum {
        NDDI_LINK_COMPONENT,
        INPUT_VECTOR_COMPONENT,
        COEFFICIENT_MATRIX_COMPONENT,
        SCALER_COMPONENT,
        FRAME_VOLUME_COMPONENT
    } component_t;

    /**
     * \brief Each charge is either a READ_ACCESS or WRITE_ACCESS to the component.
     *
     * Each charge is either a READ_ACCESS or WRITE_ACCESS to the component.
     */
    typedef enum {
        READ_ACCESS,
        WRITE_ACCESS
    } memory_access_t;
    ///\cond
    static const char* memory_access_text[] = {"\"READ_ACCESS\"", "\"WRITE_ACCESS\""};
    ///\endcond

    /**
     * \brief Used to designate which, if any, of the detailed memory component
     *        charges to record.
     *
     * The bits of the logcosts argument to the CostModel constructor indicate
     * which of the memory components to log detailed charges for. The default is
     * NO_CHARGES, which will still log simple count-based charges such as bytes
     * transmitted and commands sent. Setting any of these memory area charges
     * will generate very detailed memory accesses which can be expensive.
     */
    enum log_charges {
        NO_CHARGES = 0,
        IV_CHARGES = 1 << 0,
        CM_CHARGES = 1 << 1,
        SC_CHARGES = 1 << 2,
        FV_CHARGES = 1 << 3,
        ALL_CHARGES = IV_CHARGES | CM_CHARGES | SC_CHARGES | FV_CHARGES
    };

    /**
     * \brief Base class for a memory component charge.
     *
     * Base class for a memory component charge.
     */
    class Charge {
    private:
    public:
        /// \brief Unique sequence number for the charge.
        ///
        /// Unique sequence number for the charge.
        unsigned int sequenceNumber;

	/// \brief Empty base constructor.
	///
	/// Empty base constructor.
        Charge(unsigned int sequenceNumber) : sequenceNumber(sequenceNumber) {}

	/// \brief Empty base printer.
	///
	/// Empty base printer. Each subclass will implement a specific printer that
	/// encode the charge as JSON.
        virtual void print(ofstream &file) {}
    };

    /**
     * \brief Represents the details of data read from or written to the Input Vector.
     *
     * Represents the details of data read from or written to the Input Vector.
     */
    class InputVectorCharge : public Charge {
    public:
        /// \brief The type of memory component access.
        ///
        /// The type of memory component access.
        memory_access_t  access;

	/// \brief The start position of the access.
	///
	/// The start position of the access.
        unsigned int     start;

	/// \brief The end position of the access.
	///
	/// The end position of the access.
        unsigned int     end;
	
	/**
	 * \brief InputVectorCharge constructor.
	 *
	 * InputVectorCharge constructor.
	 *
	 * @param sequenceNumber The unique sequence number to use for this new charge.
	 * @param access The type of memory component access.
	 * @param start The start position of the access.
	 * @param end The end position of the access.
	 */
        InputVectorCharge(unsigned int sequenceNumber, memory_access_t access, unsigned int start, unsigned int end)
        : Charge(sequenceNumber), access(access), start(start), end(end) {}

	/**
	 * \brief Print the JSON for the this memory component charge.
	 *
	 * Print the JSON for the this memory component charge.
	 *
	 * @param file Output file stream for the JSON.
	 */
        void print(ofstream &file) {
            file << "{" << endl;
            file << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            file << "  \"inputVectorCharge\" : {" << endl;
            file << "    \"access\" : " << memory_access_text[access] << "," << endl;
            file << "    \"start\" : " << start << "," << endl;
            file << "    \"end\" : " << end << endl;
            file << "  }" << endl;
            file << "}" << endl;
        }
    };

    /**
     * \brief Represents the details of data read from or written to single Coefficients.
     *
     * Represents the details of data read from or written to single Coefficients.
     */
    class CoefficientCharge : public Charge {
    public:
	/// \brief The type of memory component access.
	///
	/// The type of memory component access.
        memory_access_t       access;

	/// \brief The start coordinate of the access.
	///
	/// The start coordinate of the access.
        vector<unsigned int>  start;

	/// \brief The end coordinate of the access.
	///
	/// The end coordinate of the access.
        vector<unsigned int>  end;

	/// \brief The row of the coefficient being accessed.
	///
	/// The row of the coefficient being accessed.
        int                   row;

	/// \brief The column of the coefficient being accessed.
	///
	/// The column of the coefficient being accessed.
        int                   col;

	/**
	 * \brief CoefficientCharge constructor.
	 *
	 * CoefficientCharge constructor.
	 *
	 * @param sequenceNumber The unique sequence number to use for this new charge.
	 * @param access The type of memory component access.
	 * @param start The start coordinate of the access.
	 * @param end The end coordinate of the access.
	 * @param row The row of the coefficient being accessed.
	 * @param col The column of the coefficient being accessed.
	 */
        CoefficientCharge(unsigned int sequenceNumber, memory_access_t access, vector<unsigned int> start, vector<unsigned int> end, int row, int col)
        : Charge(sequenceNumber), access(access), start(start), end(end), row(row), col(col) {}

	/**
	 * \brief Print the JSON for the this memory component charge.
	 *
	 * Print the JSON for the this memory component charge.
	 *
	 * @param file Output file stream for the JSON.
	 */
        void print(ofstream &file) {
            file << "{" << endl;
            file << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            file << "  \"coefficientCharge\" : {" << endl;
            file << "    \"access\" : " << memory_access_text[access] << "," << endl;
            file << "    \"start\" : " << "[" << start[0];
            for (int i = 1; i < start.size(); i++) { file << "," << start[i]; }
            file << "]," << endl;
            file << "    \"end\" : " << "[" << end[0];
            for (int i = 1; i < end.size(); i++) { file << "," << end[i]; }
            file << "]," << endl;
            file << "    \"row\" : " << row << "," << endl;
            file << "    \"col\" : " << col << endl;
            file << "  }" << endl;
            file << "}" << endl;
        }
    };

    /**
     * \brief Represents the details of data read from or written to Coefficient Matrices.
     *
     * Represents the details of data read from or written to Coefficient Matrices.
     */
    class CoefficientMatrixCharge : public Charge {
    public:
	/// \brief The type of memory component access.
	///
	/// The type of memory component access.
        memory_access_t        access;

	/// \brief The start coordinate of the access.
	///
	/// The start coordinate of the access.
        vector<unsigned int>   start;

	/// \brief The end coordinate of the access.
	///
	/// The end coordinate of the access.
        vector<unsigned int>   end;

	/// \brief Holds a coefficient with values masked out used to indicate which coefficients
	///        within the matrices were updated.
	///
	/// Holds a coefficient matrix with values masked out used to indicate which coefficients
	/// within the matrices were updated.
        vector< vector<int> >  cm;

	/**
	 * \brief CoefficientMatrixCharge constructor.
	 *
	 * CoefficientMatrixCharge constructor.
	 *
	 * @param sequenceNumber The unique sequence number to use for this new charge.
	 * @param access The type of memory component access.
	 * @param start The start coordinate of the access.
	 * @param end The end coordinate of the access.
	 * @param cm Holds a coefficient matrix with values masked out used to indicate which coefficients
	 *           within the matrices were updated.
	 */
        CoefficientMatrixCharge(unsigned int sequenceNumber, memory_access_t access, vector<unsigned int> start, vector<unsigned int> end, vector< vector<int> > cm)
        : Charge(sequenceNumber), access(access), start(start), end(end), cm(cm) {}

	/**
	 * \brief Print the JSON for the this memory component charge.
	 *
	 * Print the JSON for the this memory component charge.
	 *
	 * @param file Output file stream for the JSON.
	 */
        void print(ofstream &file) {
            file << "{" << endl;
            file << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            file << "  \"coefficientMatrixCharge\" : {" << endl;
            file << "    \"access\" : " << memory_access_text[access] << "," << endl;
            file << "    \"start\" : " << "[" << start[0];
            for (int i = 1; i < start.size(); i++) { file << "," << start[i]; }
            file << "]," << endl;
            file << "    \"end\" : " << "[" << end[0];
            for (int i = 1; i < end.size(); i++) { file << "," << end[i]; }
            file << "]," << endl;
            file << "    \"cm\" : " << "[";
            for (int row = 0; row < cm.size(); row++) {
                if (row != 0) { file << ","; }
                file << "[";
                for (int col = 0; col < cm[row].size(); col++) {
                    if (col != 0) { file << ","; }
                    file << cm[row][col];
                }
                file << "]";
            }
            file << "]" << endl;
            file << "  }" << endl;
            file << "}" << endl;
        }
    };

    /**
     * \brief Represents the details of data read from or written to Scalers.
     *
     * Represents the details of data read from or written to Scalers.
     */
    class ScalerCharge : public Charge {
    public:
	/// \brief The type of memory component access.
	///
	/// The type of memory component access.
        memory_access_t       access;

	/// \brief The start coordinate of the access.
	///
	/// The start coordinate of the access.
        vector<unsigned int>  start;

	/// \brief The end coordinate of the access.
	///
	/// The end coordinate of the access.
        vector<unsigned int>  end;

	/**
	 * \brief ScalerCharge constructor.
	 *
	 * ScalerCharge constructor.
	 *
	 * @param sequenceNumber The unique sequence number to use for this new charge.
	 * @param access The type of memory component access.
	 * @param start The start coordinate of the access.
	 * @param end The end coordinate of the access.
	 */
        ScalerCharge(unsigned int sequenceNumber, memory_access_t access, vector<unsigned int> start, vector<unsigned int> end)
        : Charge(sequenceNumber), access(access), start(start), end(end) {}

	/**
	 * \brief Print the JSON for the this memory component charge.
	 *
	 * Print the JSON for the this memory component charge.
	 *
	 * @param file Output file stream for the JSON.
	 */
	void print(ofstream &file) {
            file << "{" << endl;
            file << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            file << "  \"scalerCharge\" : {" << endl;
            file << "    \"access\" : " << memory_access_text[access] << "," << endl;
            file << "    \"start\" : " << "[" << start[0];
            for (int i = 1; i < start.size(); i++) { file << "," << start[i]; }
            file << "]," << endl;
            file << "    \"end\" : " << "[" << end[0];
            for (int i = 1; i < end.size(); i++) { file << "," << end[i]; }
            file << "]" << endl;
            file << "  }" << endl;
            file << "}" << endl;
        }
    };

    /**
     * \brief Represents the details of data read from or written to the Frame Volume.
     *
     * Represents the details of data read from or written to the Frame Volume.
     */
    class FrameVolumeCharge : public Charge {
    public:
	/// \brief The type of memory component access.
	///
	/// The type of memory component access.
        memory_access_t       access;

	/// \brief The start coordinate of the access.
	///
	/// The start coordinate of the access.
        vector<unsigned int>  start;

	/// \brief The end coordinate of the access.
	///
	/// The end coordinate of the access.
        vector<unsigned int>  end;

	/**
	 * \brief FrameVolumeCharge constructor.
	 *
	 * FrameVolumeCharge constructor.
	 *
	 * @param sequenceNumber The unique sequence number to use for this new charge.
	 * @param access The type of memory component access.
	 * @param start The start coordinate of the access.
	 * @param end The end coordinate of the access.
	 */
        FrameVolumeCharge(unsigned int sequenceNumber, memory_access_t access, vector<unsigned int> start, vector<unsigned int> end)
        : Charge(sequenceNumber), access(access), start(start), end(end) {}

	/**
	 * \brief Print the JSON for the this memory component charge.
	 *
	 * Print the JSON for the this memory component charge.
	 *
	 * @param file Output file stream for the JSON.
	 */
        void print(ofstream &file) {
            file << "{" << endl;
            file << "  \"sequenceNumber\" : " << sequenceNumber << "," << endl;
            file << "  \"frameVolumeCharge\" : {" << endl;
            file << "    \"access\" : " << memory_access_text[access] << "," << endl;
            file << "    \"start\" : " << "[" << start[0];
            for (int i = 1; i < start.size(); i++) { file << "," << start[i]; }
            file << "]," << endl;
            file << "    \"end\" : " << "[" << end[0];
            for (int i = 1; i < end.size(); i++) { file << "," << end[i]; }
            file << "]" << endl;
            file << "  }" << endl;
            file << "}" << endl;
        }
    };

    /**
     * \brief The main CostModel class.
     *
     * The CostModel provides a means to log the a variety of NDDI charges including commands sent,
     * bytes transmitted, pixels blended, pixels mapped, and memory component accesses.
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

        vector<unsigned int> fvDimensions;
        unsigned int inputVectorSize = 0;
        vector<Charge*> charges;

        bool headless = false;
        unsigned char logcosts = NO_CHARGES;

        ofstream logfile;

    public:

	/**
	 * \brief Basic constructor that does not include memory component sizes.
	 *
	 * Basic constructor that does include memory component sizes.
	 *
	 * @param headless Indicates if the display is headless, passed here for convenient retrieval by the application later.
	 * @param logcosts 8-bit field holding values from the log_costs enumeration used to indicate which memory component charges should be logged in detail.
	 */
        CostModel(bool headless, unsigned char logcosts)
        : headless(headless), logcosts(logcosts) {
            clearCosts();
        }

	/**
	 * \brief Main constructor which includes memory component sizes.
	 *
	 * Main constructor which includes memory component sizes.
	 *
	 * @param fvDimensions Vector holding the dimensions of the frame volume.
	 * @param inputVectorSize Integer holding the size of the input vector.
	 * @param headless Indicates if the display is headless, passed here for convenient retrieval by the application later.
	 * @param logcosts 8-bit field holding values from the log_costs enumeration used to indicate which memory component charges should be logged in detail.
	 */
        CostModel(vector<unsigned int> &fvDimensions, unsigned int inputVectorSize, bool headless, unsigned char logcosts)
        : fvDimensions(fvDimensions), inputVectorSize(inputVectorSize), headless(headless), logcosts(logcosts) {
            clearCosts();
        }

	/**
	 * \brief Destructor.
	 *
	 * Destructor.
	 */
        ~CostModel() {
            for (int i = 0; i < charges.size(); i++) {
                delete(charges[i]);
            }
        }

	/**
	 * \brief Zeroes the various counts and deletes any detailed memory component charges.
	 *
	 * Zeroes the various counts and deletes any detailed memory component charges. Used
	 * primarily after initial display setup to exclude those costs.
	 */
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
            for (int i = 0; i < charges.size(); i++) {
                delete(charges[i]);
            }
            charges.clear();
        }

	/**
	 * \brief Registers a charge for the input vector memory component.
	 *
	 * Registers a charge for the input vector memory component.
	 * Will update the global counts and optionally log a detailed charge if
	 * the cost model was constructed with the proper log_charges.
	 *
	 * @param access The type of memory component access.
	 * @param start The start position of the access.
	 * @param end The end position of the access.
	 * @param count The number of identical charges to register. Default is 1.
	 */
        void registerInputVectorMemoryCharge(
                memory_access_t access,
                unsigned int start,
                unsigned int end,
                unsigned int count = 1) {

            unsigned int bytes = (end - start + 1) * BYTES_PER_IV_VALUE * count;
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

            if (logcosts & IV_CHARGES) {
                InputVectorCharge* c = new InputVectorCharge(charges.size(), access, start, end);
                charges.push_back(c);
            }
        }

	/**
	 * \brief Registers a single coefficient charge for the coefficient planes memory component.
	 *
	 * Registers a single coefficient charge for the coefficient planes memory component.
	 * Will update the global counts and optionally log a detailed charge if
	 * the cost model was constructed with the proper log_charges.
	 *
	 * @param access The type of memory component access.
	 * @param start The start coordinate of the access.
	 * @param end The end coordinate of the access.
	 * @param cmRow The row of the coefficient being accessed.
	 * @param cmCol The column of the coefficient being accessed.
	 */
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

            if (logcosts & CM_CHARGES) {
                CoefficientCharge* c = new CoefficientCharge(charges.size(), access, start, end, cmRow, cmCol);
                charges.push_back(c);
            }
        }

	/**
	 * \brief Registers a coefficent matrix charge for the coeffcient planes memory component.
	 *
	 * Registers a coefficent matrix charge for the coefficient planes memory component.
	 * Will update the global counts and optionally log a detailed charge if
	 * the cost model was constructed with the proper log_charges.
	 *
	 * @param access The type of memory component access.
	 * @param start The start coordinate of the access.
	 * @param end The end coordinate of the access.
	 * @param coefficientMatrix Holds a coefficient matrix with values masked out used to indicate
	 *                          which coefficients within the matrices were updated.
	 */
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

            if (logcosts & CM_CHARGES) {
                CoefficientMatrixCharge* c = new CoefficientMatrixCharge(charges.size(), access, start, end, coefficientMatrix);
                charges.push_back(c);
            }
        }

	/**
	 * \brief Registers a scaler charge for the coefficient planes memory component.
	 *
	 * Registers a scaler charge for the coefficient planes memory component.
	 * Will update the global counts and optionally log a detailed charge if
	 * the cost model was constructed with the proper log_charges.
	 *
	 * @param access The type of memory component access.
	 * @param start The start coordinate of the access.
	 * @param end The end coordinate of the access.
	 */
        void registerScalerMemoryCharge(
                memory_access_t access,
                vector<unsigned int> &start,
                vector<unsigned int> &end) {

            assert(start.size() == end.size());
            unsigned int bytes = BYTES_PER_SCALER;
            for (int i = 0; i < start.size(); i++) {
                bytes *= end[i] - start[i] + 1;
            }
            // TODO(CDE): Update to use specific scaler counts. PixelBridge statistics and csv will need to be updated.
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

            if (logcosts & SC_CHARGES) {
                ScalerCharge* c = new ScalerCharge(charges.size(), access, start, end);
                charges.push_back(c);
            }
        }

	/**
	 * \brief Registers a charge for the frame volume memory component.
	 *
	 * Registers a charge for the frame volume memory component.
	 * Will update the global counts and optionally log a detailed charge if
	 * the cost model was constructed with the proper log_charges.
	 *
	 * @param access The type of memory component access.
	 * @param start The start coordinate of the access.
	 * @param end The end coordinate of the access.
	 */
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

            if (logcosts & FV_CHARGES) {
                FrameVolumeCharge* c = new FrameVolumeCharge(charges.size(), access, start, end);
                charges.push_back(c);
            }
        }

	/**
	 * \brief Registers a transmission charge.
	 *
	 * Registers a transmission charge, updating both the command count by one and the
	 * number of bytes transmitted.
	 *
	 * @param numBytes The number of bytes to charge.
	 * @param time Deprecated.
	 */
        void registerTransmissionCharge(unsigned long numBytes, unsigned long time) {
#pragma omp atomic
            linkCommandsSent++;

#pragma omp atomic
            linkBytesSent += numBytes;
        }

	/**
	 * \brief Registers a number of pixel blend charges.
	 *
	 * Registers a number of pixel blend charges.
	 *
	 * @param numBlends The number of blends to charge.
	 */
        void registerPixelBlendCharge(unsigned long numBlends) {
#pragma omp atomic
            pixelsBlended += numBlends;
        }

	/**
	 * \brief Registers a number of pixel mapping charges.
	 *
	 * Registers a number of pixel mapping charges. A pixel mapping
	 * is the set of calculations to determine which pixel to map
	 * from the frame volume to a location on the display panel.
	 *
	 * @param numMappings The number of mappings to charge.
	 */
        void registerPixelMappingCharge(unsigned long numMappings) {
#pragma omp atomic
            pixelsMapped += numMappings;
        }

	/**
	 * \brief Used to get the number of commands sent over the nddi link.
	 *
	 * Used to get the number of commands sent over the nddi link.
	 *
	 * @return Count of commands sent over the nddi link.
	 */
        unsigned long getLinkCommandsSent() {
            return linkCommandsSent;
        }

	/**
	 * \brief Used to get the number of bytes sent over the nddi link.
	 *
	 * Used to get the number of bytes sent over the nddi link.
	 *
	 * @return Count of bytes sent over the nddi link.
	 */
        unsigned long getLinkBytesTransmitted() {
            return linkBytesSent;
        }

	/**
	 * \brief Used to get the count of read accesses to a memory component.
	 *
	 * Used to get the count of read accesses to a memory component.
	 *
	 * @param component Specifies which component to get the read access count for.
	 * @return Count of read accesses to a memory component.
	 */
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

 	/**
	 * \brief Used to get the count of write accesses to a memory component.
	 *
	 * Used to get the count of write accesses to a memory component.
	 *
	 * @param component Specifies which component to get the write access count for.
	 * @return Count of write accesses to a memory component.
	 */
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

	/**
	 * \brief Used to get the count of bytes read from a memory component.
	 *
	 * Used to get the count of bytes read from a memory component.
	 *
	 * @param component Specifies which component to get the bytes read for.
	 * @return Count of bytes read from a memory component.
	 */
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

	/**
	 * \brief Used to get the count of bytes written to a memory component.
	 *
	 * Used to get the count of bytes written to a memory component.
	 *
	 * @param component Specifies which component to get the bytes written for.
	 * @return Count of bytes written to a memory component.
	 */
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

	/**
	 * \brief Used the get the count of pixels blended.
	 *
	 * Used the get the count of pixels blended.
	 *
	 * @return Count of pixels blended.
	 */
        unsigned long getPixelsBlended() {
            return pixelsBlended;
        }

	/**
	 * \brief Used the get the count of pixels mapped.
	 *
	 * Used the get the count of pixels mapped.
	 *
	 * @return Count of pixels mapped.
	 */
        unsigned long getPixelsMapped() {
            return pixelsMapped;
        }

	/**
	 * \brief Used to determine if the display owning this CostModel is headless.
	 *
	 * Used to determine if the display owning this CostModel is headless.
	 *
	 * @return True if headless, false otherwise.
	 */
        bool isHeadless() {
            return headless;
        }

	/**
	 * \brief Prints the JSON for all of the charges registered.
	 *
	 * Prints the JSON for all of the charges registered. Prints the configuration
	 * and then any detailed memory component charge that was logged. If an application
	 * wants to print the global counts, then the individual getters can be used to get
	 * the counts and the the application can display them in the desired format.
	 */
        void printCharges() {

            assert(fvDimensions.size() > 0);
            assert(inputVectorSize > 0);

            char filename[24];
            srand(time(NULL) * getpid());
            sprintf(filename, "costlog-%x.json", rand() % 0xffffff + 1);

            logfile.open(filename);
            cout << filename;

            logfile << "{\n" << endl;

            logfile << "\"config\" : {" << endl;
            logfile << "  \"bytePerPixel\": " << BYTES_PER_PIXEL << "," << endl;
            logfile << "  \"bytePerIvValue\": " << BYTES_PER_IV_VALUE << "," << endl;
            logfile << "  \"bytePerCoefficient\": " << BYTES_PER_COEFF << "," << endl;
            logfile << "  \"bytePerScaler\": " << BYTES_PER_SCALER << "," << endl;
            logfile << "  \"inputVectorSize\": " << inputVectorSize << "," << endl;
            logfile << "  \"fvDimensions\": [";
            for (int i = 0; i < fvDimensions.size(); i++) {
                if (i > 0)
                    logfile << ",";
                logfile << fvDimensions[i];
            }
            logfile << "]" << endl;
            logfile << "}," << endl;

            logfile << "\"charges\": [" << endl;
            for (int i = 0; i < charges.size(); i++) {
                if (i > 0)
                    logfile << "," << endl;
                charges[i]->print(logfile);
            }
            logfile << "]" << endl;

            logfile << "\n}" << endl;
            logfile.close();
        }
    };
}
#endif
