// Created 3-Dec-2015 by Michael Blomqvist (Laboratoire d'Astrophysique de Marseille) <michael.blomqvist@lam.fr>

#ifndef BAOFIT_BAO_DISTORTION_MATRIX
#define BAOFIT_BAO_DISTORTION_MATRIX

#include "likely/types.h"

#include <vector>

namespace baofit {
	// Represents a distortion matrix that models the continuum fitting broadband distortion
	class DistortionMatrix {
	public:
	    // Creates a distortion matrix based on values read from file. Missing entries
	    // are assigned the value zero.
	    DistortionMatrix(std::string const &distMatrixName, int distMatrixOrder, bool verbose = false);
	    virtual ~DistortionMatrix();
	    // Returns the value of the distortion matrix for the specified indices.
	    double getDistortion(int index1, int index2) const;
	    // Sets the value of the undistorted correlation function for the specified bin.
	    void setCorrelation(int bin, double value);
	    // Returns the value of the undistorted correlation function for the specified bin.
	    double getCorrelation(int bin) const;
	protected:
	private:
	    int _nbins;
	    std::vector<double> _ucf, _dist;
    }; // DistortionMatrix
} // baofit

#endif // BAOFIT_BAO_DISTORTION_MATRIX
