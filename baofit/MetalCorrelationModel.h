// Created 24-Apr-2015 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef BAOFIT_METAL_CORRELATION_MODEL
#define BAOFIT_METAL_CORRELATION_MODEL

#include "baofit/AbsCorrelationModel.h"

#include "likely/types.h"

#include <string>
#include <vector>

namespace baofit {
	// Represents an r-space model for metal line correlations
	class MetalCorrelationModel : public AbsCorrelationModel {
	public:
	    // Creates a metal correlation model. Bool parameters metalModel and metalTemplate
	    // are used for selecting the method used to model the metal correlations. The
	    // option metalModel will read a set of required Lya-metal and metal-metal templates
	    // provided by the user. The option metalTemplate uses an empirical template derived
	    // from metals added to mock data.
	    MetalCorrelationModel(std::string const &metalModelName, bool metalModel = false,
	        bool metalModelInterpolate = false, bool metalTemplate = false,
	        bool crossCorrelation = false, AbsCorrelationModel *base = 0);
	    virtual ~MetalCorrelationModel();
	    // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
	protected:
	    // Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h.
	    virtual double _evaluate(double r, double mu, double z, bool anyChanged, int index) const;
	    virtual double _evaluateKSpace(double k, double mu_k, double pk, double z) const;
	    virtual int _getIndexBase() const;
	private:
	    void _initialize(std::vector<double> &vector, std::string const &filename);
	    std::vector<double> _corrLyaSi2a0, _corrLyaSi2a2, _corrLyaSi2a4, _corrLyaSi2b0, _corrLyaSi2b2,
	        _corrLyaSi2b4, _corrLyaSi2c0, _corrLyaSi2c2, _corrLyaSi2c4, _corrLyaSi30, _corrLyaSi32,
	        _corrLyaSi34, _corrSi2aSi2a0, _corrSi2aSi2a2, _corrSi2aSi2a4, _corrSi2aSi2b0, _corrSi2aSi2b2,
	        _corrSi2aSi2b4, _corrSi2aSi2c0, _corrSi2aSi2c2, _corrSi2aSi2c4, _corrSi2bSi2b0, _corrSi2bSi2b2,
	        _corrSi2bSi2b4, _corrSi2bSi2c0, _corrSi2bSi2c2, _corrSi2bSi2c4, _corrSi2cSi2c0, _corrSi2cSi2c2,
	        _corrSi2cSi2c4, _corrSi3Si2a0, _corrSi3Si2a2, _corrSi3Si2a4, _corrSi3Si2b0, _corrSi3Si2b2,
	        _corrSi3Si2b4, _corrSi3Si2c0, _corrSi3Si2c2, _corrSi3Si2c4, _corrSi3Si30, _corrSi3Si32, _corrSi3Si34,
	        _corrLyaCIVa0, _corrLyaCIVa2, _corrLyaCIVa4, _corrSi2aCIVa0, _corrSi2aCIVa2, _corrSi2aCIVa4,
	        _corrSi2bCIVa0, _corrSi2bCIVa2, _corrSi2bCIVa4, _corrSi2cCIVa0, _corrSi2cCIVa2, _corrSi2cCIVa4,
	        _corrSi3CIVa0, _corrSi3CIVa2, _corrSi3CIVa4, _corrCIVaCIVa0, _corrCIVaCIVa2, _corrCIVaCIVa4;
	    int _indexBase, _lastLines;
	    double _rperpMin, _rparMin, _rperpMax, _rparMax;
	    bool _metalModel, _metalModelInterpolate, _metalTemplate, _crossCorrelation;
	    AbsCorrelationModel &_base;
	    likely::BiCubicInterpolatorPtr _LyaSi2a0, _LyaSi2a2, _LyaSi2a4, _LyaSi2b0, _LyaSi2b2, _LyaSi2b4,
	        _LyaSi2c0, _LyaSi2c2, _LyaSi2c4, _LyaSi30, _LyaSi32, _LyaSi34, _Si2aSi2a0, _Si2aSi2a2,
	        _Si2aSi2a4, _Si2aSi2b0, _Si2aSi2b2, _Si2aSi2b4, _Si2aSi2c0, _Si2aSi2c2, _Si2aSi2c4,
	        _Si2bSi2b0, _Si2bSi2b2, _Si2bSi2b4, _Si2bSi2c0, _Si2bSi2c2, _Si2bSi2c4, _Si2cSi2c0,
	        _Si2cSi2c2, _Si2cSi2c4, _Si3Si2a0, _Si3Si2a2, _Si3Si2a4, _Si3Si2b0, _Si3Si2b2, _Si3Si2b4,
	        _Si3Si2c0, _Si3Si2c2, _Si3Si2c4, _Si3Si30, _Si3Si32, _Si3Si34,
	        _QSOSi2a0, _QSOSi2a2, _QSOSi2a4, _QSOSi2b0, _QSOSi2b2, _QSOSi2b4,
	        _QSOSi2c0, _QSOSi2c2, _QSOSi2c4, _QSOSi30, _QSOSi32, _QSOSi34;
    }; // MetalCorrelationModel
    
	// Calculates the normalization factor for each correlation function multipole.
	void updateNormFactors(double &norm0, double &norm2, double &norm4, double biasSq, double betaAvg, double betaProd);

} // baofit

#endif // BAOFIT_METAL_CORRELATION_MODEL