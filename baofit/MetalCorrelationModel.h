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
	    // provided by the user. The option toyMetal uses an empirical toy model derived
	    // from metals added to mock data.
	    MetalCorrelationModel(std::string const &metalModelName, bool metalModel = false,
	        bool metalModelInterpolate = false, bool metalCIV = false, bool toyMetal = false,
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
	    void _initializeGrid(std::vector<double> &vector, std::string const &filename);
	    std::vector<std::vector<double> > _metaltemplates, _zgrid;
	    std::vector<likely::BiCubicInterpolatorPtr> _metalintertemplates;
	    std::vector<int> _paramindex1, _paramindex2;
	    int _indexBase, _lastLines, _nmet, _ncomb;
	    double _rperpMin, _rparMin, _rperpMax, _rparMax;
	    bool _metalModel, _metalModelInterpolate, _metalCIV, _toyMetal, _crossCorrelation;
	    AbsCorrelationModel &_base;
    }; // MetalCorrelationModel
    
	// Calculates the normalization factor for each correlation function multipole.
	void updateNormFactors(double &norm0, double &norm2, double &norm4, double biasSq, double betaAvg, double betaProd);

} // baofit

#endif // BAOFIT_METAL_CORRELATION_MODEL