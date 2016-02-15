// Created 07-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_ABS_CORRELATION_MODEL
#define BAOFIT_ABS_CORRELATION_MODEL

#include "likely/FitModel.h"

#include "cosmo/types.h"

#include <string>
#include <vector>

namespace baofit {
	class AbsCorrelationModel : public likely::FitModel {
	// Represents an abstract parameterized model of a two-point correlation function.
	public:
	    // Creates a new model with the specified name.
		AbsCorrelationModel(std::string const &name);
		virtual ~AbsCorrelationModel();
		// Returns the correlation function evaluated in redshift space where (r,mu) is
		// the pair separation and z is their average redshift. The separation r should
		// be provided in Mpc/h. Updates our current parameter values. If the parameter
		// delta-v is non-zero, the separation triangle defined by (r,mu) will be modified
		// by displacing the edge along the line of sight by an amount dpi = (dv/100)(1+z)/H(z)
		// in Mpc/h, where dv = delta-v is in km/s.
        double evaluate(double r, double mu, double z, likely::Parameters const &params, int index);
        // Returns the correlation function for the specified multipole at co-moving pair separation
        // r and average pair redshift z. Updates our current parameter values. The value of
        // delta-v has no effect.
        double evaluate(double r, cosmo::Multipole multipole, double z, likely::Parameters const &params,
            int index);
        // Sets the grid coordinates to use for the distortion matrix.
        void setCoordinates(std::vector<double> rbin, std::vector<double> mubin,
            std::vector<double> zbin);
        // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
    protected:
        friend class BaoCorrelationModel;
        friend class BaoKSpaceCorrelationModel;
        friend class BaoKSpaceFftCorrelationModel;
        friend class BaoKSpaceHybridCorrelationModel;
        friend class BroadbandModel;
        friend class MetalCorrelationModel;
        friend class NonLinearCorrectionModel;
        // The public methods above call these protected methods after making parameter values
        // and changes available via our base class' getParameterValue() and isParameterValueChanged()
        // methods. Any registered changes to parameter values are reset after calling any of these.
        virtual double _evaluate(double r, double mu, double z, bool changed, int index) const = 0;
        // We provide a default implementation of the (r,ell,z) evaluator that performs the
        // projection integral over mu weighted with LegendreP(ell) numerically.
        virtual double _evaluate(double r, cosmo::Multipole multipole, double z, bool changed, int index) const;
        // k-space
        virtual double _evaluateKSpace(double k, double mu_k, double pk, double z) const = 0;
        // Defines the standard set of linear bias parameters used by _getNormFactor below, in
        // addition to a parameter "delta-v" used by _applyVelocityShift below. Adds a second set
        // of bias and beta parameters if crossCorrelation is true. Returns the index
        // of the last parameter defined.
        int _defineLinearBiasParameters(double zref, bool crossCorrelation = false);
        // Returns the reference redshift
        double _getZRef() const;
        // Sets the reference redshift
        void _setZRef(double zref);
        // Sets the index of the Lya beta parameter to use.
        void _setBetaIndex(int betaIndex);
        // Sets the index of the Lya bias*(1+beta) parameter to use.
        void _setBbIndex(int bbIndex);
        // Sets the index of the bias power-law evolution parameter to use.
        void _setGammaBiasIndex(int gammabiasIndex);
        // Sets the index of the beta power-law evolution parameter to use.
        void _setGammaBetaIndex(int gammabetaIndex);
        // Sets the index of the Lya beta*bias parameter to use.
        void _setBetaBiasIndex(int betabiasIndex);
        // Sets the index of the velocity shift parameter to use.
        void _setDVIndex(int index);
        // Sets the index of the bias2 parameter to use.
        void _setBias2Index(int bias2Index);
        // Sets the index of the beta2*bias2 parameter to use.
        void _setBeta2Bias2Index(int beta2bias2Index);
        // Updates the internal parameters.
        void _updateInternalParameters();
        // Returns the Lya beta parameter.
        double _getBeta() const;
        // Returns the Lya bias parameter.
        double _getBias() const;
        // Returns the bias power-law evolution parameter.
        double _getGammaBias() const;
        // Returns the beta power-law evolution parameter.
        double _getGammaBeta() const;
        // Returns the beta2 parameter.
        double _getBeta2() const;
        // Returns the bias2 parameter.
        double _getBias2() const;
        // Applies a shift dpi in the parallel direction to the separation (r,mu) using a fiducial
        // cosmology to convert from dv in km/s to dpi in Mpc/h at the specified redshift z. The
        // value of dv is obtained from the "delta-v" parameter defined by _defineLinearBiasParameters.
        void _applyVelocityShift(double &r, double &mu, double z);
        // Updates the multipole normalization factors b^2(z)*C_ell(beta(z)) returned by getNormFactor(ell).
        double _getNormFactor(cosmo::Multipole multipole, double z) const;
        // Returns the radius in Mpc/h for the specified bin.
        double _getRBin(int index) const;
        // Returns the cosine of the theta angle for the specified bin.
        double _getMuBin(int index) const;
        // Returns the redshift for the specified bin.
        double _getZBin(int index) const;
    private:
        int _indexBase, _dvIndex, _betaIndex, _bbIndex, _gammabiasIndex, _gammabetaIndex, _betabiasIndex,
            _bias2Index, _beta2bias2Index, _nbins;
        bool _crossCorrelation;
        enum IndexOffset {
            BETA = 0, BB = 1, GAMMA_BIAS = 2, GAMMA_BETA = 3, DELTA_V = 4, BIAS2 = 5, BB2 = 6
        };
        double _zref, _beta, _bias, _gammaBias, _gammaBeta, _bias2, _beta2;
        std::vector<double> _rbin, _mubin, _zbin;
        mutable std::vector<double> _rbinShift, _mubinShift, _zbinShift;
	}; // AbsCorrelationModel

    inline double AbsCorrelationModel::_getZRef() const { return _zref; }
    inline void AbsCorrelationModel::_setDVIndex(int index) { _dvIndex = index; }
    inline void AbsCorrelationModel::_setBetaIndex(int betaIndex) { _betaIndex = betaIndex; }
    inline void AbsCorrelationModel::_setBbIndex(int bbIndex) { _bbIndex = bbIndex; }
    inline void AbsCorrelationModel::_setGammaBiasIndex(int gammabiasIndex) { _gammabiasIndex = gammabiasIndex; }
    inline void AbsCorrelationModel::_setGammaBetaIndex(int gammabetaIndex) { _gammabetaIndex = gammabetaIndex; }
    inline void AbsCorrelationModel::_setBetaBiasIndex(int betabiasIndex) { _betabiasIndex = betabiasIndex; }
    inline void AbsCorrelationModel::_setBias2Index(int bias2Index) { _bias2Index = bias2Index; }
    inline void AbsCorrelationModel::_setBeta2Bias2Index(int beta2bias2Index) { _beta2bias2Index = beta2bias2Index; }
    inline double AbsCorrelationModel::_getBeta() const { return _beta; }
    inline double AbsCorrelationModel::_getBias() const { return _bias; }
    inline double AbsCorrelationModel::_getGammaBias() const { return _gammaBias; }
    inline double AbsCorrelationModel::_getGammaBeta() const { return _gammaBeta; }
    inline double AbsCorrelationModel::_getBeta2() const { return _beta2; }
    inline double AbsCorrelationModel::_getBias2() const { return _bias2; }
    // Evaluates the redshift evolution p(z) of a parameter for which p(zref)=p0 according to
    // p(z) = p0*((1+z)/(1+zref))^gamma.
    double redshiftEvolution(double p0, double gamma, double z, double zref);

} // baofit

#endif // BAOFIT_ABS_CORRELATION_MODEL
