// Created 31-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_CORRELATION_ANALYZER
#define BAOFIT_CORRELATION_ANALYZER

#include "baofit/types.h"

#include "cosmo/types.h"

#include "likely/BinnedDataResampler.h"
#include "likely/FitParameter.h"

#include <iosfwd>

namespace baofit {
    // Accumulates correlation data and manages its analysis.
	class CorrelationAnalyzer {
	public:
	    // Creates a new analyzer using the specified minimization method.
	    // The range [rmin,rmax] will be used for dumping any model multipoles.
		CorrelationAnalyzer(std::string const &method, double rmin, double rmax, int covSampleSize = 0,
		    bool verbose = true, bool scalarWeights = false);
		virtual ~CorrelationAnalyzer();
		// Set the verbose level during analysis.
        void setVerbose(bool value);
		// Adds a new correlation data object to this analyzer. Reuse the covariance of a
		// previously added dataset specified by reuseCovIndex, unless it is < 0. Returns
		// the index of the newly added dataset.
        int addData(AbsCorrelationDataCPtr data, int reuseCovIndex);
        // Returns the number of data objects added to this analyzer.
        int getNData() const;
        // Sets the correlation model to use.
        void setModel(AbsCorrelationModelPtr model);
        // Returns a shared pointer to the combined correlation data added to this
        // analyzer, after it has been finalized. If verbose, prints out the number
        // of bins with data before and after finalizing the data.
        AbsCorrelationDataPtr getCombined(bool verbose = false, bool finalized = true) const;
        // Fits the combined correlation data aadded to this analyzer and returns
        // the estimated function minimum. Use the optional config script to modify
        // the initial parameter configuration used for the fit (any changes do not
        // propagate back to the model or modify subsequent fits).
        likely::FunctionMinimumPtr fitSample(AbsCorrelationDataCPtr sample,
            std::string const &config = "") const;        
        // Performs a bootstrap analysis and returns the number of fits to bootstrap
        // samples that failed. Specify a non-zero bootstrapSize to generate trials with
        // a number of observations different than getNData(). Specify a refitConfig script
        // to fit each bootstrap sample twice: first with the default model config, then
        // with the refit config script applied. In this case, a trial is only considered
        // successful if both fits succeed. If a saveFile is specified, the parameter
        // values and chi-square value from each fit will be saved to the specified filename.
        // If nsave > 0, then the best-fit model multipoles will be appended to each line
        // using the oneLine option to dumpModel(). In case of refits, the output from both
        // fits will be concatenated on each line. Setting fixCovariance to false means that
        // fits will use a covariance matrix that does not correctly account for double
        // counting. See likely::BinnedDataResampler::bootstrap for details.
        int doBootstrapAnalysis(int bootstrapTrials, int bootstrapSize, bool fixCovariance,
            likely::FunctionMinimumPtr fmin,
            likely::FunctionMinimumPtr fmin2 = likely::FunctionMinimumPtr(),
            std::string const &refitConfig = "", std::string const &saveName = "", int nsave = 0,
            double zsave = -1) const;
        // Performs a jackknife analysis and returns the number of fits that failed. Set
        // jackknifeDrop to the number of observations to drop from each sample. See
        // doBootstrapAnalysis for a description of the other parameters.
        int doJackknifeAnalysis(int jackknifeDrop, likely::FunctionMinimumPtr fmin, 
            likely::FunctionMinimumPtr fmin2 = likely::FunctionMinimumPtr(),
            std::string const &refitConfig = "", std::string const &saveName = "", int nsave = 0,
            double zsave = -1) const;
        // Performs a Markov-chain sampling of the likelihood function for the combined data with
        // the current model, using the specified function minimum to initialize the sampling.
        // Saves nchain samples, using only one per interval trials. See doBootstrapAnalysis for a
        // description of the other parameters.
        void generateMarkovChain(int nchain, int interval, likely::FunctionMinimumCPtr fmin,
            std::string const &saveName = "", int nsave = 0, double zsave = -1) const;
        // Compares each observation to the combined observations, saving one line per observation
        // to the specified filename with the format (k indexes observations):
        //
        //   k log(|C_k|)/n chi2_k chi2_k,1 chi2_k,2 ... chi2_k,nbins
        //
        // where chi2_k = (d_k - D).(C_k - C)^(-1).(d_k - D) for combined data D and covariance C,
        // which should be chi-square distributed with nbins dof if each C_k is correct. The
        // chi2_k,i are the contributions to chi2_k associated with each eigenmode of C_k-C, and
        // sum up to give chi2_k. The distribution of chi2_k,i should be chi-square with 1 dof.
        // The individual observations and combination will be finalized if requested.
        void compareEach(std::string const &saveName, bool finalized) const;
        // Fits each observation separately and returns the number of fits that failed.
        // See doBootstrapAnalysis for a description of the other parameters.
        int fitEach(likely::FunctionMinimumPtr fmin,
            likely::FunctionMinimumPtr fmin2 = likely::FunctionMinimumPtr(),
            std::string const &refitConfig = "", std::string const &saveName = "", int nsave = 0,
            double zsave = -1) const;
        // Refits the specified sample on the parameter grid specified by each parameter's binning
        // spec and saves the results to the specified file name. Returns the number of fits performed.
        int parameterScan(likely::FunctionMinimumCPtr fmin,
            AbsCorrelationDataCPtr sample, std::string const &saveName = "", int nsave = 0,
            double zsave = -1) const;
        // Generates and fits toy Monte Carlo samples and returns the number of fits that failed.
        // Samples are generated by calculating the truth corresponding to the best-fit input
        // parameters in fmin and adding noise sampled from the combined dataset covariance matrix.
        // Each sample is fit using the combined dataset covariance matrix. Saves the first generated
        // sample to the specified filename, if one is provided. The covariance matrix used for noise
        // sampling is scaled by the specified factor. See doBootstrapAnalysis for a description of
        // the other parameters.
        int doToyMCSampling(int ngen, std::string const &mcConfig, std::string const &mcSaveFile,
            double varianceScale, likely::FunctionMinimumPtr fmin, likely::FunctionMinimumPtr fmin2,
            std::string const &refitConfig, std::string const &saveName, int nsave, double zsave) const;
        // Dumps the data, prediction, and diagonal error for each bin of the specified combined
        // data set to the specified output stream. The fit result is assumed to correspond
        // to model that is currently associated with this analyzer. Use the optional script
        // to modify the parameters used in the model. By default, the gradient of each
        // bin with respect to each floating parameter is append to each output row, unless
        // dumpGradients = false.
        void dumpResiduals(std::ostream &out, likely::FunctionMinimumPtr fmin,
            AbsCorrelationDataCPtr combined, std::string const &script = "",
            bool dumpGradients = true) const;
        // Dumps the model predictions for the specified fit parameters to the specified
        // output stream. The input parameters are assumed to correspond to the model that is
        // currently associated with this analyzer. Use the optional script to modify
        // the parameters that will be used to evaluate the model. By default, values are output
        // as "rval mono quad hexa" on separate lines. With oneLine = true, values of
        // "mono quad hexa" are concatenated onto a single line. Uses the specified zdump
        // redshift to evaluate the model. If zdump < 0, then uses the redshift of the first
        // bin in the data grid.
        void dumpModel(std::ostream &out, likely::FitParameters parameters,
            int ndump, double zdump, std::string const &script = "", bool oneLine = false) const;
        // Fills the vector provided with the decorrelated weights of the specified data using
        // the specified parameter values.
        void getDecorrelatedWeights(AbsCorrelationDataCPtr data, likely::Parameters const &params,
            std::vector<double> &dweights) const;
        // Calculates and prints the redshift where the error on the parameter named scaleName
        // has a minimum, assuming that it evolves according to a parameter named "gamma-alpha".
        // Returns true if successful, or false unless both scaleName and "gamma-alpha" are
        // floating parameters of the specified function minimum. Uses the specified zref to
        // calculate the redshift evolution of the scale and its error.
        bool printScaleZEff(likely::FunctionMinimumCPtr fmin, double zref, std::string const &scaleName) const;
        // Returns a bootstrap estimate of the combined data's covariance matrix (before any final cuts)
        // using the specified number of bootstrap trials. The individual observations are combined
        // with inverse covariance weights for each resampling, but then each resampling is fed to
        // an unweighted covariance accumulator (so the total covariance of each resampling is not used).
        // See likely::BinnedDataResampler for more details.
        likely::CovarianceMatrixPtr
            estimateCombinedCovariance(int nSamples, std::string const &filename) const;
	private:
        std::string _method;
        double _rmin, _rmax;
        int _covSampleSize;
        bool _verbose;
        likely::BinnedDataResampler _resampler;
        AbsCorrelationModelPtr _model;
        
        class AbsSampler;
        class JackknifeSampler;
        class BootstrapSampler;
        class EachSampler;
        class ToyMCSampler;
        int doSamplingAnalysis(AbsSampler &sampler, std::string const &method,
            likely::FunctionMinimumPtr fmin, likely::FunctionMinimumPtr fmin2,
            std::string const &refitConfig, std::string const &saveName, int nsave, double zsave) const;
        
	}; // CorrelationAnalyzer
	
    inline void CorrelationAnalyzer::setVerbose(bool value) { _verbose = value; }
    inline int CorrelationAnalyzer::getNData() const { return _resampler.getNObservations(); }
    inline void CorrelationAnalyzer::setModel(AbsCorrelationModelPtr model) { _model = model; }

} // baofit

#endif // BAOFIT_CORRELATION_ANALYZER
