// Created 06-Apr-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef BAOFIT_RUNTIME_ERROR
#define BAOFIT_RUNTIME_ERROR

#include <stdexcept>
#include <string>

namespace baofit {
	class RuntimeError : public std::runtime_error {
	public:
		explicit RuntimeError(std::string const &reason);
		virtual ~RuntimeError() throw ();
	private:
	}; // RuntimeError
	
	inline RuntimeError::RuntimeError(std::string const &reason)
	: std::runtime_error(reason) { }

    inline RuntimeError::~RuntimeError() throw () { }
} // baofit

#endif // BAOFIT_RUNTIME_ERROR
