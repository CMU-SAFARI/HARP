/**
 * @file gf2_noz3.h
 *
 * @brief Representation of GF(2) data types as compatible with Eigen
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#ifndef GF2_NOZ3_H
#define GF2_NOZ3_H

#include "Eigen/Eigen"

/**
 * @brief representation of a GF(2) element
 * 
 * Note: DOES NOT SUPPORT MULTIPLE Z3 CONTEXTS!! IT USES A SINGLE STATIC GLOBAL ONE.
 */
class gf2_noz3
{	
private:
	bool e;

public:
	gf2_noz3() : e(0) {}
	gf2_noz3(const gf2_noz3 &v) : e(v.get_value()) {}
	gf2_noz3(bool n) : e(n) {}
	gf2_noz3(int n) : e(n == 1) { assert((n & 1) == n && "creating GF(2) value from non-boolean"); }
	// gf2_noz3(double n) : e(n == 1.0) { assert((n == 0.0 || n == 1.0) && "creating GF(2) value from non-boolean"); }
	~gf2_noz3() {}
	
    bool get_value(void) const { return e; }

	gf2_noz3 &operator=(bool n) { e = n; return *this; }
	gf2_noz3 &operator=(int n)  { assert((n & 1) == n && "creating GF(2) value from non-boolean"); e = n; return *this; }

	gf2_noz3 operator==(const gf2_noz3 &rhs) const { return this->e == rhs.get_value(); }
	gf2_noz3 operator==(const bool rhs) const { return this->e == rhs; }
	
    gf2_noz3 operator!=(const gf2_noz3 &rhs) const { return this->e != rhs.get_value(); }
	gf2_noz3 operator!=(const bool rhs) const { return this->e != rhs; }
	
	gf2_noz3 operator&&(const gf2_noz3 &rhs) const { return this->e && rhs.get_value(); }
	gf2_noz3 operator&&(const bool rhs) const { return this->e && rhs; }
	
	gf2_noz3 operator||(const gf2_noz3 &rhs) const { return this->e || rhs.get_value(); }
	gf2_noz3 operator||(const bool rhs) const { return this->e || rhs; }
	
	gf2_noz3 operator*(const gf2_noz3 &rhs) const { return this->e && rhs.get_value(); }
	gf2_noz3 operator*(const bool rhs) const { return this->e && rhs; }
	
	gf2_noz3 operator+(const gf2_noz3 &rhs) const { return this->e != rhs.get_value(); }
	gf2_noz3 operator+(const bool rhs) const { return this->e != rhs; }
	
	gf2_noz3 &operator+=(const gf2_noz3 &rhs) { this->e = this->e != rhs.get_value(); return *this; }
	gf2_noz3 &operator+=(const bool rhs) { this->e = this->e != rhs; return *this; }
	
	gf2_noz3 operator!(void) const { return !this->e; }

    operator int(void) const { return (int)this->e; }
    // operator double(void) const { return (double)this->e; }
    operator bool(void) const { return this->e; }

	friend std::ostream &operator<<(std::ostream &os, gf2_noz3 const &m);
};

/**
 * @brief Eigen support for processing GF(2) types
 */
namespace Eigen
{
	/**
	 * @brief defining the Eigen traits for the GF(2) type
	 */
	template<> struct NumTraits<gf2_noz3>
	    : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
	{
	    typedef gf2_noz3 Real; /**< defining GF(2) as a real type */
	    typedef gf2_noz3 NonInteger; /**< defining GF(2) as a non-integer type */
	    typedef gf2_noz3 Nested; /**< defining GF(2) as a nestable type */
	    enum
	    {
	        IsComplex = 0,
	        IsInteger = 1,
	        IsSigned = 0,
	        RequireInitialization = 1,
	        ReadCost = 1,
	        AddCost = 3,
	        MulCost = 3
	    };
	};

	/**
	 * @brief defining binary operations for GF(2) types on bools
	 */
	template<> 
	struct ScalarBinaryOpTraits< gf2_noz3, bool >
	{
		enum { Defined = 1 };
		typedef gf2_noz3 ReturnType; /**< defining the return type of the binary operation */
	};

	/**
	 * @brief defining binary operations for bools on GF(2) types
	 * 
	 * @tparam element type
	 */
	template<> 
	struct ScalarBinaryOpTraits< bool, gf2_noz3 >
	{
		enum { Defined = 1 };
		typedef gf2_noz3 ReturnType; /**< defining the return type of the binary operation */
	};
}


#endif /* GF2_NOZ3_H */
