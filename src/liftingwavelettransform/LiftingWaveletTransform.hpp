/*
    LiftingWaveletTransform
    Copyright (C) 2009-2013 German Tischler
    Copyright (C) 2011-2013 Genome Research Limited

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#if ! defined(LIFTINGWAVELETTRANSFORM_HPP)
#define LIFTINGWAVELETTRANSFORM_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

/**
 * class containing lifting scheme implementations of the CDF 9/7 and 5/3 filters
 * plus some tests for these filters
 **/
struct LiftingWaveletTransform
{
	/**
	 * array accessor with mirroring at borders
	 **/
	template<typename _data_iterator>
	struct MirrorAccessor
	{
		typedef _data_iterator data_iterator;
		typedef typename std::iterator_traits<data_iterator>::value_type value_type;
		
		data_iterator const A;
		int64_t n;
		
		MirrorAccessor(data_iterator rA, uint64_t const rn) : A(rA), n(rn) {}
		
		value_type operator()(int64_t i, int64_t o) const
		{
			if ( n == 1 )
				return A[0];
				
			int64_t const j = i + o;
		
			if ( j >= 0 && j < n )
				return A[j];	

			while ( o != 0 )
			{
				if ( i == 0 && o < 0 )
				{
					o = -o;
					i += 1;
					o--;
				}
				else if ( i == static_cast<int64_t>(n-1) && o > 0 )
				{
					o--;
					o = -o;
					i -= 1;
				}
				else if ( o > 0 )
				{
					i += 1;
					o -= 1;
				}
				else
				{
					i -= 1;
					o += 1;
				}
			}
			
			return A[i];
		}
	};

	/**
	 * lifting prediction step
	 **/
	template<typename data_iterator>
	static void predict(data_iterator A, uint64_t const n, typename ::std::iterator_traits<data_iterator>::value_type const coeff)
	{
		for ( uint64_t i = 1; i+1 < n; i += 2 )
			A[i] += coeff * (A[i-1]+A[i+1]);
		
		if ( (n & 1) == 0 && n > 1 )
			A[n-1] += (2.0 * coeff) * A[n-2];
	}

	/**
	 * lifting inverse prediction step
	 **/
	template<typename data_iterator>
	static void ipredict(data_iterator const A, uint64_t const n, typename ::std::iterator_traits<data_iterator>::value_type const coeff)
	{
		if ( (n & 1) == 0 && n > 1 )
			A[n-1] -= (2.0 * coeff) * A[n-2];

		for ( uint64_t i = 1; i+1 < n; i += 2 )
			A[i] -= coeff * (A[i-1]+A[i+1]);
	}

	/**
	 * lifting update step
	 **/
	template<typename data_iterator>
	static void update(data_iterator const A, uint64_t const n, typename ::std::iterator_traits<data_iterator>::value_type const coeff)
	{
		for ( uint64_t i = 2; i+1 < n; i += 2 )
			A[i] += coeff * (A[i-1]+A[i+1]);

		if ( n > 1 )
		{
			A[0] += (2.0 * coeff) * A[1];
			
			if ( (n & 1) == 1 )
				A[n-1] += (2.0*coeff) * A[n-2];
		}
	}

	/**
	 * lifting inverse update step
	 **/
	template<typename data_iterator>
	static void iupdate(data_iterator const A, uint64_t const n, typename ::std::iterator_traits<data_iterator>::value_type const coeff)
	{
		if ( n > 1 )
		{
			if ( (n & 1) == 1 )
				A[n-1] -= (2.0*coeff) * A[n-2];
				
			A[0] -= (2.0 * coeff) * A[1];
		}

		for ( uint64_t i = 2; i+1 < n; i += 2 )
			A[i] -= coeff * (A[i-1]+A[i+1]);	
	}

	/**
	 * lifting scaling step
	 **/
	template<typename data_iterator>
	static void scale(data_iterator const A, uint64_t const n, 
		typename ::std::iterator_traits<data_iterator>::value_type const lmult,
		typename ::std::iterator_traits<data_iterator>::value_type const hmult		
	)
	{
		for ( uint64_t i = 0; i < n; i += 2 )
			A[i] *= lmult;
		for ( uint64_t i = 1; i < n; i += 2 )
			A[i] *= hmult;
	}

	/**
	 * lifting inverse scaling step
	 **/
	template<typename data_iterator>
	static void iscale(
		data_iterator const A, uint64_t const n, 
		typename ::std::iterator_traits<data_iterator>::value_type const lmult,
		typename ::std::iterator_traits<data_iterator>::value_type const hmult
	)
	{
		typename ::std::iterator_traits<data_iterator>::value_type const ilmult = 1.0 / lmult;
		typename ::std::iterator_traits<data_iterator>::value_type const ihmult = 1.0 / hmult;
		
		for ( uint64_t i = 0; i < n; i += 2 )
			A[i] *= ilmult;
		for ( uint64_t i = 1; i < n; i += 2 )
			A[i] *= ihmult;
	}

	/**
	 * interleave array elements
	 **/
	template<typename data_iterator>
	static void interleave(data_iterator const A, uint64_t const n, bool const quick_reorder = false)
	{
		// quick reordering using additional space (linear time)
		if ( quick_reorder )
		{
			std::vector< typename std::iterator_traits<data_iterator>::value_type > P(n);
			for ( uint64_t i = 0; i < (n+1)/2; ++i )
				P[2*i] = A[i];
			for ( uint64_t i = 0; i < (n+0)/2; ++i )
				P[2*i+1] = A[(n+1)/2+i];
				
			std::copy(P.begin(),P.end(),A);		
		}
		// in place reordering without additional space (n log n time)
		else
		{
			uint64_t blocksize = 1;
			while ( (blocksize<<2) < n )
				blocksize <<= 1;

			for ( ; blocksize; blocksize >>= 1)
			{
				uint64_t low = 0;
				
				while ( low+2*blocksize < n )
				{
					uint64_t const a = low;
					uint64_t const b = a+blocksize;
					uint64_t const c = b+blocksize;
					uint64_t const r = std::min(blocksize<<1,n-c);
					uint64_t const r_0 = (r+1)/2;
					uint64_t const r_1 = r-r_0;
					
					uint64_t const d = c + r_0;
					uint64_t const e = d + r_1;
					
					std::reverse(A+b,A+d);
					std::reverse(A+c,A+d);			
					std::reverse(A+b,A+c);
					
					low = e;
				}
			}
		}
	}

	/**
	 * deinterleave array elements
	 **/
	template<typename data_iterator>
	static void deinterleave(data_iterator const A, uint64_t const n, bool const quick_reorder = false)
	{
		// quick reordering using additional space (linear time)
		if ( quick_reorder )
		{
			std::vector< typename std::iterator_traits<data_iterator>::value_type > P(n);
			for ( uint64_t i = 0; i < (n+1)/2; ++i )
				P [ i ] = A[2*i];
			for ( uint64_t i = 0; i < (n+0)/2; ++i )
				P [ i+(n+1)/2 ] = A[2*i+1];
			std::copy(P.begin(),P.end(),A);
		}
		// in place reordering without additional space (n log n time)
		else
		{
			for ( uint64_t blocksize = 1; 2*blocksize < n; blocksize <<= 1 )
			{
				uint64_t low = 0;
				
				while ( low+2*blocksize < n )
				{
					uint64_t const a = low;
					uint64_t const b = a+blocksize;
					uint64_t const c = b+blocksize;
					uint64_t const r = std::min(blocksize<<1,n-c);
					uint64_t const r_0 = (r+1)/2;
					uint64_t const r_1 = r-r_0;
					uint64_t const d = c + r_0;
					uint64_t const e = d + r_1;
					
					std::reverse(A+b,A+c);
					std::reverse(A+c,A+d);
					std::reverse(A+b,A+d);
					
					low = e;
				}
			}
		}
	}

	/**
	 * Cohen, Daubechies, Feauveau 9/7 tap filter implemented using Swelden's lifting scheme
	 *
	 * runtime is linear if quick_reorder==true, n log n otherwise
	 **/
	template<typename data_iterator>
	static void cdf97(data_iterator const A, uint64_t const n, bool quick_reorder = false)
	{
		predict<data_iterator>(A,n,-1.586134342);
		update<data_iterator>(A,n,-0.05298011854);
		predict<data_iterator>(A,n,0.8829110762);
		update<data_iterator>(A,n,0.4435068522);
		scale<data_iterator>(A,n,1.0 / 1.23017410558578,1.0 / 1.62578613134411);
		deinterleave(A,n,quick_reorder);
	}

	/**
	 * Inverse Cohen, Daubechies, Feauveau 9/7 tap filter implemented using Swelden's lifting scheme
	 *
	 * runtime is linear if quick_reorder==true, n log n otherwise
	 **/
	template<typename data_iterator>
	static void icdf97(data_iterator const A, uint64_t const n, bool const quick_reorder = false)
	{
		interleave(A,n,quick_reorder);
		iscale<data_iterator>(A,n,1.0 / 1.23017410558578,1.0 / 1.62578613134411);
		iupdate<data_iterator>(A,n,0.4435068522);
		ipredict<data_iterator>(A,n,0.8829110762);
		iupdate<data_iterator>(A,n,-0.05298011854);
		ipredict<data_iterator>(A,n,-1.586134342);
	}

	/**
	 * Cohen, Daubechies, Feauveau 5/3 tap filter implemented using Swelden's lifting scheme
	 *
	 * runtime is linear if quick_reorder==true, n log n otherwise
	 **/
	template<typename data_iterator>
	static void cdf53(data_iterator const A, uint64_t const n, bool quick_reorder = false)
	{
		predict<data_iterator>(A,n,-0.5);
		update<data_iterator>(A,n,0.25);
		deinterleave(A,n,quick_reorder);
	}

	/**
	 * Inverse Cohen, Daubechies, Feauveau 5/3 tap filter implemented using Swelden's lifting scheme
	 *
	 * runtime is linear if quick_reorder==true, n log n otherwise
	 **/
	template<typename data_iterator>
	static void icdf53(data_iterator const A, uint64_t const n, bool const quick_reorder = false)
	{
		interleave(A,n,quick_reorder);
		iupdate<data_iterator>(A,n,0.25);
		ipredict<data_iterator>(A,n,-0.5);
	}
	
	/**
	 * 5 tap low pass coefficients
	 **/
	static std::vector<double> get53LowPassFilter()
	{
		std::vector<double> V(5);
		V[0] = V[4] = -1.0/8.0;
		V[1] = V[3] =  2.0/8.0;
		V[2] =         6.0/8.0;
		return V;
	}

	/**
	 * 3 tap high pass coefficients
	 **/
	static std::vector<double> get53HighPassFilter()
	{
		std::vector<double> V(3);
		V[0] = V[2] = -1.0/2.0;
		V[1] =         1.0/1.0;
		return V;
	}

	/**
	 * 9 tap low pass coefficients
	 **/
	static std::vector<double> get97LowPassFilter()
	{
		std::vector<double> V(9);
		
		V[0] = V[8] =  .02674875741080976000;
		V[1] = V[7] = -.01686411844287495000;
		V[2] = V[6] = -.07822326652898785000;
		V[3] = V[5] =  .26686411844287230000;
		V[4] =         .60294901823635790000;
		
		return V;
	}

	/**
	 * 7 tap low pass coefficients
	 **/
	static std::vector<double> get97HighPassFilter()
	{
		std::vector<double> V(7);
		
		V[0] = V[6] =  .04563588155712474000;
		V[1] = V[5] = -.02877176311424978500;
		V[2] = V[4] = -.29563588155712350000;
		V[3] =         .55754352622849700000;
		
		return V;
	}
	
	/**
	 * generic convolution function (generating single point at i)
	 **/
	template<typename accessor_type>
	static typename accessor_type::value_type convolve(
		accessor_type const & M,
		std::vector<double> const & V,
		uint64_t const i
	)
	{
		typedef typename accessor_type::value_type value_type;
		uint64_t const v2 = V.size()/2;
	
		value_type cval = value_type();
		for ( int64_t j = 0; j < V.size(); ++j )
			cval += M(i,j-v2) * V[j];
			
		return cval;
	}
	
	/**
	 * generic convolution function using additional array
	 **/
	template<typename data_iterator>
	static void convolve(data_iterator const A, uint64_t const n, std::vector<double> const & filter)
	{
		typedef typename std::iterator_traits<data_iterator>::value_type value_type;
		std::vector<value_type> B(n);
		MirrorAccessor< data_iterator > M(A,n);
		
		for ( uint64_t i = 0; i < n; ++i )
			B[i] = convolve(M,filter,i);
			
		std::copy(B.begin(),B.end(),A);
	}
	
	template<typename data_iterator>
	static std::vector< typename std::iterator_traits<data_iterator>::value_type > filter53Low(data_iterator const A, uint64_t const n)
	{
		std::vector< typename std::iterator_traits<data_iterator>::value_type > V(A,A+n);
		convolve(V.begin(),V.size(),get53LowPassFilter());
		return V;
	}

	template<typename data_iterator>
	static std::vector< typename std::iterator_traits<data_iterator>::value_type > filter53High(data_iterator const A, uint64_t const n)
	{
		std::vector< typename std::iterator_traits<data_iterator>::value_type > V(A,A+n);
		convolve(V.begin(),V.size(),get53HighPassFilter());
		return V;
	}

	template<typename data_iterator>
	static std::vector< typename std::iterator_traits<data_iterator>::value_type > filter97Low(data_iterator const A, uint64_t const n)
	{
		std::vector< typename std::iterator_traits<data_iterator>::value_type > V(A,A+n);
		convolve(V.begin(),V.size(),get97LowPassFilter());
		return V;
	}

	template<typename data_iterator>
	static std::vector< typename std::iterator_traits<data_iterator>::value_type > filter97High(data_iterator const A, uint64_t const n)
	{
		std::vector< typename std::iterator_traits<data_iterator>::value_type > V(A,A+n);
		convolve(V.begin(),V.size(),get97HighPassFilter());
		return V;
	}
	
	template<typename data_iterator>
	static void haar(data_iterator const A, uint64_t const n)
	{
		typedef typename std::iterator_traits<data_iterator>::value_type value_type;
	
		for ( uint64_t i = 0; i < n; i += 2 )
		{
			value_type const a0 = A[i];
			value_type const a1 = A[i+1];
			
			A[i]   = (a0+a1)/2;
			A[i+1] = (a0-a1)/2;
		}
	}
	
	template<typename container_type>
	static container_type haar(container_type const & C)
	{
		container_type R(C.begin(),C.end());
		haar(R.begin(),R.size());
		return R;
	}

	template<typename data_iterator>
	static void ihaar(data_iterator const A, uint64_t const n)
	{
		typedef typename std::iterator_traits<data_iterator>::value_type value_type;
	
		for ( uint64_t i = 0; i < n; i += 2 )
		{
			value_type const a0 = A[i];
			value_type const a1 = A[i+1];
			
			A[i]   = (a0 + a1);
			A[i+1] = (a0 - a1);
		}
	}

	template<typename container_type>
	static container_type ihaar(container_type const & C)
	{
		container_type R(C.begin(),C.end());
		ihaar(R.begin(),R.size());
		return R;
	}
	
	/*
	 * compare 5/3 filter lifting implementation with convolution based method
	 */
	template<typename data_iterator>
	static bool testConvolution53(data_iterator const A, uint64_t const n)
	{
		typedef typename std::iterator_traits<data_iterator>::value_type value_type;
		std::vector<value_type> R(A,A+n);

		LiftingWaveletTransform::cdf53(A,n);
				
		MirrorAccessor< typename std::vector<value_type>::const_iterator > M(R.begin(),n);
		
		bool ok = true;
		
		std::vector<double> const L = get53LowPassFilter();
		for ( uint64_t i = 0; i < n; i += 2 )
		{
			value_type cval = convolve(M,L,i);
			value_type const lval = A[i/2];
			
			if ( std::abs(cval-lval) > 1e-4 )
			{
				std::cerr << "L[" << i << "]=(" << cval << "," << lval << ")" << std::endl;
				ok = false;
			}
		}

		std::vector<double> const H = get53HighPassFilter();
		for ( uint64_t i = 1; i < n; i += 2 )
		{
			value_type const cval = convolve(M,H,i);
			value_type const lval = A[ (n+1)/2 + i/2];

			if ( std::abs(cval-lval) > 1e-4 )
			{
				std::cerr << "H[" << i << "]=(" << cval << "," << lval << ")" << std::endl;
				ok = false;
			}
		}

		return ok;
	}

	/*
	 * compare 9/7 filter lifting implementation with convolution based method
	 */
	template<typename data_iterator>
	static bool testConvolution97(data_iterator A, uint64_t const n)
	{
		typedef typename std::iterator_traits<data_iterator>::value_type value_type;
		std::vector<value_type> R(A,A+n);

		LiftingWaveletTransform::cdf97(A,n);
			
		MirrorAccessor< typename std::vector<value_type>::const_iterator > M(R.begin(),n);

		bool ok = true;

		std::vector<double> const L = get97LowPassFilter();
		for ( uint64_t i = 0; i < n; i += 2 )
		{	
			value_type const cval = convolve(M,L,i);
			value_type const lval = A[i/2];

			if ( std::abs(cval-lval) > 1e-4 )
			{
				std::cerr << "L[" << i << "]=(" << cval << "," << lval << ")" << std::endl;
				ok = false;
			}
		}

		std::vector<double> const H = get97HighPassFilter();
		uint64_t const h2 = H.size()/2;
		for ( uint64_t i = 1; i < n; i += 2 )
		{
			value_type const cval = convolve(M,H,i);
			value_type const lval = A[ (n+1)/2 + i/2];
			
			if ( std::abs(cval-lval) > 1e-4 )
			{
				std::cerr << "H[" << i << "]=(" << cval << "," << lval << ")" << std::endl;
				ok = false;
			}
		}
		
		return ok;
	}

	/**
	 * compare results for lifting and convolution based computations for two vectors
	 **/
	static bool testConvolution()
	{
		bool ok = true;

		double A[] = { 1,2,3,4,3,2,1,0,7,5,6,2 };
		ok = ok && testConvolution53(&A[0],sizeof(A)/sizeof(A[0]));
		ok = ok && testConvolution97(&A[0],sizeof(A)/sizeof(A[0]));
		double B[] = { 1,2,3,4,3,2,1,0,7,5,6,2,11 };
		ok = ok && testConvolution53(&B[0],sizeof(B)/sizeof(B[0]));
		ok = ok && testConvolution97(&B[0],sizeof(B)/sizeof(B[0]));
		
		return ok;
	}

	/**
	 * run forward and reverse tests with sets of random numbers
	 **/
	static bool testRandom()
	{
		bool ok = true;
		
		/*
		 * generate random double vectors of length n
		 * and process them with forward and reverse transforms
		 * after each cycle compare them to the original vector
		 */
		for ( uint64_t n = 1; n <= 128; ++n )
		{
			for ( uint64_t i = 0; i < 1024; ++i )
			{
				std::vector<double> V;
				for ( uint64_t j = 0; j < n; ++j )
					V.push_back(rand());

				std::vector<double> const R = V;
					
				LiftingWaveletTransform::cdf53(V.begin(),V.size());
				LiftingWaveletTransform::icdf53(V.begin(),V.size());
				
				for ( uint64_t j = 0; j < n; ++j )
					if ( std::abs(V[j]-R[j]) > 1e-4 )
					{
						std::cerr << std::abs(V[j]-R[j]) << std::endl;
						ok = false;
					}

				LiftingWaveletTransform::cdf53(V.begin(),V.size(),true);
				LiftingWaveletTransform::icdf53(V.begin(),V.size(),true);
				
				for ( uint64_t j = 0; j < n; ++j )
					if ( std::abs(V[j]-R[j]) > 1e-4 )
					{
						std::cerr << std::abs(V[j]-R[j]) << std::endl;
						ok = false;
					}

				LiftingWaveletTransform::cdf97(V.begin(),V.size());
				LiftingWaveletTransform::icdf97(V.begin(),V.size());

				for ( uint64_t j = 0; j < n; ++j )
					if ( std::abs(V[j]-R[j]) > 1e-4 )
					{
						std::cerr << std::abs(V[j]-R[j]) << std::endl;
						ok = false;
					}

				LiftingWaveletTransform::cdf97(V.begin(),V.size(),true);
				LiftingWaveletTransform::icdf97(V.begin(),V.size(),true);

				for ( uint64_t j = 0; j < n; ++j )
					if ( std::abs(V[j]-R[j]) > 1e-4 )
					{
						std::cerr << std::abs(V[j]-R[j]) << std::endl;
						ok = false;
					}

				LiftingWaveletTransform::cdf97(V.begin(),V.size(),false);
				LiftingWaveletTransform::icdf97(V.begin(),V.size(),true);

				for ( uint64_t j = 0; j < n; ++j )
					if ( std::abs(V[j]-R[j]) > 1e-4 )
					{
						std::cerr << std::abs(V[j]-R[j]) << std::endl;
						ok = false;
					}

				LiftingWaveletTransform::cdf97(V.begin(),V.size(),true);
				LiftingWaveletTransform::icdf97(V.begin(),V.size(),false);

				for ( uint64_t j = 0; j < n; ++j )
					if ( std::abs(V[j]-R[j]) > 1e-4 )
					{
						std::cerr << std::abs(V[j]-R[j]) << std::endl;
						ok = false;
					}

				LiftingWaveletTransform::haar(V.begin(),V.size());
				LiftingWaveletTransform::ihaar(V.begin(),V.size());

				for ( uint64_t j = 0; j < n; ++j )
					if ( std::abs(V[j]-R[j]) > 1e-4 )
					{
						std::cerr << std::abs(V[j]-R[j]) << std::endl;
						ok = false;
					}
			}
		}
		
		return ok;
	}

	/**
	 * call tests
	 **/	
	static bool test()
	{
		bool ok = true;
		ok = ok && testConvolution();
		ok = ok && testRandom();
		return ok;	
	}
};
#endif
