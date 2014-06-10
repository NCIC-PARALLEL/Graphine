#ifndef _GRE_Util_Stat_HPP
#define _GRE_Util_Stat_HPP
#include <cstring>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
namespace GRE{
	namespace Util{
	class Statistics{
	public:
		enum StatType {Stat_Rate, Stat_Value};
	public:
		template <typename T>
		void calculate(std::vector<T>& data_)
		{
			calculate(data_, 0, data_.size());
		}
		template <typename T>
		void calculate(std::vector<T>& data_, const int begin, const int n)
		{
			assert(begin+n <= data_.size());
			std::vector<double> data;
			data.resize(n);
			for(int i=0; i<n; i++){
				data[i]=static_cast<double>(data_[begin+i]);
			}
			calculate1(data);
		}
		template <typename T>
		void calculate(T* data_, const int begin, const int n)
		{
			std::vector<double> data;
			data.resize(n);
			for(int i=0; i<n; i++){
				data[i]=static_cast<double>(data_[begin+i]);
			}
			calculate1(data);
		}
		//
		void calculate1(std::vector<double>& data)
		{
			const int n=data.size();
			//
			long double s, mean;
			double t;
			int k;
			//
			stats.clear();
			stats.resize(nStats);
			/* Quartiles */
			sort(data.begin(), data.end());
		  	stats[0] = data[0];
			t = (n+1) / 4.0;
		  	k = (int) t;
			if (t == k)
		    	stats[1] = data[k];
			else
    			stats[1] = 3*(data[k]/4.0) + data[k+1]/4.0;
	  		t = (n+1) / 2.0;
  			k = (int) t;
  			if (t == k)
    			stats[2] = data[k];
	  		else
    			stats[2] = data[k]/2.0 + data[k+1]/2.0;
			t = 3*((n+1) / 4.0);
			k = (int) t;
			if (t == k)
			    stats[3] = data[k];
			else
				stats[3] = data[k]/4.0 + 3*(data[k+1]/4.0);
	  		stats[4] = data[n-1];
	
			s = data[n-1];
			for (k = n-1; k > 0; --k)
    			s += data[k-1];
			mean = s/n;
			stats[5] = mean;
			s = data[n-1] - mean;
			s *= s;
			for (k = n-1; k > 0; --k) {
    			long double tmp = data[k-1] - mean;
    			s += tmp * tmp;
  			}
			stats[6] = sqrt (s/(n-1));

  			s = (data[0]? 1.0L/data[0] : 0);
  			for (k = 1; k < n; ++k)
		    	s += (data[k]? 1.0L/data[k] : 0);
  			stats[7] = n/s;
	  		mean = s/n;

		  /* 
		   * Nilan Norris, The Standard Errors of the Geometric and Harmonic
		   * Means and Their Application to Index Numbers, 1940.
		   * http://www.jstor.org/stable/2235723
		   * */
			s = (data[0]? 1.0L/data[0] : 0) - mean;
			s *= s;
			for (k = 1; k < n; ++k) {
				long double tmp = (data[k]? 1.0L/data[k] : 0) - mean;
				s += tmp * tmp;
			}
			s = (sqrt (s)/(n-1)) * stats[7] * stats[7];
			stats[8] = s;
		}
		void printAll(StatType st = Stat_Rate)
		{
    		printf ("min: %20.17e\n", stats[0]);	
    		printf ("firstquartile: %20.17e\n", stats[1]);
    		printf ("median: %20.17e\n", stats[2]);
    		printf ("thirdquartile: %20.17e\n", stats[3]);
    		printf ("max: %20.17e\n", stats[4]);
    		if (st==Stat_Rate) {
      			printf ("mean: %20.17e\n", stats[5]);
      			printf ("stddev: %20.17e\n", stats[6]);
    		} else {
      			printf ("harmonic_mean: %20.17e\n", stats[7]);
      			printf ("harmonic_stddev: %20.17e\n", stats[8]);
    		}
		}
		double getMin() const
		{
			return stats[0];
		}
		double getMax() const
		{
			return stats[4];
		}
		double get1stQuartile() const
		{
			return stats[1];
		}
		double getMedian() const
		{
			return stats[2];
		}
		double get3rdQuartile() const
		{
			return stats[3];
		}
		double getMean() const
		{
			return stats[5];
		}
		double getHarmonicMean() const
		{
			return stats[7];
		}
		double getStdDev() const
		{
			return stats[6];
		}
		double getHarmonicStdDev() const
		{
			return stats[8];
		}
	private:
		static const int nStats = 9; 
		std::vector<double> stats; 
	};
	}//end Util
}
#endif

