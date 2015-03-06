#ifndef __aspect_impact_h
#define __aspect_impact_h

#include <deal.II/base/point.h>


using namespace dealii;
namespace aspect
{
	template<int dim>
	class ImpactFunction
	{
	public:
		ImpactFunction(){};
		void Initialize(const std::string &filename, double time0, double R0);
		double deltaT(const Point<dim> &p, double pressure, double temperature);
		void update_time(double time);
	private:
		struct Impact_Data
		{
			Point<dim>  position;
			double      time;
			double      radius;
			double      velocity;
		};
		std::vector<struct Impact_Data> Impacts_all;
		std::vector<struct Impact_Data> Impacts_active;
		double old_time;
    double R0;
	};
};
#endif
