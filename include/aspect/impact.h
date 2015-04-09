#ifndef __aspect_impact_h
#define __aspect_impact_h

#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx1x/array.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/tensor.h>

using namespace dealii;
namespace aspect
{
  namespace impact
  {
    namespace internal
    {
      using namespace dealii;
      class Rotate
      {
        public:
          Rotate(){};

          void build_rotation_matrix(const Tensor<1,2> &surface_point_one,
                                     const Tensor<1,2> &surface_point_two);

          template <int dim>
          Point<3> get_new_coord(const Point<dim> &p) const;

          template <int dim>
          void screen_output(const Tensor<1,2> &surface_point_one,
                             const Tensor<1,2> &surface_point_two,
                             const ConditionalOStream &pcout) const;
        private:
          /**
           * Roation function copied from glpate module of ASPECT
           */

          Tensor<1,3>
          rotate_around_axis (const Tensor<1,3> &position,
                              const Tensor<1,3> &rotation_axis,
                              const double angle) const;

          std_cxx1x::array<double,3>
          angles_from_matrix (const Tensor<2,3> &rotation_matrix) const;

          double
          rotation_axis_from_matrix (Tensor<1,3> &rotation_axis,
                                     const Tensor<2,3> &rotation_matrix) const;

          Tensor<2,3>
          rotation_matrix_from_axis (const Tensor<1,3> &rotation_axis,
                                     const double rotation_angle) const;

          template <int in, int out>
          Tensor<1,out> convert_tensor (const Tensor<1,in> &old_tensor) const;

          Tensor<1,3>
          cartesian_surface_coordinates(const Tensor<1,3> &sposition) const;
          Tensor<2,3> rotation_matrix;
      };
    }

    template<int dim>
    class ImpactFunction
    {
      public:
        ImpactFunction();
        void Initialize(const std::string &filename,
                        double time0, 
                        double R0,
                        const Tensor<1,2> &surface_point_one,
                        const Tensor<1,2> &surface_point_two,
                        ConditionalOStream &pcout);
        double deltaT(const Point<dim> &p, double pressure, double temperature);
        void update_time(double time);
      private:
        struct Impact_Data
        {
          //Point<dim>  position;
          Point<3>  position;
          double      time;
          double      radius;
          double      velocity;
        };
        std::vector<struct Impact_Data> Impacts_all;
        std::vector<struct Impact_Data> Impacts_active;
        double old_time;
        double R0;
        internal::Rotate rotate;
    };
  }
}
#endif
