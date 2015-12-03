/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/visualization/matrix_compression.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      MatrixCompression<dim>::
      MatrixCompression ()
        :
        DataPostprocessorScalar<dim> ("matrix_compression",
                                      update_q_points)
      {}



      template <int dim>
      void
      MatrixCompression<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &,
                                         const std::vector<std::vector<Tensor<1,dim> > > &,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = evaluation_points.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());

        const Postprocess::MeltStatistics<dim> * melt_statistics
          = this->template find_postprocessor<const Postprocess::MeltStatistics<dim> >();
        AssertThrow (melt_statistics!=NULL, 
            ExcMessage("The matrix compression postprocessor has to work with melt_statistics postprocessor."));
        double time_step = this->get_timestep();

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = melt_statistics->get_compression(evaluation_points[q])/time_step*year_in_seconds;
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MatrixCompression,
                                                  "matrix compression",
                                                  "A visualization output object that generates output "
                                                  "for the matrix compression due to melt extraction. Units: m/year")
    }
  }
}
