/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/melt_fraction2.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/melt.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      MeltFraction2<dim>::
      MeltFraction2 ()
        :
        DataPostprocessorScalar<dim> ("melt_fraction2",
                                      update_values | update_q_points)
      {}



      template <int dim>
      void
      MeltFraction2<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                         const std::vector<Point<dim> >                  &normals,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const MaterialModel::Melt<dim> *material_model =
            dynamic_cast<const MaterialModel::Melt<dim> *> (&(this->get_material_model()));
        AssertThrow(material_model!=0, ExcMessage("This postprocess only works with melt material model."));
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const double pressure    = uh[q][this->introspection().component_indices.pressure];
            const double temperature = uh[q][this->introspection().component_indices.temperature];
            double corrected_temperature,
                   corrected_pressure;
            std::vector<double> composition(this->n_compositional_fields());

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              composition[c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
            //if(pressure<10.e9 && pressure>0)
            // Correct the temperature and pressure for incompressible case
            if(material_model->is_compressible())
            {
              corrected_temperature = temperature;
              corrected_pressure    = pressure;
            }
            else
            {
              corrected_temperature = temperature + this->get_adiabatic_conditions().temperature(evaluation_points[q])
                                                  - this->get_adiabatic_surface_temperature();
              corrected_pressure    = this->get_adiabatic_conditions().pressure(evaluation_points[q]);
            }
              
            if(pressure>0)
              //computed_quantities[q](0) = material_model->melt_fraction(corrected_temperature,corrected_pressure,composition,evaluation_points[q]);
              computed_quantities[q](0) = material_model->melt_extraction(corrected_temperature,corrected_pressure,composition,evaluation_points[q]);
            else
              computed_quantities[q](0) = 0.;
          }
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltFraction2,
                                                  "melt fraction2", // TODO write down equations here
                                                  "A visualization output object that generates output "
                                                  "for the melt fraction at the temperature and "
                                                  "pressure of the current point (batch melting). "
                                                  "Does not take into account latent heat. "
                                                  "If there are no compositional fields, this postprocessor "
                                                  "will visualize the melt fraction of peridotite "
                                                  "(calculated using the anhydrous model of Katz, 2003). ")
    }
  }
}

