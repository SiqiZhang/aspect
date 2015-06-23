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


#include <aspect/gravity_model/layer_density.h>

#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    Tensor<1,dim>
    LayerDensity<dim>::gravity_vector (const Point<dim> &p) const
    {
      double mass = 0.;
      // Gravity constant
      const double G = 6.674e-11;
      const double r = p.norm();
      for(unsigned i=0;i<n_layers;i++)
      {
        if(r>radiuses[i])
        {
          if(i==0)
            mass += 4./3.*M_PI*densities[i]*radiuses[i]*radiuses[i]*radiuses[i];
          else
            mass += 4./3.*M_PI*densities[i]*(radiuses[i]*radiuses[i]*radiuses[i]
                               -radiuses[i-1]*radiuses[i-1]*radiuses[i-1]);
        }
        else
        {
          if(i==0)
            mass += 4./3.*M_PI*r*r*r*densities[i];
          else
            mass += 4./3.*M_PI*densities[i]*(r*r*r-radiuses[i-1]*radiuses[i-1]*radiuses[i-1]);
        }
      }
      return - (G * mass / r / r) * p/r;
    }



    template <int dim>
    void
    LayerDensity<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Layer density");
        {
          prm.declare_entry ("Number of layers", "0",
                             Patterns::Integer (0),
                             "Number of density layers. ");
          prm.declare_entry ("Radiuses","",
                             Patterns::List(Patterns::Double(0)),
                             "Radiuses of different layers, from inner to outter "
                             "(Unit: m)");
          prm.declare_entry ("Densities","",
                             Patterns::List(Patterns::Double(0)),
                             "Densities of different layers, from inner to outter "
                             "(Unit: kg/m^3");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    LayerDensity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Layer density");
        {
          n_layers  = prm.get_integer ("Number of layers");
          radiuses  = dealii::Utilities::string_to_double
                      (dealii::Utilities::split_string_list(prm.get("Radiuses")));
          densities = dealii::Utilities::string_to_double
                      (dealii::Utilities::split_string_list(prm.get("Densities")));
          AssertThrow (radiuses.size() == n_layers && densities.size() == n_layers,
                       ExcMessage ("Number of layers doesn't match the data in Radiuses and Densities."));
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(LayerDensity,
                                  "layer density",
                                  "A gravity model in which the gravity direction is radially inward "
                                  "and the magnitude is calculated with 1D layered density profile. ")

  }
}
