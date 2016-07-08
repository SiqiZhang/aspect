/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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



#ifndef __aspect__heating_model_impact_h
#define __aspect__heating_model_impact_h

#include <aspect/simulator_access.h>
#include <aspect/heating_model/interface.h>

#include <deal.II/base/parsed_function.h>

//TODO: Remove the dependency on the old impact code
#include <aspect/impact.h>

namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a heating model based on radioactive decay.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class Impact : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Impact ();


        /**
         * Function to compute the heating terms introduced by impact events
         */
        virtual
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual
        void
        initialize ();

        /**
         * A function that is called at the beginning of each time step.
         * This function updates the active impactor list during each 
         * timestep.
         */
        virtual
        void
        update();

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        double deltaT(const Point<dim> &p, 
                      double pressure, 
                      double temperature, 
                      const std::vector<double> &compositional_fields) const;
        struct Impact_Data
        {
          Point<3>  position;
          double    time;
          double    radius;
          double    velocity;
        };
        std::vector<struct Impact_Data> Impacts_all;
        std::vector<struct Impact_Data> Impacts_active;
        double R0;
        aspect::impact::internal::Rotate rotate;
        Tensor<1,2> surface_point_one;
        Tensor<1,2> surface_point_two;
        std::string filename;
        double      super_solidus_ratio;
    };
  }
}


#endif
