/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id: heat_flux_statistics.h 1433 2012-12-08 08:24:55Z bangerth $  */


#ifndef __aspect__postprocess_dynamic_core_statistics_h
#define __aspect__postprocess_dynamic_core_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the dynamic core.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class DynamicCoreStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some heat_flux statistics.
         **/
        DynamicCoreStatistics();

        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        double get_CMB_heat_flux() const;
        void   set_CMB_heat_flux(double heat_flux);
        bool   is_initialized() const;
      private:
        double CMB_heat_flux;
        bool   initialized;
    };
  }
}


#endif
