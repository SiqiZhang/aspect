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


#include <aspect/termination_criteria/end_walltime.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    time_t
    EndWalltime<dim>::start_walltime = std::time(NULL);

    template <int dim>
    bool
    EndWalltime<dim>::execute()
    {
      // Measure the time at process 0 to detarmin termination.
      // This makes all the processors terminate together.
      int global_do_terminate = ((std::time(NULL)-start_walltime) >= end_wall_time);
      MPI_Bcast(&global_do_terminate, 1, MPI_INT, 0, this->get_mpi_communicator());
      return (global_do_terminate == 1);
    }


    template <int dim>
    void
    EndWalltime<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        prm.declare_entry ("Wall time",
                           "0",
                           Patterns::Double (),
                           "The wall time of the simulation. The default value is zero, "
                           "which disabled this termination criteria. Unit: seconds.");
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    EndWalltime<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        end_wall_time = prm.get_double ("Wall time");
        start_walltime = std::time(NULL);
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace TerminationCriteria
  {
    ASPECT_REGISTER_TERMINATION_CRITERION(EndWalltime,
                                          "wall time",
                                          "Terminate the simulation once the wall time limit has reached.");
  }
}
