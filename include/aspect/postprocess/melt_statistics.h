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
/*  $Id: temperature_statistics.h 1433 2012-12-08 08:24:55Z bangerth $  */


#ifndef __aspect__postprocess__melt_statistics_h
#define __aspect__postprocess__melt_statistics_h

#include <vector>

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the temperature.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class MeltStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some temperature statistics.
         **/
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

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
        
        static Tensor<1,dim> transfer_coord(const Point<dim> &p);

        /** 
         * Calculate matrix compression at given point
         **/
        double get_compression(const Point<dim> &p) const;
        
        /**
         * Additional grid to store melt fraction and matrix compression
         **/
        //TODO: Only works in 2D at the moment
        class MeltGrid
        {
          public:
            MeltGrid();

            void set_grid(double R0, double R1, unsigned nr, unsigned nh);
            
            void add_melt(double h0, double h1, double r0, double r1, double f);

            bool is_grid_set() const;

            void calculate_compression();

            double get_compression(const Point<dim> &p) const;

            void output_vtk(unsigned int time_step_num) const;

            void sum_melt_fraction_parallel(const MPI_Comm &  mpi_communicator);
            
            struct MeltCell
            {
              double r0;
              double r1;
              double h0;
              double h1;
            };
          private:
            double               R0;
            double               R1;
            double               dr;
            double               dh;
            unsigned             nr;
            unsigned             nh;
            std::vector<double>  melt_fraction;
            std::vector<double>  matrix_compression;
            bool                 grid_set;

            struct MeltCell get_overlapping_cell(struct MeltCell a, struct MeltCell b) const;

            double modify_radian(double ref_angle, double angle) const;

            double get_volume(struct MeltCell a) const;

            unsigned cell_index(unsigned i_h, unsigned i_r) const
            {
              return i_h*nr+i_r;
            }

            unsigned point_index(unsigned i_h, unsigned i_r) const
            {
              return i_h*(nr+1)+i_r;
            }

        };

      private:
        MeltGrid melt_grid;
        int       melt_grid_nh;
        int       melt_grid_nr;
        bool      melt_grid_output;
    };
  }
}


#endif
