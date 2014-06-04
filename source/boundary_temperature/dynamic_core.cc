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



#include <aspect/boundary_temperature/dynamic_core.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/dynamic_core_statistics.h>
#include <aspect/simulator.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {

    template <int dim>
    double
    Dynamic_core<dim>::
    temperature (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &location) const
    {
      // verify that the geometry is in fact a spherical shell since only
      // for this geometry do we know for sure what boundary indicators it
      // uses and what they mean
      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*>(&geometry_model)
              != 0,
              ExcMessage ("This boundary model is only implemented if the geometry is "
                          "in fact a spherical shell."));

      switch (boundary_indicator)
        {
          case 0:
            return inner_temperature;
          case 1:
            return outer_temperature;
          default:
            Assert (false, ExcMessage ("Unknown boundary indicator."));
            return std::numeric_limits<double>::quiet_NaN();
        }
    }


    template <int dim>
    double
    Dynamic_core<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return std::min (inner_temperature, outer_temperature);
    }



    template <int dim>
    double
    Dynamic_core<dim>::
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return std::max (inner_temperature, outer_temperature);
    }



    template <int dim>
    void
    Dynamic_core<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Dynamic core");
        {
          prm.declare_entry ("Outer temperature", "0",
                             Patterns::Double (),
                             "Temperature at the outer boundary (lithosphere water/air). Units: K.");
          prm.declare_entry ("Inner temperature", "6000",
                             Patterns::Double (),
                             "Temperature at the inner boundary (core mantle boundary) at the "
                             "beginning. Units: K.");
          prm.declare_entry ("Core density", "12.5e3",
                             Patterns::Double (),
                             "Density of the core. Units: kg/m^3");
          prm.declare_entry ("Gravity acceleration", "9.8",
                             Patterns::Double (),
                             "Gravitation acceleration at CMB. Units: m/s^2.");
          prm.declare_entry ("CMB pressure", "0.14e12",
                             Patterns::Double (),
                             "Pressure at CMB. Units: Pa.");
          prm.declare_entry ("Initial light composition", "0.01",
                             Patterns::Double (0),
                             "Initial light composition (eg. S,O) concentration "
                             "in weight fraction.");
          prm.declare_entry ("Max iteration", "30000",
                             Patterns::Integer (0),
                             "The max iterations for nonliner core energy solver.");
          prm.declare_entry ("Core heat capacity", "840",
                             Patterns::Double (0),
                             "Heat capacity of the core. Unit: J/kg/K");
          prm.declare_entry ("K0", "4.111e11",
                             Patterns::Double (0),
                             "Core compressibility at zero pressure. "
                             "Referring to Nimmo et al. (2004) for more detials.");
          prm.declare_entry ("Rho_0", "7.019e3",
                             Patterns::Double (0),
                             "Core density at zero pressure. Unit: kg/m^3"
                             "Referring to Nimmo et al. (2004) for more details.");
          prm.declare_entry ("Alpha", "1.35e-5",
                             Patterns::Double (0),
                             "Core thermal expansivity. Unit: 1/K");
          prm.declare_entry ("Lh", "750e3",
                             Patterns::Double (0),
                             "The latent heat of core freeze. Unit: J/kg");
          prm.declare_entry ("Beta_c", "1.1",
                             Patterns::Double (0),
                             "Compositional expansion coefficient."
                             "Referring to Nimmo et al. (2004) for more details.");
          prm.declare_entry ("Delta","0.5",
                             Patterns::Double (0,1),
                             "Partition coefficient of the light element.");
          prm.declare_entry ("k_c", "60",
                             Patterns::Double (0),
                             "Core heat conductivity. Unit: W/m/K");
          prm.enter_subsection("Geotherm parameters");
          {
              prm.declare_entry ("Tm0","1695",
                      Patterns::Double (0),
                      "Melting cure (Nimmo et al. [2004] eq. (40)) parameter Tm0. Unit: K");
              prm.declare_entry ("Tm1","10.9",
                      Patterns::Double (),
                      "Melting cure (Nimmo et al. [2004] eq. (40)) parameter Tm1. Unit: 1/Tpa");
              prm.declare_entry ("Tm2","-8.0",
                      Patterns::Double (),
                      "Melting cure (Nimmo et al. [2004] eq. (40)) parameter Tm2. Unit: 1/TPa^2");
              prm.declare_entry ("Theta","0.11",
                      Patterns::Double (),
                      "Melting cure (Nimmo et al. [2004] eq. (40)) parameter Theta.");
               prm.declare_entry ("Composition dependency","true",
                       Patterns::Bool (),
                       "If melting cure dependent on composition.");
              prm.declare_entry ("Ta1","3.5",
                      Patterns::Double (),
                      "Adiabatic temperature cure (Nimmo et al. [2004] eq. (41)) parameter Ta1. Unit: 1/TPa");
              prm.declare_entry ("Ta2","-1.8",
                      Patterns::Double (),
                      "Adiabatic temperature cure (Nimmo et al. [2004] eq. (41)) parameter Ta2. Unit: 1/TPa^2");
          }
          prm.leave_subsection ();
          prm.enter_subsection("Radioactive heat source");
          {
              prm.declare_entry ("Number of radioactive heating elements","0",
                      Patterns::Integer (0),
                      "Number of different radioactive heating elements in core");
              prm.declare_entry ("Heating rates","",
                      Patterns::List (Patterns::Double ()),
                      "Heating rates of different elements (W/kg)");
              prm.declare_entry ("Half life times","",
                      Patterns::List (Patterns::Double ()),
                      "Half decay times of different elements (Ga)");
              prm.declare_entry ("Initial concentrations","",
                      Patterns::List (Patterns::Double ()),
                      "Initial concentrations of different elements (ppm)");
          }
          prm.leave_subsection ();

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Dynamic_core<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Dynamic core");
        {
          inner_temperature = prm.get_double ("Inner temperature");
          outer_temperature = prm.get_double ("Outer temperature");
          Rho_cen           = prm.get_double ("Core density");
          g                 = prm.get_double ("Gravity acceleration");
          P_CMB             = prm.get_double ("CMB pressure");
          X_init            = prm.get_double ("Initial light composition");
          max_steps         = prm.get_integer ("Max iteration");
          Cp                = prm.get_double ("Core heat capacity");
          CpRho             = Cp*Rho_cen;

          //Nimmo et al. [2004]
          K0                = prm.get_double ("K0");
          Alpha             = prm.get_double ("Alpha");
          Rho_0             = prm.get_double ("Rho_0");
          Lh                = prm.get_double ("Lh");
          Beta_c            = prm.get_double ("Beta_c");
          k_c               = prm.get_double ("k_c");
          Delta             = prm.get_double ("Delta");
          
          prm.enter_subsection("Geotherm parameters");
          {
              Tm0           =  prm.get_double ("Tm0");
              Tm1           =  prm.get_double ("Tm1");
              Tm2           =  prm.get_double ("Tm2");
              Theta         =  prm.get_double ("Theta");
              Ta1           =  prm.get_double ("Ta1");
              Ta2           =  prm.get_double ("Ta2");
              composition_dependency 
                            =prm.get_bool("Composition dependency");
          }
          prm.leave_subsection ();
          
          prm.enter_subsection("Radioactive heat source");
          {
              n_radioheating_elements = prm.get_integer ("Number of radioactive heating elements");
              heating_rate = Utilities::string_to_double
                  (Utilities::split_string_list
                   (prm.get("Heating rates")));
              AssertThrow(n_radioheating_elements==heating_rate.size(),
                      ExcMessage("Number of heating rate entities does not match "
                          "the number of radioactive elements."));
              half_life = Utilities::string_to_double                          
                  (Utilities::split_string_list                                   
                   (prm.get("Half life times")));                                   
              AssertThrow(n_radioheating_elements==half_life.size(),           
                      ExcMessage("Number of half life time entities does not match "
                          "the number of radioactive elements."));                
              initial_concentration = Utilities::string_to_double                          
                  (Utilities::split_string_list                                   
                   (prm.get("Initial concentrations")));                                   
              AssertThrow(n_radioheating_elements==initial_concentration.size(),           
                      ExcMessage("Number of initial concentration entities does not match "
                          "the number of radioactive elements."));                            
          }
          prm.leave_subsection ();

          L=sqrt(3*K0*(log(Rho_cen/Rho_0)+1)/(2*M_PI*G*Rho_0*Rho_cen));
          D=sqrt(3*Cp/(2*M_PI*Alpha*Rho_cen*G));
          Mc=get_Mass(Rc);

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
    
    template <int dim>
    Dynamic_core<dim>::Dynamic_core():
    //Parameters
    G(6.67221937e-11)
    {
        is_first_call=true;
    }

    template <int dim>
    double
    Dynamic_core<dim>::get_Ri(double X,double T,double R) const
    {
        
        double const1,const2,a,b,c,p_ic1,p_ic2,P_ic;
        if(composition_dependency)
            const1 = Tm0*(1 - Theta*X);
        else
            const1 = Tm0*(1 - Theta);
        const2 = T/(1+ Ta1*P_CMB + Ta2*(P_CMB*P_CMB));
        a = const1*Tm2 - const2*Ta2;
        b = const1*Tm1 - const2*Ta1;
        c = const1 - const2;
        p_ic1=  (-b + sqrt(b*b - 4*a*c))/(2*a);
        p_ic2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
        P_ic = p_ic1;
        if(P_ic<P_Core)
            return (get_r_from_p(P_ic,R));
        else
            return 0.;
    }


    template <int dim>
    double
    Dynamic_core<dim>::get_a_from_R(double R)
    {
        return( 1. - get_Qlg(core_data.Ti,R)/(core_data.Q-core_data.Qr)/core_data.dt );
    }
            
    template <int dim>
    void
    Dynamic_core<dim>::solve_time_step(double &X, double &T, double &R)
    {
        // Well solving the change in core-mantle boundary temperature T, inner core radius R, and 
        //    light component (e.g. S, O, Si) composition X, the following relations has to be respected:
        // 1. At the inner core boundary the adiabatic temperature should be equal to solidu temperature
        // 2. The following energy production rate should be balanced in core:
        //    Heat flux at core-mantle boundary         Q
        //    Specific heat                             Qs*dT/dt
        //    Radioactive heating                       Qr
        //    Gravitational contribution                Qg*dR/dt
        //    Latent heat                               Ql*dR/dt
        //    So that         Q+Qs*dT/dt+Qr*H+Qg*dR/dt*Ql*dR/dt=0
        // 3. The light component composition X depends on inner core radius (See function get_X() ),
        //    and core solidus may dependent on X as well
        // This becomes a small nonliner problem. Directly iterate through the above three system doesn't 
        //   converge well. Alternatively, we define the energy distribution between core cooling and inner 
        //   core growth as a fraction factor 'a'. By solving factor 'a' using bisection method, we get easy 
        //   converge.

        int steps=1;
        double a0,a1,a2;
        double T_0,T_1,T_2;
        double R_0,R_1,R_2;
        double X_0,X_1,X_2;
        double dT0,dT1,dT2;
        a0=0;
        dT0=get_dT(X_0,T_0,R_0,a0);
        a2=1;
        dT2=get_dT(X_2,T_2,R_2,a2);

        while(!(dT0==0 || dT2==0 || steps>max_steps))
        {
            // If solution is out of the interval, then something is wrong. 
            if(dT0*dT2>0)
            {
              cout<<"Step: "<<steps<<endl
                <<" X=["<<X_0<<","<<X_2<<"]"
                <<" T=["<<T_0<<","<<T_2<<"]"<<"(K)"
                <<" R=["<<R_0/1e3<<","<<R_2/1e3<<"]"<<"(km)"
                <<endl;
              cout<<core_data.Q<<endl;
              AssertThrow(dT0*dT2<=0,ExcMessage("No single solution for inner core!"));
            }
            // Get the middle point of the interval
            a1=(a0+a2)/2;
            dT1=get_dT(X_1,T_1,R_1,a1);
            // If the returning R_0 is large than Rc the core has to be complete frozen
            // in this time step. Get factor 'a' base on that.
            if(R_0>Rc)
            {
                a0=get_a_from_R(Rc);
                dT0=get_dT(X_0,T_0,R_0,a0);
                dT0=0;
            }
            // If the returning R_0 is negative the core has to be complete molten 
            // in this time step. Get factor 'a' base on that       
            if(R_0<0)
            {
                a0=get_a_from_R(0.);
                dT0=get_dT(X_0,T_0,R_0,a0);
                dT0=0;
            }
            // If the solution is in left interval
            if(dT1*dT0<=0)
            {
                a2=a1;
                dT2=dT1;
                T_2=T_1;
                R_2=R_1;
                X_2=X_1;
            }
            // If the solution is in right interval
            else if(dT1*dT2<=0)
            {
                a0=a1;
                dT0=dT1;
                T_0=T_1;
                R_0=R_1;
                X_0=X_1;
            }
            steps++;
        }
        // Choose which side of the interval would be the solution
        if(dT2==0)
        {
            T=T_2;
            R=R_2;
            X=X_2;
        }
        else
        {
            T=T_0;
            R=R_0;
            X=X_0;
        }
    }
    
    template <int dim>
    double
    Dynamic_core<dim>::get_dT(double &X, double &T, double &R,double a)
    {
        double p,Ta,Ts;
        T=core_data.Ti+(core_data.Q+core_data.Qr)*a/core_data.Qs*core_data.dt;
        R=get_R_loop(T,(1-a)*(core_data.Q+core_data.Qr)*core_data.dt);
        X=get_X(R);
        p=get_Pressure(R);
        Ta=get_adiabatic(T,p);
        Ts=get_solidus(X,p);
        if(R==0 && Ts<Ta)
            return 0.;
        else
            return(Ta-Ts);
    }

    template <int dim>
    double
    Dynamic_core<dim>::get_R_loop(double T, double Q)
    {
        double R0,R1,R2;
        double Q0,Q1,Q2;
        if(Q==0.)
            return core_data.Ri;
        else
        {
            R0=0.;
            R2=Rc;
            Q0=get_Qlg(T,R0);
            Q2=get_Qlg(T,R2);
            if(Q2<Q)
                return(Rc+1.e3);//Even core complete freeze is not enough for Q, pass R>Rc
            if(Q<Q0)
                return(-1.e3);//Even core complete molten is not enough for Q, pass R<0
            double steps=1;
            while(!( (Q0-Q)*(Q2-Q)==0 || steps>max_steps) )
            {
                // If no solution in the interval, something is wrong.
                AssertThrow((Q0-Q)*(Q2-Q)<=0,ExcMessage("No single solution for core radius!"));
                R1=(R0+R2)/2.;
                Q1=get_Qlg(T,R1);
                if(Q1==Q)
                    return(R1);
                else if(Q1>Q)
                {
                    R2=R1;
                    Q2=Q1;
                }
                else
                {
                    R0=R1;
                    Q0=Q1;
                }
                steps++;
            }
            if(Q0==Q)
                return(R0);
            else if(Q2==Q)
                return(R2);
            else
                return((R0+R2)/2.);
        }
    }
    template <int dim>
    double
    Dynamic_core<dim>::get_Qlg(double T, double R)
    {
        double Qg,Eg,Ql,El;
        // Using the mid-point of current inner core radius and
        // radius of last time step as the radius for gravity and 
        // latent heat factor calculations.
        double Ri=(core_data.Ri+R)/2.;
        double X=get_X(Ri);
        get_gravity_heating(T,Ri,X,Qg,Eg);
        get_latent_heating(T,Ri,El,Ql);
        return((Ql+Qg)*(R-core_data.Ri));
    }

    template <int dim>
    void
    Dynamic_core<dim>::update_core_data()
    {
        get_specific_heating(core_data.Ti,core_data.Qs,core_data.Es);
        get_radio_heating(core_data.Ti,core_data.Qr,core_data.Er);
        get_gravity_heating(core_data.Ti,core_data.Ri,core_data.Xi,core_data.Qg,core_data.Eg);
        get_adiabatic_heating(core_data.Ti,core_data.Ek,core_data.Qk);
        get_latent_heating(core_data.Ti,core_data.Ri,core_data.El,core_data.Ql);
    }

    template <int dim>
    const struct _Core_Data*
    Dynamic_core<dim>::get_core_data() const
    {
        return &core_data;
    }

    template <int dim>
    double
    Dynamic_core<dim>::get_adiabatic(double Tc, double p) const
    {
        return( Tc * (1+Ta1*p+Ta2*pow(p,2)) / (1+Ta1*P_CMB+Ta2*pow(P_CMB,2)) );
    }

    template <int dim>
    double
    Dynamic_core<dim>::get_solidus(double X,double p) const
    {
        if(composition_dependency)
            return(Tm0*(1-Theta*X)*(1+Tm1*p+Tm2*pow(p,2)));
        else
            return(Tm0*(1-Theta)*(1+Tm1*p+Tm2*pow(p,2)));
    }

    template <int dim>
    double
    Dynamic_core<dim>::get_X(double r) const
    {
        double xi_3=pow(r/Rc,3);
        return X_init/(1-xi_3+Delta*xi_3);
    }

    template <int dim>
    void
    Dynamic_core<dim>::update()
    {
        const Postprocess::DynamicCoreStatistics<dim> * dynamic_core_statistics
          = this->template find_postprocessor<const Postprocess::DynamicCoreStatistics<dim> >();
        AssertThrow(dynamic_core_statistics!=NULL,
            ExcMessage ("Dynamic core boundary condition has to work with dynamic core statistics postprocessor."));
        core_data.Q=dynamic_core_statistics->get_CMB_heat_flux();
        core_data.dt=this->get_timestep();
        core_data.H=get_radioheating_rate();
        if(is_first_call==true)
        {
            const GeometryModel::SphericalShell<dim>* shperical_shell_geometry =
                dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&(this->get_geometry_model()));
            AssertThrow (shperical_shell_geometry != NULL,
                    ExcMessage ("This boundary model is only implemented if the geometry is "
                                "in fact a spherical shell."));
            Rc=shperical_shell_geometry->R0;
            P_Core=get_Pressure(0);

            core_data.Ti=inner_temperature;
            core_data.Ri=get_Ri(X_init,inner_temperature,0);
            core_data.Xi=get_X(core_data.Ri);
            core_data.Q=0.;
            core_data.dt=0.;
            core_data.dT_dt=0.;
            core_data.dR_dt=0.;
            core_data.dX_dt=0.;

            is_first_call=false;

        }
        else if(core_data.Q!=0.)
        {
            double X1,R1,T1;
            solve_time_step(X1,T1,R1);
            if(core_data.dt!=0)
            {
                core_data.dR_dt=(R1-core_data.Ri)/core_data.dt;
                core_data.dT_dt=(T1-core_data.Ti)/core_data.dt;
                core_data.dX_dt=(X1-core_data.Xi)/core_data.dt;
            }
            else
            {
                core_data.dR_dt=0.;
                core_data.dT_dt=0.;
                core_data.dX_dt=0.;
            }
            core_data.Xi=X1;
            core_data.Ri=R1;
            core_data.Ti=T1;
            inner_temperature=T1;
        }
        update_core_data();
    }
    template <int dim>
    double
    Dynamic_core<dim>::
    get_Mass(double r) const
    {
        return 4.*M_PI*Rho_cen*(-pow(L,2)/2.*r*exp(-pow(r/L,2))+pow(L,3)/4.*sqrt(M_PI)*erf(r/L));
    }

    template <int dim>
    double
    Dynamic_core<dim>::
    fun_Sn(double B,double R,double n) const
    {
        double S=R/(2.*sqrt(M_PI));
        for(unsigned i=1;i<=n;i++)
            S+=B/sqrt(M_PI)*exp(-pow(i,2)/4.)/(i*sinh(i*R/B));
        return S;
    }
    template <int dim>
    double
    Dynamic_core<dim>::
    get_Pressure(double r) const
    {
        return P_CMB-(4*M_PI*G*pow(Rho_cen,2))/3*((3*pow(r,2)/10.-pow(L,2)/5)*exp(-pow(r/L,2))-(3*pow(Rc,2)/10-pow(L,2)/5)*exp(-pow(Rc/L,2)));
    }
    
    template <int dim>
    double
    Dynamic_core<dim>::
    get_dp_dr(double r) const
    {
        return -(5-3*pow(r/L,2))*r/5*exp(-pow(r/L,2));
    }
    template <int dim>
    double
    Dynamic_core<dim>::
    get_r_from_p(double p,double r0) const
    {
        double r_ini=r0;
        if(p>get_Pressure(0))
            return 0.;
        else if(p<P_CMB)
            return Rc;
        else
        {
            double r_error=1e-3;//R_error*1e-2;
            double i_max=1000;//Max_steps;
            double dr=r_error*10;
            for(unsigned i=0;i<i_max;i++)
            {
                if(r0<=0)r0=r_error;
                dr=(p-get_Pressure(r0))/get_dp_dr(r0);
                r0+=dr;
                if(r0>Rc)r0=Rc;
                if(fabs(dr)<r_error)
                    return r0;
            }
            cout<<"p="<<p<<",r0="<<r_ini<<",r="<<r0<<endl;
            AssertThrow(fabs(dr)<r_error,ExcMessage("Get r from pressure not converge!"));
        }
        return 0.;
    }

    template <int dim>
    double
    Dynamic_core<dim>::
    get_Rho(double r) const
    {
        return Rho_cen*exp(-pow(r/L,2));
    }

    template <int dim>
    double
    Dynamic_core<dim>::
    get_g(double r) const
    {
        return (4*M_PI/3)*G*Rho_cen*r*(1-3*pow(r,2)/(5*pow(L,2)));
    }

    template <int dim>
    double
    Dynamic_core<dim>::
    get_T(double Tc,double r) const
    {
        return Tc*exp((pow(Rc,2)-pow(r,2))/pow(D,2));
    }

    template <int dim>
    double
    Dynamic_core<dim>::
    get_gravity_potential(double r) const
    {
        return 2./3.*M_PI*G*Rho_cen*(pow(r,2)*(1.-3.*pow(r,2)/(10.*pow(L,2)))-pow(Rc,2)*(1.-3.*pow(Rc,2)/(10.*pow(L,2))));
    }

    template <int dim>
    void
    Dynamic_core<dim>::
    get_specific_heating(double Tc, double &Qs,double &Es)
    {
        double A=sqrt(1./(pow(L,-2)+pow(D,-2)));
        double Is=4.*M_PI*get_T(Tc,0.)*Rho_cen*(-pow(A,2)*Rc/2.*exp(-pow(Rc/A,2))+pow(A,3)*sqrt(M_PI)/4.*erf(Rc/A));

        Qs=-Cp/Tc*Is;
        Es=Cp/Tc*(Mc-Is/Tc);
    }

    template <int dim>
    void
    Dynamic_core<dim>::
    get_radio_heating(double Tc, double &Qr, double &Er)
    {
        double B,It;
        if(D>L)
        {
            B=sqrt(1/(1/pow(L,2)-1/pow(D,2)));
            It=4*M_PI*Rho_cen/get_T(Tc,0)*(-pow(B,2)*Rc/2*exp(-pow(Rc/B,2))+pow(B,3)/sqrt(M_PI)/4*erf(Rc/B));
        }
        else
        {
            B=sqrt(1/(pow(D,-2)-pow(L,-2)));
            It=4*M_PI*Rho_cen/get_T(Tc,0)*(pow(B,2)*Rc/2*exp(pow(Rc/B,2))-pow(B,2)*fun_Sn(B,Rc,100)/2);
        }
        Qr=Mc*core_data.H;
        Er=(Mc/Tc-It)*core_data.H;

    }
    
    template <int dim>
    void
    Dynamic_core<dim>::
    get_gravity_heating(double Tc, double r,double X,double &Qg,double &Eg)
    {
        double Cc=4*M_PI*pow(r,2)*get_Rho(r)*X/(Mc-get_Mass(r));
        double C_2=3./16.*pow(L,2)-0.5*pow(Rc,2)*(1.-3./10.*pow(Rc/L,2));
        Qg=(8./3.*pow(M_PI*Rho_cen,2)*G*(
                    ((3./20.*pow(Rc,5)-pow(L,2)*pow(Rc,3)/8.-C_2*pow(L,2)*Rc)*exp(-pow(Rc/L,2))
                       +C_2/2.*pow(L,3)*sqrt(M_PI)*erf(Rc/L))
                   -((3./20.*pow(r,5)-pow(L,2)*pow(r,3)/8.-C_2*pow(L,2)*r)*exp(-pow(r/L,2))
                       +C_2/2.*pow(L,3)*sqrt(M_PI)*erf(r/L)))
                -(Mc-get_Mass(r))*get_gravity_potential(r))*Beta_c*Cc;
        Eg=Qg/Tc;

    }

    template <int dim>
    void
    Dynamic_core<dim>::
    get_adiabatic_heating(double Tc, double &Ek, double &Qk)
    {
        Ek=16*M_PI*k_c*pow(Rc,5)/5/pow(D,4);
        Qk=8*M_PI*pow(Rc,3)*k_c*Tc/pow(D,2);
    }
    template <int dim>
    void
    Dynamic_core<dim>::
    get_latent_heating(double Tc, double r, double &El, double &Ql)
    {
        Ql=4.*M_PI*pow(r,2)*Lh*get_Rho(r);
        El=Ql*(get_T(Tc,r)-Tc)/(Tc*get_T(Tc,r));
    }
  
    template <int dim>
    double
    Dynamic_core<dim>::
    get_radioheating_rate() const
    {
      double time=this->get_time()+0.5*this->get_timestep();
      double Ht=0;
      for(unsigned i=0;i<n_radioheating_elements;i++)
        Ht+=heating_rate[i]*initial_concentration[i]*1e-6*pow(0.5,time/half_life[i]/year_in_seconds/1e9);
      return Ht;
    }

  }

}


// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(Dynamic_core,
                                               "Dynamic core",
                                               "A model in which the initial temperature is chosen on "
                                               "the inner and outer boundaries of a spherical shell "
                                               "and CMB temperature is changing with core energy balance."
                                               "Parameters are read from subsection 'Dynamic core'.")
  }
}
