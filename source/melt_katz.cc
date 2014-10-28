#include "aspect/melt_katz.h"
#include "aspect/bisection.h"
#include<cmath>

namespace Melt_Katz
{
    double
    Melt_Katz::get_modified_temperature(double T, double P, double Mcpx, double X, double &melt)
    {
        double delta_T=melt_fraction.get_deltaT(melt_fraction.get_X(X,0.),P);
        double Ts=melt_fraction.T_solidus(P)-delta_T;
        this->T=T-273.15;
        this->P=P;
        this->Mcpx=Mcpx;
        this->X=X;

        //double return_T=this->solve();
        double return_T=Bisection::bisecion_solve(*this,Ts,this->T,1e-3);
        melt=melt_fraction.get_melt_fraction(return_T,P,Mcpx,X);
        return return_T;
    }

    double
    Melt_Katz::Melt_fraction::get_deltaT(double XH2O, double P) const
    {
        if(XH2O==0)
            return 0.;
        else
        {
            return parameters.K*pow(std::min(get_Xs(P),XH2O),parameters.gamma);
        }
    }

    double
    Melt_Katz::Melt_fraction::T_solidus(double P) const
    {
        return parameters.A1+parameters.A2*(P*1e-9)+parameters.A3*pow(P*1e-9,2);
    }

    double
    Melt_Katz::Melt_fraction::T_lherz(double P) const
    {
        return parameters.B1+parameters.B2*(P*1e-9)+parameters.B3*pow(P*1e-9,2);
    }

    double
    Melt_Katz::Melt_fraction::T_liquidus(double P) const
    {
        return parameters.C1+parameters.C2*(P*1e-9)+parameters.C3*pow(P*1e-9,2);
    }

    double
    Melt_Katz::Melt_fraction::get_Xs(double P) const
    {
        return parameters.X1*pow(P*1e-9,parameters.lambda)+parameters.X2*(P*1e-9);
    }

    double
    Melt_Katz::Melt_fraction::get_X(double X_bulk,double melt_fraction) const
    {
        return X_bulk/(parameters.D_H2O+melt_fraction*(1.-parameters.D_H2O));
    }

    double
    Melt_Katz::Melt_fraction::melt_fraction_x_melt(double T, double P, double Mcpx, double X_H2O_melt) const
    {
        double delta_T=get_deltaT(X_H2O_melt,P);
        double Ts=T_solidus(P)-delta_T;
        double Tl=T_liquidus(P)-delta_T;
        double Tll=T_lherz(P)-delta_T;
        double Fcpx=pow(std::max(T-Ts,0.)/(Tll-Ts),parameters.beta1);
        double Fcpx_out=Mcpx/(parameters.r1+parameters.r2*P*1e-9);
        if(Fcpx>Fcpx_out)
            // cpx depleted melt
        {
            double Tcpx_out=pow(Fcpx_out,1./parameters.beta1)*(Tll-Ts)+Ts;
            return Fcpx_out+(1.-Fcpx_out)*pow((T-Tcpx_out)/(Tl-Tcpx_out),parameters.beta2);
        }
        else if(X_H2O_melt==0)
            // Anhydrous Melting
        {
            return Fcpx;
        }
        else
            // Hydrous Melting
        {
            double dTs=T-Ts;
            if(dTs<0.)
                return 0.;
            else
                return pow((T-Ts)/(Tll-Ts),parameters.beta1);
        }
    }

    double
    Melt_Katz::Melt_fraction::get_melt_fraction(double T, double P, double Mcpx, double X_H2O)

    {
        this->T=T;
        this->P=P;
        this->Mcpx=Mcpx;
        this->X_H2O=X_H2O;
        return Bisection::bisecion_solve(*this,0,1,1e-6);
    }
}
