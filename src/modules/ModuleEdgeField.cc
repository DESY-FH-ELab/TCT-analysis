/**
 * \file
 * \brief Implementation of TCT::ModuleEdgeField methods.
 */

#include "modules/ModuleEdgeField.h"
#include "TStyle.h"

#ifdef USE_GUI
#include "QVBoxLayout"
#include "QLabel"
#endif

namespace TCT {

/// Analyse data (should be reimplemented by developer)
bool ModuleEdgeField::Analysis() {

    //creating of the folders Velocity, Velocity_Normed and Velocity_Diode
    //containing velocity profiles not normed (through solving numerical
    //equation with bias voltage), normed with photodiode data and
    //calculated from physical assumptions(not finished yet)
    TDirectory *dir_vel = f_rootfile->mkdir("Velocity");
    TDirectory *dir_vel_normed;
    if(config->CH_PhDiode()) dir_vel_normed = f_rootfile->mkdir("Velocity_Normed");
    TDirectory *dir_vel_diode;
    if(config->CH_PhDiode()) dir_vel_diode = f_rootfile->mkdir("Velocity_Diode");
    dir_vel->cd();

    //memory allocation and setting of variables corresponding to IDs of the axis
    Int_t numVolt,numS;
    Float_t Ss,Sc0;

    Int_t ChNumber=(config->CH_Det())-1;              // select the oscilloscope channel
    Int_t volt_source=(config->VoltSource());            // select
    Int_t scanning_axis=config->ScAxis()-1;         // select scanning axis (0=x,1=y,2=z)
    Float_t *voltages;

    switch(volt_source)
      {
      case 1: numVolt=stct->NU1; voltages=stct->U1.fArray; break; //u1
      case 2: numVolt=stct->NU2; voltages=stct->U2.fArray; break; //u2
      }

    //memory allocation for charge, velocity and field profiles
    TGraph **cc = new TGraph*[numVolt];
    TGraph **cc_norm; // charge collection graph
    TGraph **ph_charge = new TGraph*[numVolt];
    TGraph **field = new TGraph*[numVolt];
    TGraph **velocity_electrons = new TGraph*[numVolt];
    TGraph **velocity_holes = new TGraph*[numVolt];

    TGraph **field_normed = new TGraph*[numVolt];
    TGraph **velocity_electrons_normed = new TGraph*[numVolt];
    TGraph **velocity_holes_normed = new TGraph*[numVolt];

    TGraph **velocity_diode = new TGraph*[numVolt];
    TGraph **field_diode = new TGraph*[numVolt];

    // setting the scanning axis
    SwitchAxis(scanning_axis,numS,Ss,Sc0);


    //calculate charge profiles for different voltages
    CalculateCharges(ChNumber,volt_source+2,numVolt,scanning_axis,numS,cc,config->FTlow(),config->FTlow()+GetEV_Time());

    //find integration ranges (for the maximal voltage in scan, hope that sensor is fully depleted
    //otherwise may be a problem with right edge
    TCTWaveform *wf0;
    switch(volt_source)
      {
      case 1: wf0 = stct->Projection(ChNumber,scanning_axis,0,0,0,numVolt-1,0,numS); break;
      case 2: wf0 = stct->Projection(ChNumber,scanning_axis,0,0,0,0,numVolt-1,numS); break;
      }
    TGraph *charge_max_bias = wf0->CCE(config->FTlow(),config->FThigh());
    delete wf0;
    Double_t left_edge,right_edge;
    FindEdges(charge_max_bias,numS,Ss,left_edge,right_edge);
    std::cout<<std::endl;
    std::cout<<"\tFound detector thickness is: "<<right_edge-left_edge<<" micrometers"<<std::endl;
    std::cout<<"\tIntegrating in the range "<<left_edge<<" - "<<right_edge<<std::endl;
    std::cout<<std::endl;


    //looking for integration ranges (Int)
    Double_t *temp_max_bias_charge = charge_max_bias->GetX();
    Int_t ix1=-1;
    Int_t ix2=-1;
    for(Int_t i=0;i<numS;i++) {
        if(temp_max_bias_charge[i]<=left_edge) ix1=i;
        if(temp_max_bias_charge[i]<=right_edge) ix2=i;
    }
    if(ix1==-1) ix1=0;
    if(ix2==-1) ix2=numS-1;
    delete temp_max_bias_charge;

    //calculate laser charge profiles
    if(config->CH_PhDiode()) CalculateCharges(config->CH_PhDiode()-1,volt_source+2,numVolt,scanning_axis,numS,ph_charge,config->FDLow(),config->FDHigh());

    //calculating the velocity profile asuuming integral(E)dx = Vbias
    Double_t eps = 1e-3;
    Double_t *temp_sensor;
    Double_t *temp_vel_h=new Double_t[numS];
    Double_t *temp_vel_el=new Double_t[numS];
    Double_t *temp_field=new Double_t[numS];
    for(int j=0;j<numVolt;j++) { // FIXME constants seems to prefer smaller values towards low electric fields , why?
        temp_sensor = cc[j]->GetY();

        Double_t a=1;
        if(GraphIntegral(cc[j],left_edge,right_edge)<0) a = -a;
        Double_t dx = 0.5*a;
        Double_t sum = 0;
        while(abs(voltages[j]-sum)>eps) {

            for(int i=0; i<numS; i++) {
                temp_field[i] = BiSectionMethod(eps,-5e3,1e6,temp_sensor[i],a);
            }
            sum=0;
            for(int i=ix1;i<=ix2;i++) sum+=temp_field[i];
            sum*=1e-4*Ss; // conversion from um to cm
            if(voltages[j]-sum>0) a-=dx;
            else a+=dx;
            dx=dx/2;
        }
        std::cout<<"U = "<<voltages[j]<<" norm const: "<<a<<std::endl;
        for(int i=0;i<numS;i++) {
            temp_field[i] = 1e-4*temp_field[i];
            temp_vel_h[i] = 1e4*temp_field[i]*Mu(temp_field[i],1)/(1+Mu(temp_field[i],1)*1e4*temp_field[i]/config->v_sat());
            temp_vel_el[i] = 1e4*temp_field[i]*Mu(temp_field[i],0)/(1+Mu(temp_field[i],0)*1e4*temp_field[i]/config->v_sat());
        }
        // building of the graphs
        field[j] = GraphBuilder(numS,cc[0]->GetX(),temp_field,"scanning distance [#mum]", "Electric Field, [V/#mum]","Electric Field");
        velocity_holes[j] = GraphBuilder(numS,cc[0]->GetX(),temp_vel_h,"scanning distance [#mum]", "Velocity [cm/s]","Holes Velocity Profile");
        velocity_electrons[j] = GraphBuilder(numS,cc[0]->GetX(),temp_vel_el,"scanning distance [#mum]", "Velocity [cm/s]","Electrons Velocity Profile");

    }

    //calculating the velocity profile asuuming integral(E)dx = Vbias with normed data
    Float_t *normcoeff = new Float_t[numVolt];
    if(config->CH_PhDiode()) {
        //calculating the mean charge from the photodetector
        cc_norm=NormedCharge(cc,ph_charge,numVolt);

        Double_t *temp_sensor;

        for(int j=0;j<numVolt;j++) {
            temp_sensor = cc_norm[j]->GetY();
            Double_t a=1;
            if(GraphIntegral(cc_norm[j],left_edge,right_edge)<0) a = -a;
            Double_t dx = 0.5*a;
            Double_t sum = 0;
            while(abs(voltages[j]-sum)>eps) {

                for(int i=0; i<numS; i++) {
                    temp_field[i] = BiSectionMethod(eps,-5e3,1e6,temp_sensor[i],a);
                }
                sum=0;
                for(int i=ix1;i<=ix2;i++) sum+=temp_field[i];
                sum*=1e-4*Ss;
                if(voltages[j]-sum>0) a-=dx;
                else a+=dx;
                dx=dx/2;
            }
            if(GraphIntegral(cc_norm[j],left_edge,right_edge)<0) normcoeff[j]=-a;
            else normcoeff[j]=a;
            //std::cout<<"U = "<<voltages[j]<<" norm const: "<<a<<std::endl;
            for(int i=0;i<numS;i++) {
                temp_field[i] = 1e-4*temp_field[i];
                temp_vel_h[i] = 1e4*temp_field[i]*Mu(temp_field[i],1)/(1+Mu(temp_field[i],1)*1e4*temp_field[i]/config->v_sat());
                temp_vel_el[i] = 1e4*temp_field[i]*Mu(temp_field[i],0)/(1+Mu(temp_field[i],0)*1e4*temp_field[i]/config->v_sat());
            }
            field_normed[j] = GraphBuilder(numS,cc[0]->GetX(),temp_field,"scanning distance [#mum]", "Electric Field, [V/#mum]","Electric Field");
            velocity_holes_normed[j] = GraphBuilder(numS,cc[0]->GetX(),temp_vel_h,"scanning distance [#mum]", "Velocity [cm/s]","Holes Velocity Profile");
            velocity_electrons_normed[j] = GraphBuilder(numS,cc[0]->GetX(),temp_vel_el,"scanning distance [#mum]", "Velocity [cm/s]","Electrons Velocity Profile");

        }
    }

    delete temp_vel_h;
    delete temp_vel_el;
    delete temp_field;


    //calculating the velocity profile using photodiode for estimate of N_e,h, Factor 100 to many :(
    if(config->CH_PhDiode()) {

        Double_t temp_sr_general = 0;
        Double_t *temp_diode;
        Double_t temp_sr = 0;
        for(int j=0;j<numVolt;j++) {
            temp_diode = ph_charge[j]->GetY();
            temp_sr = 0;
            for(int i=0;i<numS;i++) {
                temp_sr += temp_diode[i];
            }
            temp_sr = temp_sr/(numS);
            temp_sr_general+=temp_sr;
        }
        temp_sr_general=temp_sr_general/(numVolt);

        Double_t ampl = config->ampl(); // FIXME needs accurate update!
        Double_t Res_sensor = config->R_sensor();
        Double_t Res_photo = config->R_diode();
        Double_t Response_photo = config->RespPhoto();
        Double_t Eweight = 1./(config->SampleThickness()*1e-4);
        Double_t E_pair = config->E_pair();
        Double_t diode_multi = config->light_split();

        Double_t *temp_sensor_charge;
        Double_t *temp_diode_charge;
        Double_t *temp_velocity_avg = new Double_t[numS];
        Double_t *temp_field = new Double_t[numS];

        Double_t Neh = 0.624*diode_multi*temp_sr_general/Res_photo/Response_photo/E_pair;
        Neh = Neh*0.02146; // Integral over one strip from energy depostion a.f.o. depth for 80 um strip, and lambda = 11 1/cm
        std::cout<<"Photodiode charge: "<<temp_sr_general<<std::endl;
        std::cout<<"Npairs: "<<1.e7*Neh<<std::endl;

        for(int j=0;j<numVolt;j++) {
            temp_sensor_charge = cc[j]->GetY();
            temp_diode_charge = ph_charge[j]->GetY();
            for(int i=0;i<numS;i++) {
                temp_velocity_avg[i] = -1e9*0.624*temp_sr_general*temp_sensor_charge[i]/(GetEV_Time()*temp_diode_charge[i]*Eweight*ampl*Res_sensor*Neh);
                temp_field[i] = 1e-4*temp_velocity_avg[i]/(config->mu0_els()+config->mu0_holes());
            }
            velocity_diode[j] = GraphBuilder(numS,cc[0]->GetX(),temp_velocity_avg,"scanning distance [#mum]", "Velocity [cm/s]","Velocity Profile");
            field_diode[j] = GraphBuilder(numS,cc[0]->GetX(),temp_field,"scanning distance [#mum]", "Electric Field, [V/#mum]","Electric Field Profile");
        }

        delete temp_velocity_avg;
        delete temp_field;

    }

    //write separate charges
    if(config->FSeparateCharges()) {
        dir_vel->cd();
        GraphSeparate(numVolt,velocity_holes,"velocity_profiles_holes","scanning distance, [#mum]","Velocity [cm/s]","Holes Velocity Profile","U = ",voltages);
        GraphSeparate(numVolt,velocity_electrons,"velocity_profiles_electrons","scanning distance, [#mum]","Velocity [cm/s]","Electrons Velocity Profile","U = ",voltages);
        GraphSeparate(numVolt,field,"field_profiles","scanning distance, [#mum]","Electric Field, [V/#mum]","Electric Field Profile","U = ",voltages);
        if(config->CH_PhDiode()) {
            dir_vel_normed->cd();
            GraphSeparate(numVolt,velocity_holes_normed,"velocity_profiles_holes","scanning distance, [#mum]","Velocity [cm/s]","Holes Velocity Profile","U = ",voltages);
            GraphSeparate(numVolt,velocity_electrons_normed,"velocity_profiles_electrons","scanning distance, [#mum]","Velocity [cm/s]","Electrons Velocity Profile","U = ",voltages);
            GraphSeparate(numVolt,field_normed,"field_profiles","scanning distance, [#mum]","Electric Field, [V/#mum]","Electric Field Profile","U = ",voltages);
            dir_vel_diode->cd();
            GraphSeparate(numVolt,ph_charge,"photodiode","scanning distance, [#mum]","Charge,[arb.]","Laser Charge vs Distance","U = ",voltages);
            GraphSeparate(numVolt,velocity_diode,"velocity_profiles","scanning distance, [#mum]","Velocity [cm/s]","Sum Velocity Profile","U = ",voltages);
            GraphSeparate(numVolt,field_diode,"field_profiles","scanning distance, [#mum]","Electric Field, [V/#mum]","Electric Field Profile","U = ",voltages);
            dir_vel->cd();
        }
    }

    //plot graphs for different voltages
    MultiGraphWriter(numVolt,field,"scanning distance [#mum]","Electric Field, [V/#mum]","Electric Field","Electric_Field");
    MultiGraphWriter(numVolt,velocity_holes,"scanning distance [#mum]","Velocity, [cm/s]","Holes Velocity Profile","Vel_Holes");
    MultiGraphWriter(numVolt,velocity_electrons,"scanning distance [#mum]","Velocity, [cm/s]","Electrons Velocity Profile","Vel_Electrons");

    //plot graphs for different voltages with normed data
    if(config->CH_PhDiode()) {
        dir_vel_normed->cd();
        MultiGraphWriter(numVolt,field_normed,"scanning distance [#mum]","Electric Field, [V/#mum]","Electric Field","Electric_Field");
        MultiGraphWriter(numVolt,velocity_holes_normed,"scanning distance [#mum]","Velocity, [cm/s]","Holes Velocity Profile","Vel_Holes");
        MultiGraphWriter(numVolt,velocity_electrons_normed,"scanning distance [#mum]","Velocity, [cm/s]","Electrons Velocity Profile","Vel_Electrons");

        //plots normalization coefficient
        TGraph *coeff = GraphBuilder(numVolt,voltages,normcoeff,"Voltage, [V]","A","Norm Coefficient for Different voltages");
        TF1 *fff = new TF1("log0","[0]*log(x)",10,120);
        coeff->Fit("log0","R");
        coeff->Write("CoeffNorm");
        delete coeff;
        delete fff;

        dir_vel_diode->cd();
        MultiGraphWriter(numVolt,velocity_diode,"scanning distance [#mum]","Velocity [cm/s]","Velocity Profiles","VelocityVsDist_diode");
        MultiGraphWriter(numVolt,field_diode,"scanning distance [#mum]","Electric Field, [V/#mum]","Electric Field","FieldVsDist_diode");
        dir_vel->cd();
    }

    //freeing the memory
    delete normcoeff;

    delete velocity_holes;
    delete velocity_electrons;
    delete velocity_electrons_normed;
    delete velocity_holes_normed;
    delete field_normed;
    delete field_diode;
    delete velocity_diode;
    delete field;

    delete cc;
    if(config->CH_PhDiode()) delete cc_norm;
    delete ph_charge;

    return true;
}

/// Check before analysis (should be reimplemented by developer)
bool ModuleEdgeField::CheckModuleData() {

    std::cout<<"\t- DO "<<GetTitle()<<"(type - "<<GetType()<<")\n";

    Int_t NSource;
    Int_t NSc;
    std::cout<<"\t\t-- Voltage Source: "<< config->VoltSource() <<std::endl;
    switch(config->VoltSource()) {
    case 1: NSource = stct->NU1; break;
    case 2: NSource = stct->NU2; break;
    }
    if(NSource>=1) std::cout<<"\t\t\tVoltage scan contains "<<NSource<<" points. OK"<<std::endl;
    std::cout<<"\t\t-- Scanning Axis: "<< config->ScAxis() <<std::endl;
    switch(config->ScAxis()) {
    case 1: NSc = stct->Nx; break;
    case 2: NSc = stct->Ny; break;
    case 3: NSc = stct->Nz; break;
    }
    if(NSc>10) std::cout<<"\t\t\tScanning axis contains "<<NSc<<" points. OK"<<std::endl;
    else {
        std::cout<<"\t\t\tScanning axis contains only "<<NSc<<" point. Not enough for electric field profile."<<std::endl;
        return false;
    }

    std::cout<<GetTitle()<<" Data Test Passed. Processing..."<<std::endl;
    return true;
}

#ifdef USE_GUI
void ModuleEdgeField::PrintConfig(std::ofstream &conf_file) {
    conf_file<<"\n#Averaging the current for electric field profile from F_TLow to F_TLow+EV_Time";
    conf_file<<"\nEV_Time\t=\t"<<GetEV_Time();
}
void ModuleEdgeField::AddParameters(QVBoxLayout *layout) {

    QHBoxLayout* hlayout = new QHBoxLayout;
    hlayout->addStretch();
    hlayout->addWidget(new QLabel("Integration time"));
    ev_time = new QDoubleSpinBox();
    ev_time->setMinimum(0);
    ev_time->setMaximum(1000);
    ev_time->setSingleStep(0.01);
    ev_time->setDecimals(4);
    hlayout->addWidget(ev_time);
    hlayout->addWidget(new QLabel("ns"));
    layout->addLayout(hlayout);

}
void ModuleEdgeField::FillParameters() {
    ev_time->setValue(GetEV_Time());
}
void ModuleEdgeField::ToVariables() {
    SetEV_Time(ev_time->value());
}

#endif

/// Mobility with energy
Double_t ModuleEdgeField::Mu(Double_t E, Int_t Type) {

    /** \param[in] E Energy
         *  \param[in] type Charge carrier type. 1 - hole, 0 - electron.
          */

    if(Type==1) return config->mu0_holes()/(1.+config->mu0_holes()*E/config->v_sat());
    if(Type==0) return config->mu0_els()/sqrt(1.+config->mu0_els()*E/config->v_sat()*config->mu0_els()*E/config->v_sat());

}

/// Function used in the bisection method
Double_t ModuleEdgeField::ff(Double_t E, Double_t Uuu, Double_t a) {
    /**   Based on formula \f$U = A*(\mu_{els}+\mu_{holes})*E\f$.
         *    Where U - measured charge, E - electric field.
          */
    return Uuu-a*(Mu(E,1)+Mu(E,0))*E;
}

/// Implements bisection method for solving the equation
Double_t ModuleEdgeField::BiSectionMethod(Double_t eps, Double_t x1, Double_t x2, Double_t Uuu, Double_t a) {
    if(ff(x1,Uuu,a)==0) return x1;
    if(ff(x2,Uuu,a)==0) return x2;
    Double_t dx = x2-x1;
    Double_t xi;
    while(abs(ff(xi,Uuu,a))>eps) {
        dx = dx/2;
        xi = x1+dx;
        if(sgn(ff(x1,Uuu,a))!=sgn(ff(xi,Uuu,a))) continue;
        else x1=xi;
    }
    //std::cout<<"x = "<<xi<<" Fe = "<<ff(xi,Uuu,a)<<std::endl;
    return xi;
}

/// Sgn implementation
template <typename T> int ModuleEdgeField::sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

}


