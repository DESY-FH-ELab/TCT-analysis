#include "base.h"
#include <fstream>

//Qt includes
#include "ui_base.h"
#include "QFileDialog"
#include "ui_form_parameters.h"
#include "ui_form_sample.h"
#include "ui_form_folders.h"
#include "gui_folders.h"
#include "gui_sample.h"
#include "QMessageBox"
#include "QDateTime"

//TCT includes
#include <acquisition.h>
#include <scanning.h>
#include <sample.h>
#include <util.h>

//ROOT includes
#include <TSystem.h>

base::base(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::base)
{
    ui->setupUi(this);

#ifdef USE_GUI
    freopen ("execution.log","a",stdout);
    QDateTime datetime;
    std::cout<<"<-------------- NEW RUN -------------->"<<std::endl;
    std::cout<<"Run Started At: "<<datetime.currentDateTime().toString().toStdString()<<std::endl;
#endif

    //on_comboBox_activated(0);
    ui->buttonGroup_mode->setId(ui->mode_top,0);
    ui->buttonGroup_mode->setId(ui->mode_edge,1);
    ui->buttonGroup_mode->setId(ui->mode_bottom,2);
    on_buttonGroup_mode_buttonClicked(3);

    ui->buttonGroup_scanning->setId(ui->sc_x,1);
    ui->buttonGroup_scanning->setId(ui->sc_y,2);
    ui->buttonGroup_scanning->setId(ui->sc_z,3);

    ui->buttonGroup_optical->setId(ui->opt_x,1);
    ui->buttonGroup_optical->setId(ui->opt_y,2);
    ui->buttonGroup_optical->setId(ui->opt_z,3);

    ui->cc_ch_1->toggle();
    ui->cc_ch_2->toggle();
    ui->cc_ch_ph->toggle();
    ui->cc_ch_trig->toggle();

    TCT::util basic_config;
    std::ifstream basic_config_file("default.conf");
    if(!basic_config_file.is_open()) {
        error("default.conf file is absent. Please create file and put the DefaultFile variable in it.");
        exit(1);
    }
    basic_config.parse(basic_config_file);
    for( auto i : basic_config.ID_val()){
      if(i.first == "DefaultFile")		_DefConfigName = i.second;
    }
    if(_DefConfigName=="def") {
        error("default.conf doesn't contain the <b>DefaultFile</b> name (the default config file address).");
        exit(1);
    }
    read_config(_DefConfigName.c_str());
    fill_config();

}

base::~base()
{
#ifdef USE_GUI
    QDateTime datetime;
    std::cout<<"Run Finished At: "<<datetime.currentDateTime().toString().toStdString()<<std::endl;
    fclose(stdout);
#endif
    delete ui;
}

void base::on_buttonGroup_mode_buttonClicked(int index)
{
    ui->group_0->setVisible(!(bool)abs(index-0));
    ui->group_0->setGeometry(300,230,360,300);
    ui->group_1->setVisible(!(bool)abs(index-1));
    ui->group_1->setGeometry(300,230,360,300);
    ui->group_2->setVisible(!(bool)abs(index-2));
    ui->group_2->setGeometry(300,230,360,300);
}
void base::read_config(const char *config_file) {

    if(config_mode) delete config_mode;
    if(config_analysis) delete config_analysis;
    if(config_sample) delete config_sample;

    TCT::util analysis_card;
    std::ifstream ana_file(config_file);

    if(!ana_file.is_open()) {
        error("Analysis Config file open failed. Please check the filename in default.conf.");
        exit(1);
    }

    analysis_card.parse(ana_file);
    config_mode = new TCT::mode_selector(analysis_card.ID_val());
    config_analysis = new TCT::tct_config(analysis_card.ID_val());

    TCT::util sample_card;
    std::ifstream sample_file(config_mode->SampleCard().c_str());

    if(!sample_file.is_open()) {
        error("Sample Card file open failed. Please specify the existing filename.");
    }

    sample_card.parse(sample_file);
    config_sample = new TCT::sample(sample_card.ID_val());

    config_analysis->SetSampleThickness(config_sample->Thickness());
    config_analysis->SetOutSample_ID(config_sample->SampleID());

    ui->statusBar->showMessage(tr("Config File Opened"), 2000);

}
void base::fill_config() {

    // tct config

    ui->cc_ch_1->setChecked(config_analysis->CH_Det());
    if(config_analysis->CH_Det()) {
        ui->cc_ch_num_1->setCurrentIndex(config_analysis->CH_Det()-1);
        ui->tlow_1->setText(QString::number(config_analysis->FTlow()));
        ui->thigh_1->setText(QString::number(config_analysis->FThigh()));

    }
    ui->cc_ch_ph->setChecked(config_analysis->CH_PhDiode());
    if(config_analysis->CH_PhDiode()) {
        ui->cc_ch_num_ph->setCurrentIndex(config_analysis->CH_PhDiode()-1);
        ui->tlow_diode->setText(QString::number(config_analysis->FDLow()));
        ui->thigh_diode->setText(QString::number(config_analysis->FDHigh()));

    }
    ui->cc_ch_trig->setChecked(config_analysis->CH_Trig());
    if(config_analysis->CH_Trig()) {
        ui->cc_ch_num_trig->setCurrentIndex(config_analysis->CH_Trig()-1);
    }

    ui->buttonGroup_scanning->button(config_analysis->ScAxis())->setChecked(true);
    ui->buttonGroup_optical->button(config_analysis->OptAxis())->setChecked(true);

    ui->cc_volt_source->setCurrentIndex(config_analysis->VoltSource()-1);

    ui->save_sep_chrges->setChecked(config_analysis->FSeparateCharges());
    ui->save_sep_waveforms->setChecked(config_analysis->FSeparateWaveforms());
    ui->timebetween->setValue(config_analysis->Movements_dt());

    ui->buttonGroup_mode->button(config_analysis->TCT_Mode())->setChecked(true);
    on_buttonGroup_mode_buttonClicked(config_analysis->TCT_Mode());

    ui->top_focus->setChecked(config_analysis->DO_focus());
    ui->edge_focus->setChecked(config_analysis->DO_focus());
    ui->edge_depl->setChecked(config_analysis->DO_EdgeDepletion());
    ui->edge_profiles->setChecked(config_analysis->DO_EdgeVelocity());

    // oscilloscope config

}

void base::tovariables_config() {

    // tct part

    if(ui->cc_ch_1->isChecked()) {
        config_analysis->SetCH_Det(ui->cc_ch_num_1->currentText().toInt());
        config_analysis->SetFTlow(ui->tlow_1->text().toFloat());
        config_analysis->SetFThigh(ui->thigh_1->text().toFloat());
    }
    else config_analysis->SetCH_Det(0);

    if(ui->cc_ch_ph->isChecked()) {
        config_analysis->SetCH_PhDiode(ui->cc_ch_num_ph->currentText().toInt());
        config_analysis->SetFDLow(ui->tlow_diode->text().toFloat());
        config_analysis->SetFDHigh(ui->thigh_diode->text().toFloat());
    }
    else config_analysis->SetCH_PhDiode(0);

    if(ui->cc_ch_trig->isChecked()) {
        config_analysis->SetCH_Trig(ui->cc_ch_num_trig->currentText().toInt());
    }
    else config_analysis->SetCH_Trig(0);

    config_analysis->SetScAxis(ui->buttonGroup_scanning->checkedId());
    config_analysis->SetOptAxis(ui->buttonGroup_optical->checkedId());
    config_analysis->SetVoltSource(ui->cc_volt_source->currentText().toInt());

    config_analysis->SetFSeparateCharges(ui->save_sep_chrges->isChecked());
    config_analysis->SetFSeparateWaveforms(ui->save_sep_waveforms->isChecked());

    config_analysis->SetMovements_dt(ui->timebetween->value());

    config_analysis->SetTCT_Mode(ui->buttonGroup_mode->checkedId());

    config_analysis->SetDO_focus(ui->edge_focus->isChecked());
    config_analysis->SetDO_EdgeDepletion(ui->edge_depl->isChecked());
    config_analysis->SetDO_EdgeVelocity(ui->edge_profiles->isChecked());

    // oscilloscope part
}

void base::on_window_folders_clicked()
{
    Folders *folders = new Folders(this);

    folders->ui->cc_data->setText(config_analysis->DataFolder().c_str());
    folders->ui->cc_out->setText(config_analysis->OutFolder().c_str());

    if(folders->exec()) {
        config_analysis->SetDataFolder(folders->ui->cc_data->text().toStdString());
        config_analysis->SetOutFolder(folders->ui->cc_out->text().toStdString());
    }

    delete folders;
}

void base::on_window_sample_clicked()
{

    Sample *sampleUi = new Sample(this);

    sampleUi->ui->cc_sensor->setText(config_mode->SampleCard().c_str());
    sampleUi->ui->folder->setText(config_sample->Folder().c_str());
    sampleUi->ui->id->setText(config_sample->SampleID().c_str());
    sampleUi->ui->thickness->setText(QString::number(config_sample->Thickness()));

    if(sampleUi->exec()) {
        config_sample->SetFolder(sampleUi->ui->folder->text().toStdString());
        config_sample->SetSampleID(sampleUi->ui->id->text().toStdString());
        config_sample->SetThickness(sampleUi->ui->thickness->text().toFloat());
        config_analysis->SetSampleThickness(config_sample->Thickness());
        config_analysis->SetOutSample_ID(config_sample->SampleID());
        config_mode->SetSampleCard(sampleUi->ui->cc_sensor->text().toStdString());
    }

    delete sampleUi;

}

void base::on_window_parameters_clicked()
{
    QDialog *parameters = new QDialog(0,0);

    Ui_Parameters *parametersUi = new Ui_Parameters;
    parametersUi->setupUi(parameters);

    parametersUi->mobility_els->setText(QString::number(config_analysis->mu0_els()));
    parametersUi->mobility_holes->setText(QString::number(config_analysis->mu0_holes()));
    parametersUi->sat_vel->setText(QString::number(config_analysis->v_sat()));
    parametersUi->res_input->setText(QString::number(config_analysis->R_sensor()));
    parametersUi->res_input_diode->setText(QString::number(config_analysis->R_diode()));
    parametersUi->response_diode->setText(QString::number(config_analysis->RespPhoto()));
    parametersUi->splitter->setText(QString::number(config_analysis->light_split()));
    parametersUi->ampl->setText(QString::number(config_analysis->ampl()));
    parametersUi->epair->setText(QString::number(config_analysis->E_pair()));

    if(parameters->exec()) {
        config_analysis->Setmu0_els(parametersUi->mobility_els->text().toFloat());
        config_analysis->Setmu0_holes(parametersUi->mobility_holes->text().toFloat());
        config_analysis->Setv_sat(parametersUi->sat_vel->text().toFloat());
        config_analysis->SetR_sensor(parametersUi->res_input->text().toFloat());
        config_analysis->SetR_diode(parametersUi->res_input_diode->text().toFloat());
        config_analysis->SetRespPhoto(parametersUi->response_diode->text().toFloat());
        config_analysis->Setlight_split(parametersUi->splitter->text().toFloat());
        config_analysis->Setampl(parametersUi->ampl->text().toFloat());
        config_analysis->SetE_pair(parametersUi->epair->text().toFloat());
    }

    delete parametersUi;
    delete parameters;

}


void base::on_start_clicked()
{

    tovariables_config();

    if(ui->tct_single->isChecked()) {
        QStringList names;
        /*names = QFileDialog::getOpenFileNames(this,
                                                  tr("Select Files"), QString(config_analysis->DataFolder().c_str()), tr("TCT Files (*.tct)"));*/
        QFileDialog dialog;
        dialog.setNameFilter(tr("TCT Files (*.tct)"));
        dialog.setFileMode(QFileDialog::ExistingFiles);
        dialog.setDirectory(QString(config_analysis->DataFolder().c_str()));
        if(dialog.exec()) {
            names = dialog.selectedFiles();
            ui->statusBar->showMessage(QString("Number of files selected: %1").arg(names.length()),500);
        }
        else {
             ui->statusBar->showMessage(tr("No files selected"));
             return;
        }

        for(int i=0;i<names.length();i++) {
            ui->statusBar->showMessage(QString("Reading file %1").arg(names.at(i).split("/").last()),1000);
            char pathandfile[250];
            strcpy(pathandfile,names.at(i).toStdString().c_str());
            TCT::Scanning daq_data;
            bool read = daq_data.ReadTCT(pathandfile,config_analysis);
            if(!read) {
                ui->statusBar->showMessage(QString("Processing of file %1 failed. Skipping").arg(names.at(i).split("/").last()));
                continue;
            }
        }

    }
    else {
        if(config_analysis->DataFolder() == "def") {
            error("No data folder was specified in analysis card. Check your analysis card, that \"DataFolder = ...\" is specified correctly");
            //exit(1);
        }
        ui->statusBar->showMessage(tr("Searching for Data in Folder"));


        void *dirp = gSystem->OpenDirectory(config_analysis->DataFolder().c_str());
        if (!dirp) {
            error("Data Folder not found. Check \"DataFolder\" in analysis card.");
            //exit(1);
        }
        std::string path = config_analysis->DataFolder();
        if (path.length() > 0) {
            std::string::iterator it = path.end() - 1;
            if (*it != '/') {
                path.append("/");
            }
        }
        config_analysis->SetDataFolder(path);

        const char *infile;
        while((infile = gSystem->GetDirEntry(dirp))) {

            if (strstr(infile,".tct")) {
                char pathandfile[250];
                strcpy(pathandfile,config_analysis->DataFolder().c_str());
                strcat(pathandfile,infile);
                ui->statusBar->showMessage(QString("Reading file %1").arg(infile));

                TCT::Scanning daq_data;
                bool read = daq_data.ReadTCT(pathandfile,config_analysis);
                if(!read) {ui->statusBar->showMessage(QString("Processing of file %1 failed. Skipping").arg(infile));  continue;}
            }

        }
    }
    ui->statusBar->showMessage(tr("Finished"));

}

void base::on_actionChange_config_triggered()
{

    QString cpart = "";
    cpart = QFileDialog::getOpenFileName(this,
                                              tr("Open Config File"), ".", tr("Text Files (*.txt)"));
    if(cpart=="") return;
    else {
        read_config(cpart.toStdString().c_str());
        fill_config();
    }

}

void base::on_actionSave_config_triggered()
{
    QString cpart = "";
    cpart = QFileDialog::getSaveFileName(this,
                                              tr("Save Config File"), tr("."));
    if(cpart=="") return;

    std::ofstream conf_file(cpart.toStdString().c_str(),std::fstream::trunc);

    if(!conf_file.is_open()) {
        error("Can not open the text file to write config. Aborting operation.");
        return;
    }

    conf_file<<"#comments have to start with a\n";
    conf_file<<"#put group key words in []\n";
    conf_file<<"#always specify ID, TAB, \"=\", TAB, value";

    conf_file<<"\n\n[General]";
    conf_file<<"\nProjectFolder\t=\t./"; //FIXME
    conf_file<<"\nDataFolder\t=\t"<<config_analysis->DataFolder();
    conf_file<<"\nOutfolder\t=\t"<<config_analysis->OutFolder();
    conf_file<<"\n\n\
#Set acq mode.\n\
# 0 - taking the sets of single measurements (*.txt or *.raw files by oscilloscope). Settings are in [Analysis]\n\
# 1 - taking the data from *.tct file produced by DAQ software. Settings are in [Scanning]";
    conf_file<<"\nMode\t=\t"<<config_mode->Mode();

    conf_file<<"\n\n[Scanning]";
    conf_file<<"\n#Channels of oscilloscope connected to detector, photodiode, trigger. Put numbers 1,2,3,4 - corresponding to channels, no such device connected put 0.";
    conf_file<<"\nCH_Detector\t=\t"<<config_analysis->CH_Det();
    conf_file<<"\n#Turning on of the Photodiode channel also adds normalisation to all scans";
    conf_file<<"\nCH_Photodiode\t=\t"<<config_analysis->CH_PhDiode();
    conf_file<<"\nCH_Trigger\t=\t"<<config_analysis->CH_Trig();
    conf_file<<"\n#Set optical Axis. 1-x,2-y,3-z";
    conf_file<<"\nOptical_Axis\t=\t"<<config_analysis->OptAxis();
    conf_file<<"\n#Set scanning Axis. 1-x,2-y,3-z";
    conf_file<<"\nScanning_Axis\t=\t"<<config_analysis->ScAxis();
    conf_file<<"\n#Set voltage source number (1 or 2)";
    conf_file<<"\nVoltage_Source\t=\t"<<config_analysis->VoltSource();
    conf_file<<"\n#Time between stage movements in seconds.";
    conf_file<<"\nMovements_dt\t=\t"<<config_analysis->Movements_dt();

    conf_file<<"\n#Perform next operations. Analysis will start only if all needed data is present:";
    conf_file<<"\n# 0-top,1-edge,2-bottom";
    conf_file<<"\nTCT_Mode\t=\t"<<config_analysis->TCT_Mode();
    conf_file<<"\n\n#Scanning over optical and perpendiculr to strip axes (or along the detector depth in case of edge-tct), fitting the best position.";

    conf_file<<"\nFocus_Search\t=\t"<<config_analysis->DO_focus();
    conf_file<<"\n#search for depletion voltage";
    conf_file<<"\nEdgeDepletionVoltage\t=\t"<<config_analysis->DO_EdgeDepletion();
    conf_file<<"\n#extracting the velocity and electric field profiles";
    conf_file<<"\nEdgeVelocityProfile\t=\t"<<config_analysis->DO_EdgeVelocity();

    conf_file<<"\n\n#Integrate sensor signal from TimeSensorLow to TimeSensorHigh - ns";
    conf_file<<"\nTimeSensorLow\t=\t"<<config_analysis->FTlow();
    conf_file<<"\nTimeSensorHigh\t=\t"<<config_analysis->FThigh();
    conf_file<<"\n#Integrate photodiode signal from TimeDiodeLow to TimeDiodeHigh - ns";
    conf_file<<"\nTimeDiodeLow\t=\t"<<config_analysis->FDLow();
    conf_file<<"\nTimeDiodeHigh\t=\t"<<config_analysis->FDHigh();
    conf_file<<"\n\n#Save charge, normed charge and photodiode charge for each Z, voltage";
    conf_file<<"\nSaveSeparateCharges\t=\t"<<config_analysis->FSeparateCharges();
    conf_file<<"\n#Save waveforms for each position and voltage";
    conf_file<<"\nSaveSeparateWaveforms\t=\t"<<config_analysis->FSeparateWaveforms();
    conf_file<<"\n#Averaging the current for electric field profile from F_TLow to F_TLow+EV_Time";
    conf_file<<"\nEV_Time\t=\t"<<config_analysis->EV_Time();


    conf_file<<"\n\n[Parameters]";
    conf_file<<"\n#low-field mobility for electrons, cm2*V^-1*s^-1";
    conf_file<<"\nMu0_Electrons\t=\t"<<config_analysis->mu0_els();
    conf_file<<"\n#low-field mobility for holes, cm2*V^-1*s^-1";
    conf_file<<"\nMu0_Holes\t=\t"<<config_analysis->mu0_holes();
    conf_file<<"\n#saturation velocity cm/s";
    conf_file<<"\nSaturationVelocity\t=\t"<<config_analysis->v_sat();
    conf_file<<"\n# amplifier amplification";
    conf_file<<"\nAmplification\t=\t"<<config_analysis->ampl();
    conf_file<<"\n# factor between charge in sensor and photodiode due to light splitting: Nsensor/Ndiode";
    conf_file<<"\nLightSplitter\t=\t"<<config_analysis->light_split();
    conf_file<<"\n# resistance of the sensor and diode output, Ohm";
    conf_file<<"\nResistanceSensor\t=\t"<<config_analysis->R_sensor();
    conf_file<<"\nResistancePhotoDetector\t=\t"<<config_analysis->R_diode();
    conf_file<<"\n# pohotodetector responce for certain wavelength, A/W";
    conf_file<<"\nResponcePhotoDetector\t=\t"<<config_analysis->RespPhoto();
    conf_file<<"\n# electron-hole pair creation energy, eV";
    conf_file<<"\nEnergyPair\t=\t"<<config_analysis->E_pair();

    conf_file<<"\n\n[Sensor]";
    conf_file<<"\nSampleCard\t=\t"<<config_mode->SampleCard();

    conf_file.close();

    ui->statusBar->showMessage(tr("Config File Saved"), 2000);
}

void base::error(const char *text) {
    QMessageBox *messageBox = new QMessageBox(this);
    messageBox->critical(0,"Error",text);
    messageBox->setFixedSize(500,200);
    delete messageBox;
}
void base::not_run() {
    std::cout<<"before"<<std::endl;
    std::vector<TCT::acquisition_single> AllAcqs;
    TCT::acquisition_avg AcqAvg(AllAcqs[0].Nsamples());
    std::cout<<"before"<<std::endl;
}
