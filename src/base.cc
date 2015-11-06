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
#include "qdebugstream.h"
#include "QProgressDialog"
#include "qdebug.h"
#include "QDirIterator"
#include "QProcess"

//TCT includes
#include "acquisition.h"
#include "measurement.h"
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

    ui->buttonGroup_input_type->setId(ui->mode_txt,0);
    ui->buttonGroup_input_type->setId(ui->mode_raw,1);

#ifndef USE_LECROY_RAW
    ui->mode_raw->setEnabled(false);
    ui->mode_raw->setToolTip("The code was built without LeCroy RAW data converter library.");
#endif

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
    if(config_mode) delete config_mode;
    if(config_tct) delete config_tct;
    if(config_analysis) delete config_analysis;
    if(config_sample) delete config_sample;
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
    if(config_tct) delete config_tct;
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
    config_tct = new TCT::tct_config(analysis_card.ID_val());
    config_analysis = new TCT::analysis(analysis_card.ID_val());

    TCT::util sample_card;
    std::ifstream sample_file(config_mode->SampleCard().c_str());

    if(!sample_file.is_open()) {
        error("Sample Card file open failed. Please specify the existing filename.");
    }

    sample_card.parse(sample_file);
    config_sample = new TCT::sample(sample_card.ID_val());

    config_tct->SetSampleThickness(config_sample->Thickness());
    config_tct->SetOutSample_ID(config_sample->SampleID());

    ui->statusBar->showMessage(tr("Config File Opened"), 2000);

}
void base::fill_config() {

    ui->modes->setCurrentIndex(config_mode->Mode());

    // tct config

    ui->cc_ch_1->setChecked(config_tct->CH_Det());
    if(config_tct->CH_Det()) {
        ui->cc_ch_num_1->setCurrentIndex(config_tct->CH_Det()-1);
        ui->tlow_1->setText(QString::number(config_tct->FTlow()));
        ui->thigh_1->setText(QString::number(config_tct->FThigh()));

    }
    ui->cc_ch_ph->setChecked(config_tct->CH_PhDiode());
    if(config_tct->CH_PhDiode()) {
        ui->cc_ch_num_ph->setCurrentIndex(config_tct->CH_PhDiode()-1);
        ui->tlow_diode->setText(QString::number(config_tct->FDLow()));
        ui->thigh_diode->setText(QString::number(config_tct->FDHigh()));

    }
    ui->cc_ch_trig->setChecked(config_tct->CH_Trig());
    if(config_tct->CH_Trig()) {
        ui->cc_ch_num_trig->setCurrentIndex(config_tct->CH_Trig()-1);
    }

    ui->buttonGroup_scanning->button(config_tct->ScAxis())->setChecked(true);
    ui->buttonGroup_optical->button(config_tct->OptAxis())->setChecked(true);

    ui->cc_volt_source->setCurrentIndex(config_tct->VoltSource()-1);

    ui->save_sep_chrges->setChecked(config_tct->FSeparateCharges());
    ui->save_sep_waveforms->setChecked(config_tct->FSeparateWaveforms());
    ui->timebetween->setValue(config_tct->Movements_dt());

    ui->buttonGroup_mode->button(config_tct->TCT_Mode())->setChecked(true);
    on_buttonGroup_mode_buttonClicked(config_tct->TCT_Mode());

    ui->top_focus->setChecked(config_tct->DO_focus());
    ui->edge_focus->setChecked(config_tct->DO_focus());
    ui->edge_depl->setChecked(config_tct->DO_EdgeDepletion());
    ui->edge_profiles->setChecked(config_tct->DO_EdgeVelocity());

    // oscilloscope config

    ui->MaxAcqs->setValue(config_analysis->MaxAcqs());
    ui->Noise_Cut->setValue(config_analysis->Noise_Cut());
    ui->NoiseEnd_Cut->setValue(config_analysis->NoiseEnd_Cut());

    ui->S2n_Cut->setValue(config_analysis->S2n_Cut());
    ui->S2n_Ref->setValue(config_analysis->S2n_Ref());
    ui->AmplNegLate_Cut->setValue(config_analysis->AmplNegLate_Cut());
    ui->AmplPosLate_Cut->setValue(config_analysis->AmplPosLate_Cut());
    ui->AmplNegEarly_Cut->setValue(config_analysis->AmplNegEarly_Cut());
    ui->AmplPosEarly_Cut->setValue(config_analysis->AmplPosEarly_Cut());
    ui->dosmearing->setChecked(config_analysis->DoSmearing());
    ui->AddNoise->setValue(config_analysis->AddNoise());
    ui->AddJitter->setValue(config_analysis->AddJitter());
    ui->savetofile->setChecked(config_analysis->SaveToFile());
    ui->savesingles->setChecked(config_analysis->SaveSingles());
    ui->PrintEvent->setValue(config_analysis->PrintEvent());

    ui->buttonGroup_input_type->button(config_analysis->LeCroyRAW())->click();
#ifndef USE_LECROY_RAW
    ui->mode_txt->click();
#endif

    ui->Nsamples_start->setEnabled(false);
    ui->Nsamples_end->setEnabled(false);
    ui->PrePulseInterval->setEnabled(false);

         /*
            Nsamples_start	=	50	# in Samples
            Nsamples_end	=	50	# in Samples
            PrePulseInterval	= 	50	# in ns.
            */

}

void base::tovariables_config() {

    config_mode->SetMode(ui->modes->currentIndex());

    // tct part

    if(ui->cc_ch_1->isChecked()) {
        config_tct->SetCH_Det(ui->cc_ch_num_1->currentText().toInt());
        config_tct->SetFTlow(ui->tlow_1->text().toFloat());
        config_tct->SetFThigh(ui->thigh_1->text().toFloat());
    }
    else config_tct->SetCH_Det(0);

    if(ui->cc_ch_ph->isChecked()) {
        config_tct->SetCH_PhDiode(ui->cc_ch_num_ph->currentText().toInt());
        config_tct->SetFDLow(ui->tlow_diode->text().toFloat());
        config_tct->SetFDHigh(ui->thigh_diode->text().toFloat());
    }
    else config_tct->SetCH_PhDiode(0);

    if(ui->cc_ch_trig->isChecked()) {
        config_tct->SetCH_Trig(ui->cc_ch_num_trig->currentText().toInt());
    }
    else config_tct->SetCH_Trig(0);

    config_tct->SetScAxis(ui->buttonGroup_scanning->checkedId());
    config_tct->SetOptAxis(ui->buttonGroup_optical->checkedId());
    config_tct->SetVoltSource(ui->cc_volt_source->currentText().toInt());

    config_tct->SetFSeparateCharges(ui->save_sep_chrges->isChecked());
    config_tct->SetFSeparateWaveforms(ui->save_sep_waveforms->isChecked());

    config_tct->SetMovements_dt(ui->timebetween->value());

    config_tct->SetTCT_Mode(ui->buttonGroup_mode->checkedId());

    config_tct->SetDO_focus(ui->edge_focus->isChecked());
    config_tct->SetDO_EdgeDepletion(ui->edge_depl->isChecked());
    config_tct->SetDO_EdgeVelocity(ui->edge_profiles->isChecked());

    // oscilloscope part

    config_analysis->SetMaxAcqs(ui->MaxAcqs->value());
    config_analysis->SetNoise_Cut(ui->Noise_Cut->value());
    config_analysis->SetNoiseEnd_Cut(ui->NoiseEnd_Cut->value());
    config_analysis->SetS2n_Cut(ui->S2n_Cut->value());
    config_analysis->SetS2n_Ref(ui->S2n_Ref->value());
    config_analysis->SetAmplNegLate_Cut(ui->AmplNegLate_Cut->value());
    config_analysis->SetAmplPosLate_Cut(ui->AmplPosLate_Cut->value());
    config_analysis->SetAmplNegEarly_Cut(ui->AmplNegEarly_Cut->value());
    config_analysis->SetAmplPosEarly_Cut(ui->AmplPosEarly_Cut->value());
    config_analysis->SetDoSmearing(ui->dosmearing->isChecked());
    config_analysis->SetAddNoise(ui->AddNoise->value());
    config_analysis->SetAddJitter(ui->AddJitter->value());
    config_analysis->SetSaveToFile(ui->savetofile->isChecked());
    config_analysis->SetSaveSingles(ui->savesingles->isChecked());
    config_analysis->SetPrintEvent(ui->PrintEvent->value());
    config_analysis->SetLeCroyRAW((bool)(ui->buttonGroup_input_type->checkedId()));
}

void base::on_window_folders_clicked()
{
    Folders *folders = new Folders(this);

    folders->ui->cc_data->setText(config_tct->DataFolder().c_str());
    folders->ui->cc_out->setText(config_tct->OutFolder().c_str());

    if(folders->exec()) {
        config_tct->SetDataFolder(folders->ui->cc_data->text().toStdString());
        config_tct->SetOutFolder(folders->ui->cc_out->text().toStdString());
        if(config_mode->Mode()==0) {
            config_analysis->SetDataFolder(config_tct->DataFolder());
            config_analysis->SetOutFolder(config_tct->OutFolder());
        }
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
        config_tct->SetSampleThickness(config_sample->Thickness());
        config_tct->SetOutSample_ID(config_sample->SampleID());
        config_mode->SetSampleCard(sampleUi->ui->cc_sensor->text().toStdString());
    }

    delete sampleUi;

}

void base::on_window_parameters_clicked()
{
    QDialog *parameters = new QDialog(0,0);

    Ui_Parameters *parametersUi = new Ui_Parameters;
    parametersUi->setupUi(parameters);

    parametersUi->mobility_els->setText(QString::number(config_tct->mu0_els()));
    parametersUi->mobility_holes->setText(QString::number(config_tct->mu0_holes()));
    parametersUi->sat_vel->setText(QString::number(config_tct->v_sat()));
    parametersUi->res_input->setText(QString::number(config_tct->R_sensor()));
    parametersUi->res_input_diode->setText(QString::number(config_tct->R_diode()));
    parametersUi->response_diode->setText(QString::number(config_tct->RespPhoto()));
    parametersUi->splitter->setText(QString::number(config_tct->light_split()));
    parametersUi->ampl->setText(QString::number(config_tct->ampl()));
    parametersUi->epair->setText(QString::number(config_tct->E_pair()));

    if(parameters->exec()) {
        config_tct->Setmu0_els(parametersUi->mobility_els->text().toFloat());
        config_tct->Setmu0_holes(parametersUi->mobility_holes->text().toFloat());
        config_tct->Setv_sat(parametersUi->sat_vel->text().toFloat());
        config_tct->SetR_sensor(parametersUi->res_input->text().toFloat());
        config_tct->SetR_diode(parametersUi->res_input_diode->text().toFloat());
        config_tct->SetRespPhoto(parametersUi->response_diode->text().toFloat());
        config_tct->Setlight_split(parametersUi->splitter->text().toFloat());
        config_tct->Setampl(parametersUi->ampl->text().toFloat());
        config_tct->SetE_pair(parametersUi->epair->text().toFloat());
    }

    delete parametersUi;
    delete parameters;

}


void base::on_start_clicked()
{
    if(ui->modes->currentIndex() == 0) start_osc();
    if(ui->modes->currentIndex() == 1) start_tct();
}

void base::start_tct() {
    tovariables_config();
    QStringList names;

    if(ui->tct_single->isChecked()) {

        QFileDialog dialog(this);
        dialog.setNameFilter(tr("TCT Files (*.tct)"));
        dialog.setFileMode(QFileDialog::ExistingFiles);
        dialog.setDirectory(QString(config_tct->DataFolder().c_str()));
        dialog.setModal(true);
        if(dialog.exec()) {
            names = dialog.selectedFiles();
        }
        else {
             ui->statusBar->showMessage(tr("No files selected"));
             return;
        }
    }
    else {
        if(config_tct->DataFolder() == "def") {
            error("No data folder was specified in analysis card. Check your analysis card, that \"DataFolder = ...\" is specified correctly");
            return;
        }

        QDir dirp(config_tct->DataFolder().c_str());
        if(!dirp.exists()) {
            error("Data Folder not found. Check \"DataFolder\" in analysis card.");
            return;
        }
        ui->statusBar->showMessage(tr("Searching for Data in Folder"));
        QStringList filters;
        filters << "*.tct";
        QStringList filenames = dirp.entryList(filters);
        for(int i=0;i<filenames.length();i++) names<<dirp.absoluteFilePath(filenames.at(i));
    }
    if(names.length()==0) {
        error("No files in specified folder!");
        return;
    }

    ui->statusBar->showMessage(QString("Number of files selected: %1").arg(names.length()),500);
    progress = new Ui::ConsoleOutput(names.length(),this);
    connect(progress,SIGNAL(OpenTBrowser()),this,SLOT(on_tbrowser_clicked()));

    // begin redirecting output
    std::ofstream log_file("execution.log",std::fstream::app);
    QDebugStream *debug = new QDebugStream(std::cout, log_file, progress->Console());
    print_run(true);

    for(int i=0;i<names.length();i++) {
        progress->setValue(i);
        if(progress->wasCanceled()) break;
        ui->statusBar->showMessage(QString("Reading file %1").arg(names.at(i).split("/").last()),1000);
        char pathandfile[250];
        strcpy(pathandfile,names.at(i).toStdString().c_str());
        TCT::Scanning daq_data;
        bool read = daq_data.ReadTCT(pathandfile,config_tct);
        if(!read) {
            ui->statusBar->showMessage(QString("Processing of file %1 failed. Skipping").arg(names.at(i).split("/").last()));
            continue;
        }
    }
    progress->setValue(names.length());
    //progress->finished(names.length());
    print_run(false);
    delete debug;
    connect(progress,SIGNAL(canceled()),this,SLOT(deleteprogress()));
    log_file.close();
    // end redirecting output

    ui->statusBar->showMessage(tr("Finished"));
}

void base::start_osc()
{
    tovariables_config();

    if(config_analysis->DataFolder() == "def") {
        error("No data folder was specified in analysis card. Check your analysis card, that \"DataFolder = ...\" is specified correctly");
        return;
    }

    progress_osc = new Ui::ConsoleOsc(this);
    connect(progress_osc,SIGNAL(OpenTBrowser()),this,SLOT(on_tbrowser_clicked()));
    progress_osc->show();
    // begin redirecting output
    std::ofstream log_file("execution.log",std::fstream::app);
    QDebugStream *debug = new QDebugStream(std::cout, log_file, progress_osc->Console());
    print_run(true);

    ui->statusBar->showMessage(tr("Searching data in folders"));

    std::cout << "The sample's data folder is = " << config_analysis->DataFolder() << ", searching data in subfolder(s) " <<  std::endl;
    std::map<std::string, std::vector<std::string>> folder_struc;
    std::vector<std::string> dirs;
    std::vector<std::string> dirs2;	// push folder for every subfolder (linearisation of folder matrix)
    std::vector<std::string> subdirs2;	// push folder for every subfolder (linearisation of folder matrix)
    std::vector<std::string> pathndirs;
    std::vector<int> Nsubdirs;

    void *dirp = gSystem->OpenDirectory(config_analysis->DataFolder().c_str());
    if (!dirp) {
        error("Data Folder not found. Check \"DataFolder\" in analysis card.");
        return;
    }
    char *direntry;
    uint32_t counterdir   = 0;
    uint32_t countersubdir= 0;
    while ((direntry = (char*)gSystem->GetDirEntry(dirp))) {
        if ( strstr(direntry,"..") || strstr(direntry,".") ) continue;
        counterdir++;
        //std::cout << "- " << direntry << std::endl;
        dirs.push_back((std::string)direntry);


        std::string subDataFolder	= config_analysis->DataFolder() + "/" + direntry;
        //std::cout << "Subdir = " << subDataFolder << std::endl;
        std::vector<std::string> subdirs;
        void *subdirp = gSystem->OpenDirectory(subDataFolder.c_str());
        char *subdirentry;
        while ((subdirentry = (char*)gSystem->GetDirEntry(subdirp))) {
            if ( strstr(subdirentry,"..") || strstr(subdirentry,"." ) ) continue;
            countersubdir++;
            //std::cout << "-- " << subdirentry << std::endl;
            subdirs.push_back(subdirentry);
            dirs2.push_back((std::string)direntry);
            subdirs2.push_back((std::string)subdirentry);
            std::string fullpath = config_analysis->DataFolder() + "/" + direntry + "/" + subdirentry + "/";
            //std::cout << fullpath << std::endl;
            pathndirs.push_back(fullpath);

        }
        Nsubdirs.push_back(countersubdir);
        countersubdir = 0;

        folder_struc[dirs[counterdir-1]] = subdirs;


    }

    for (auto i : Nsubdirs) {
        countersubdir +=i;
        //std::cout << " i = " << i << std::endl;
    }

    std::cout << " Found the following folder structure: " << std::endl;
    for(auto i : folder_struc) {
        for(auto j : i.second)
            std::cout << i.first << " " <<  j << " " << "\n";
    }

    ui->statusBar->showMessage(QString("In total, found %1 subfolder(s)").arg(countersubdir));
    std::cout << " In total, found " << countersubdir << " subfolder(s) " << std::endl;

    uint32_t counter;
    //if(countersubdir > 0) {
    counter=0;
    while(1){


#ifdef DEBUG
        std::cout << " Start with subfolder # " << counter << std::endl;
#endif
        //int i = 0; i < countersubdir; i++) {
        // create vec with acq_singles in it
        std::vector<TCT::acquisition_single> AllAcqs;

        if(countersubdir == 0) {
            // check if DataFolder() ends on "/", if not, add it
            std::string path = config_analysis->DataFolder();
            if (path.length() > 0) {
                std::string::iterator it = path.end() - 1;
                if (*it != '/') {
                    path.append("/");
                }
            }
            config_analysis->SetDataFolder(path);
            // push DataFolder to pathndirs (this is a hack to cope with Nsubdir == 0...)
            pathndirs.push_back(config_analysis->DataFolder());
        }
        // create measurement object from one subdir for each cycle
        TCT::measurement meas(pathndirs[counter]);

        if(!meas.AcqsLoader(&AllAcqs, config_analysis->MaxAcqs(),config_analysis->LeCroyRAW())) {
            std::cout << " Folder empty! Skipping folder" << std::endl;
            counter++;
            if( countersubdir == 0) break;
            if(counter == countersubdir) break;
            continue;
        };

        // now create instance of avg acquisition using Nsamples from loaded files
        TCT::acquisition_avg AcqAvg(AllAcqs[0].Nsamples());
        AcqAvg.SetPolarity(AllAcqs[0].Polarity());

        //now analyse all acquisitions
        int Nselected = 0;

#ifdef DEBUG
        std::cout << "Size of AllAcqs = " << AllAcqs.size() << std::endl;
#endif

        AcqAvg.SetNanalysed(AllAcqs.size());
        for(uint32_t i_acq = 0; i_acq < AllAcqs.size(); i_acq++){

#ifdef DEBUG
            std::cout << " - Start with Acq #" << i_acq << std::endl;
#endif

            TCT::acquisition_single* acq = &AllAcqs[i_acq];
            if(config_analysis->DoSmearing()) config_analysis->AcqsSmearer(acq, config_analysis->AddNoise(), false);
            config_analysis->AcqsAnalyser(acq, i_acq, &AcqAvg);
            if(config_analysis->DoSmearing()) config_analysis->AcqsSmearer(acq, false, config_analysis->AddJitter()); // AcqsAnalyser removes jitter by determining each acqs delay. Hence, to add jitter, delay has to be manipulated after AcqsAnalyser (and before filling of profile

#ifdef DEBUG
            std::cout << *acq << std::endl;
#endif

            if( config_analysis->AcqsSelecter(acq) ) {
                Nselected++;
                acq->SetSelect(true);
            }
            AcqAvg.SetNselected(Nselected);
            config_analysis->AcqsProfileFiller(acq, &AcqAvg);

        } // end fot AllAcqs.size()

        //std::cout << "Mean s2nval = " << AcqAvg.M_V_S2nval() << std::endl;

        config_analysis->SetOutSample_ID(config_sample->SampleID());
        if(countersubdir > 0) {
            config_analysis->SetOutSubFolder(dirs2[counter]);
            config_analysis->SetOutSubsubFolder(subdirs2[counter]);
        } else {
            config_analysis->SetOutSubFolder("def");
            config_analysis->SetOutSubsubFolder("def");
        }


        if(config_analysis->SaveToFile())
            if(countersubdir > 0) config_analysis->AcqsWriter(&AllAcqs, &AcqAvg, true);
            else config_analysis->AcqsWriter(&AllAcqs, &AcqAvg, false);

        std::cout << "   Nselected = " << Nselected << std::endl;
        std::cout << "   ratio of selected acqs = " << Nselected << " / " << AllAcqs.size() << " = " << (float)Nselected/AllAcqs.size()*100. << "%\n\n" << std::endl;

        // now take care of memory management
        // delete remaning TH1Fs in acquisition_single and then clear AllAcqs
        for(int j = 0; j < AllAcqs.size(); j++) {
            AllAcqs[j].Clear();
        }
        AllAcqs.clear();

        counter++;
        if( countersubdir == 0) break;
        if(counter == countersubdir) break;
#ifdef DEBUG
        std::cout << "   AllAcqs has " << AllAcqs.size() << " objects left" << std::endl;
#endif
    }
    progress_osc->SetButtonEnabled();
    print_run(false);
    delete debug;
    connect(progress_osc,SIGNAL(rejected()),this,SLOT(deleteprogress_osc()));
    log_file.close();
    // end redirecting output

    ui->statusBar->showMessage(tr("Finished"));
}

void base::deleteprogress() {
    delete progress;
    progress = NULL;
}
void base::deleteprogress_osc() {
    delete progress_osc;
    progress_osc = NULL;
}

void base::on_actionChange_config_triggered()
{

    QString cpart = "";
    cpart = QFileDialog::getOpenFileName(this,
                                              tr("Open Config File"), "../testanalysis", tr("Text Files (*.txt)"));
    if(cpart=="") return;
    else {
        read_config(cpart.toStdString().c_str());
        fill_config();
    }

}

void base::on_actionSave_config_triggered()
{
    tovariables_config();

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
    conf_file<<"\nDataFolder\t=\t"<<config_tct->DataFolder();
    conf_file<<"\nOutfolder\t=\t"<<config_tct->OutFolder();
    conf_file<<"\n\n\
#Set acq mode.\n\
# 0 - taking the sets of single measurements (*.txt or *.raw files by oscilloscope). Settings are in [Analysis]\n\
# 1 - taking the data from *.tct file produced by DAQ software. Settings are in [Scanning]";
    conf_file<<"\nMode\t=\t"<<config_mode->Mode();

    conf_file<<"\n\n[Analysis]";
    conf_file<<"\nMaxAcqs\t=\t"<<config_analysis->MaxAcqs();
    conf_file<<"\nNoise_Cut\t=\t"<<config_analysis->Noise_Cut();
    conf_file<<"\nNoiseEnd_Cut\t=\t"<<config_analysis->NoiseEnd_Cut();
    conf_file<<"\nS2n_Cut\t=\t"<<config_analysis->S2n_Cut();
    conf_file<<"\nS2n_Ref\t=\t"<<config_analysis->S2n_Ref();
    conf_file<<"\nAmplNegLate_Cut\t=\t"<<config_analysis->AmplNegLate_Cut();
    conf_file<<"\nAmplPosLate_Cut\t=\t"<<config_analysis->AmplPosLate_Cut();
    conf_file<<"\nAmplNegEarly_Cut\t=\t"<<config_analysis->AmplNegEarly_Cut();
    conf_file<<"\nAmplPosEarly_Cut\t=\t"<<config_analysis->AmplPosEarly_Cut();
    conf_file<<"\nDoSmearing\t=\t"<<config_analysis->DoSmearing();
    conf_file<<"\nAddNoise\t=\t"<<config_analysis->AddNoise();
    conf_file<<"\nAddJitter\t=\t"<<config_analysis->AddJitter();
    conf_file<<"\nSaveToFile\t=\t"<<config_analysis->SaveToFile();
    conf_file<<"\nSaveSingles\t=\t"<<config_analysis->SaveSingles();
    conf_file<<"\nPrintEvent\t=\t"<<config_analysis->PrintEvent();
    conf_file<<"\nLeCroyRAW\t=\t"<<config_analysis->LeCroyRAW();

    conf_file<<"\n\n[Scanning]";
    conf_file<<"\n#Channels of oscilloscope connected to detector, photodiode, trigger. Put numbers 1,2,3,4 - corresponding to channels, no such device connected put 0.";
    conf_file<<"\nCH_Detector\t=\t"<<config_tct->CH_Det();
    conf_file<<"\n#Turning on of the Photodiode channel also adds normalisation to all scans";
    conf_file<<"\nCH_Photodiode\t=\t"<<config_tct->CH_PhDiode();
    conf_file<<"\nCH_Trigger\t=\t"<<config_tct->CH_Trig();
    conf_file<<"\n#Set optical Axis. 1-x,2-y,3-z";
    conf_file<<"\nOptical_Axis\t=\t"<<config_tct->OptAxis();
    conf_file<<"\n#Set scanning Axis. 1-x,2-y,3-z";
    conf_file<<"\nScanning_Axis\t=\t"<<config_tct->ScAxis();
    conf_file<<"\n#Set voltage source number (1 or 2)";
    conf_file<<"\nVoltage_Source\t=\t"<<config_tct->VoltSource();
    conf_file<<"\n#Time between stage movements in seconds.";
    conf_file<<"\nMovements_dt\t=\t"<<config_tct->Movements_dt();

    conf_file<<"\n#Perform next operations. Analysis will start only if all needed data is present:";
    conf_file<<"\n# 0-top,1-edge,2-bottom";
    conf_file<<"\nTCT_Mode\t=\t"<<config_tct->TCT_Mode();
    conf_file<<"\n\n#Scanning over optical and perpendiculr to strip axes (or along the detector depth in case of edge-tct), fitting the best position.";

    conf_file<<"\nFocus_Search\t=\t"<<config_tct->DO_focus();
    conf_file<<"\n#search for depletion voltage";
    conf_file<<"\nEdgeDepletionVoltage\t=\t"<<config_tct->DO_EdgeDepletion();
    conf_file<<"\n#extracting the velocity and electric field profiles";
    conf_file<<"\nEdgeVelocityProfile\t=\t"<<config_tct->DO_EdgeVelocity();

    conf_file<<"\n\n#Integrate sensor signal from TimeSensorLow to TimeSensorHigh - ns";
    conf_file<<"\nTimeSensorLow\t=\t"<<config_tct->FTlow();
    conf_file<<"\nTimeSensorHigh\t=\t"<<config_tct->FThigh();
    conf_file<<"\n#Integrate photodiode signal from TimeDiodeLow to TimeDiodeHigh - ns";
    conf_file<<"\nTimeDiodeLow\t=\t"<<config_tct->FDLow();
    conf_file<<"\nTimeDiodeHigh\t=\t"<<config_tct->FDHigh();
    conf_file<<"\n\n#Save charge, normed charge and photodiode charge for each Z, voltage";
    conf_file<<"\nSaveSeparateCharges\t=\t"<<config_tct->FSeparateCharges();
    conf_file<<"\n#Save waveforms for each position and voltage";
    conf_file<<"\nSaveSeparateWaveforms\t=\t"<<config_tct->FSeparateWaveforms();
    conf_file<<"\n#Averaging the current for electric field profile from F_TLow to F_TLow+EV_Time";
    conf_file<<"\nEV_Time\t=\t"<<config_tct->EV_Time();


    conf_file<<"\n\n[Parameters]";
    conf_file<<"\n#low-field mobility for electrons, cm2*V^-1*s^-1";
    conf_file<<"\nMu0_Electrons\t=\t"<<config_tct->mu0_els();
    conf_file<<"\n#low-field mobility for holes, cm2*V^-1*s^-1";
    conf_file<<"\nMu0_Holes\t=\t"<<config_tct->mu0_holes();
    conf_file<<"\n#saturation velocity cm/s";
    conf_file<<"\nSaturationVelocity\t=\t"<<config_tct->v_sat();
    conf_file<<"\n# amplifier amplification";
    conf_file<<"\nAmplification\t=\t"<<config_tct->ampl();
    conf_file<<"\n# factor between charge in sensor and photodiode due to light splitting: Nsensor/Ndiode";
    conf_file<<"\nLightSplitter\t=\t"<<config_tct->light_split();
    conf_file<<"\n# resistance of the sensor and diode output, Ohm";
    conf_file<<"\nResistanceSensor\t=\t"<<config_tct->R_sensor();
    conf_file<<"\nResistancePhotoDetector\t=\t"<<config_tct->R_diode();
    conf_file<<"\n# pohotodetector responce for certain wavelength, A/W";
    conf_file<<"\nResponcePhotoDetector\t=\t"<<config_tct->RespPhoto();
    conf_file<<"\n# electron-hole pair creation energy, eV";
    conf_file<<"\nEnergyPair\t=\t"<<config_tct->E_pair();

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
void base::print_run(bool start) {

    QDateTime datetime;
    if(start) {
        std::cout<<"<-------------- NEW RUN -------------->"<<std::endl;
        std::cout<<"Run Started At: "<<datetime.currentDateTime().toString().toStdString()<<std::endl;
    }
    else std::cout<<"Run Finished At: "<<datetime.currentDateTime().toString().toStdString()<<std::endl;

}

void base::on_tbrowser_clicked()
{
    if(browserProcess != NULL && browserProcess->isOpen()) browserProcess->kill();
    QString program = QDir::currentPath()+"/tbrowser";
    browserProcess = new QProcess(this);
    browserProcess->setWorkingDirectory(config_tct->OutFolder().c_str());
    browserProcess->setProcessChannelMode(QProcess::ForwardedChannels);
    browserProcess->start(program);
}
void base::kill_tbrowser() {
    if(browserProcess != NULL && browserProcess->isOpen()) browserProcess->kill();
}

void base::on_actionAbout_triggered()
{
    QMessageBox::about(this, tr("About, version 5.11.2015"),
            tr("<h2>TCT Data Analysis Framework. Graphical Version.</h2>"
               "<p>-> <b>Mykyta Haranko, 2015</b>"
               "<p>-> Oscilloscope data analysis by <b>Hendrik Jansen</b>"
               "<p>-> TCT Data files readout system by <b>particulars.si</b>"
               "<p><center><img src=\"../icons/knu.png\" width=\"87\"/> <img src=\"../icons/desy.png\" width=\"87\"/> <img src=\"../icons/particulars.png\" width=\"50\"/></center>"));
}
