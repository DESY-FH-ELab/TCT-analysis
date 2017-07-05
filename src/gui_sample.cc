/**
 * \file
 * \brief Implementation of Ui::Sample methods
 */

#include "gui_sample.h"
#include "ui_form_sample.h"
#include "QFileDialog"
#include "util.h"
#include "base.h"

Sample::Sample(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Sample)
{
    ui->setupUi(this);
}

Sample::~Sample()
{
    delete ui;
}
void Sample::on_chdir_sensor_clicked()
{
    QString cpart = "";
    cpart = QFileDialog::getOpenFileName(this,
                                         tr("Open Sample Card"), ui->folder->text(), tr("Text Files (*.txt)"));

    if(cpart=="") return;
    ui->cc_sensor->setText(cpart);


    TCT::util sample_card;
    std::ifstream sample_file(cpart.toStdString().c_str());
    sample_card.parse(sample_file);
    temp_config = new TCT::sample(sample_card.ID_val());

    ui->folder->setText(temp_config->Folder().c_str());
    ui->id->setText(temp_config->SampleID().c_str());
    ui->thickness->setText(QString::number(temp_config->Thickness()));

    delete temp_config;
}
