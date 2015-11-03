#include "gui_folders.h"
#include "ui_form_folders.h"
#include "base.h"
#include "QFileDialog"
#include "QString"

Folders::Folders(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Folders)
{
    ui->setupUi(this);
}

Folders::~Folders()
{
    delete ui;
}

void Folders::on_chdir_data_clicked()
{
    QString cpart = "";
    cpart = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                              ui->cc_data->text(),
                                              QFileDialog::ShowDirsOnly
                                              | QFileDialog::DontResolveSymlinks);
    if(cpart=="") return;
    ui->cc_data->setText(cpart);
}

void Folders::on_chdir_output_clicked()
{
    QString cpart = "";
    cpart = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                              ui->cc_out->text(),
                                              QFileDialog::ShowDirsOnly
                                              | QFileDialog::DontResolveSymlinks);
    if(cpart=="") return;
    ui->cc_out->setText(cpart);
}
