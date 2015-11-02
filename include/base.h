#ifndef BASE_H
#define BASE_H

#include <QMainWindow>
#include "stdint.h"
#include "util.h"
#include "tct_config.h"
#include "sample.h"
namespace Ui {
class base;
}

class base : public QMainWindow
{
    Q_OBJECT

public:
    explicit base(QWidget *parent = 0);
    ~base();

private slots:

    void on_chdir_data_clicked();

    void on_chdir_output_clicked();

    void on_cc_ch_1_toggled(bool checked);

    void on_cc_ch_2_toggled(bool checked);

    void on_cc_ch_ph_toggled(bool checked);

    void on_cc_ch_trig_toggled(bool checked);

    void on_buttonGroup_buttonClicked(int index);

    void on_window_sample_clicked();

    void on_windows_parameters_clicked();

    void on_chdir_sensor_clicked();

    void on_start_clicked();

    void on_actionChange_config_triggered();

    void on_actionSave_config_triggered();

private:
    Ui::base *ui;
    std::string _DefConfigName = "def";
    TCT::tct_config *config_analysis = NULL;
    TCT::sample *config_sample = NULL;
    TCT::mode_selector *config_mode = NULL;
    void read_config(const char*);
    void fill_config();
    void tovariables_config();
    void error(const char*);
    void not_run();
};

#endif // BASE_H
