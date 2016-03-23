/**
 * \file
 * \brief Definition of the Ui::base class.
 * \details This class is the Main GUI File.
 */

#ifndef BASE_H
#define BASE_H

#include <QMainWindow>
#include "stdint.h"
#include "util.h"
#include "tct_config.h"
#include "sample.h"
#include "analysis.h"
#include "gui_consoleoutput.h"
#include "QProcess"

namespace Ui {
class base;
}

class base : public QMainWindow
{
    Q_OBJECT

public:
    explicit base(QWidget *parent = 0);
    /*std::vector<QWidget*> top_widgets;
    std::vector<QWidget*> edge_widgets;
    std::vector<QWidget*> bottom_widgets;*/
    std::map<QWidget*,int> top_widgets;
    std::map<QWidget*,int> edge_widgets;
    std::map<QWidget*,int> bottom_widgets;
    ~base();

private slots:

    void on_buttonGroup_mode_buttonClicked(int index);

    void on_window_folders_clicked();

    void on_window_sample_clicked();

    void on_window_parameters_clicked();

    void on_start_clicked();

    void on_actionChange_config_triggered();

    void on_actionSave_config_triggered();

    void deleteprogress();

    void deleteprogress_osc();

    void on_tbrowser_clicked();

    void kill_tbrowser();

    void on_actionAbout_triggered();

private:
    Ui::base *ui;
    std::string _DefConfigName;
    TCT::tct_config *config_tct;
    TCT::analysis *config_analysis;
    TCT::sample *config_sample;
    TCT::mode_selector *config_mode;
    Ui::ConsoleOutput *progress;
    Ui::ConsoleOsc *progress_osc;
    QProcess *browserProcess;
    void read_config(const char*);
    void fill_config();
    void tovariables_config();
    void error(const char*);
    void print_run(bool start);
    void start_tct();
    void start_osc();
    void FillBoxes();

};

#endif // BASE_H
