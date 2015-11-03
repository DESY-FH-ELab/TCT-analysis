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

    void on_buttonGroup_mode_buttonClicked(int index);

    void on_window_folders_clicked();

    void on_window_sample_clicked();

    void on_window_parameters_clicked();

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
