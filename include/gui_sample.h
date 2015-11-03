#ifndef GUI_SAMPLE_H
#define GUI_SAMPLE_H

#include <sample.h>
#include <QDialog>

namespace Ui {
class Sample;
}

class Sample : public QDialog
{
    Q_OBJECT

public:
    explicit Sample(QWidget *parent = 0);
    Ui::Sample *ui;
    ~Sample();

private slots:
    void on_chdir_sensor_clicked();

private:
    TCT::sample *temp_config;

};

#endif // GUI_SAMPLE_H
