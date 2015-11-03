#ifndef FOLDERS_H
#define FOLDERS_H

#include <QDialog>

namespace Ui {
class Folders;
}

class Folders : public QDialog
{
    Q_OBJECT

public:
    explicit Folders(QWidget *parent = 0);
    Ui::Folders *ui;
    ~Folders();

private slots:
    void on_chdir_data_clicked();

    void on_chdir_output_clicked();

private:
};

#endif // FOLDERS_H
