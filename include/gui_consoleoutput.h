/**
 * \file
 * \brief Definition of the Ui::ConsoleOutput and Ui::ConsoleOsc classes. Responsible for progress bars in GUI.
 */

#ifndef GUI_CONSOLEOUTPUT_H
#define GUI_CONSOLEOUTPUT_H

#include <QProgressDialog>
#include <QLayout>
#include <QTextEdit>
#include <QPushButton>
#include <QCommandLinkButton>

namespace Ui {

class ConsoleOutput : public QProgressDialog
{
    Q_OBJECT
public:
    ConsoleOutput(int max, QWidget* parent=0):
        QProgressDialog(parent){
        setMaximum(max);
        setCancelButtonText("Abort");

        console = new QTextEdit();
        tbrowser = new QPushButton("Run TBrowser");
        //tbrowser->setGeometry(10,300,120,31);
        tbrowser->setMaximumHeight(23);
        tbrowser->setEnabled(false);

        layout1 = new QGridLayout();
        layout1->addWidget(console,0,0,9,4);
        layout1->addWidget(tbrowser,9,0,1,1);
        layout1->setVerticalSpacing(50);
        //layout1->setContentsMargins(10,10,10,100);
        setLayout(layout1);

        setWindowTitle("Analysis Progress");
        //setFixedWidth(500);
        //setFixedHeight(400);
        setMinimumWidth(600);
        setMinimumHeight(400);
        setWindowModality(Qt::WindowModal);
        setAutoClose(false);
        setAutoReset(true);
        setMinimumDuration(0);
        connect(this,SIGNAL(finished(int)),this,SLOT(rename()));
        connect(tbrowser,SIGNAL(clicked()),this,SLOT(on_tbrowser_clicked()));
    };
    ~ConsoleOutput() {
        delete console;
        delete tbrowser;
        delete layout1;
    }
    QTextEdit* Console() {return console;}
    void finished1(int result) {finished(result);}

signals:
    void OpenTBrowser();

private:
    QTextEdit *console;
    QPushButton *tbrowser;
    QGridLayout *layout1;
private slots:
    void on_tbrowser_clicked() {
        emit OpenTBrowser();
        this->reject();
    }
    void rename() {
        tbrowser->setEnabled(true);
        setCancelButtonText("Close");
    }

};

class ConsoleOsc : public QDialog
{
    Q_OBJECT
public:
    ConsoleOsc(QWidget* parent=0):
        QDialog(parent){

        layout1 = new QGridLayout();
        console = new QTextEdit();
        close = new QPushButton("Ok");
        tbrowser = new QCommandLinkButton("Run TBrowser");
        layout1->addWidget(console,0,0,1,4);
        layout1->addWidget(tbrowser,1,0,1,1);
        layout1->addWidget(close,1,3,1,1);
        layout1->setContentsMargins(10,10,10,10);
        setLayout(layout1);

        setWindowTitle("Analysis Progress");
        setMinimumWidth(500);
        setMinimumHeight(400);
        setWindowModality(Qt::WindowModal);
        close->setEnabled(false);
        connect(close,SIGNAL(clicked()),this,SLOT(close()));
        connect(tbrowser,SIGNAL(clicked()),this,SLOT(on_tbrowser_clicked()));

    };
    ~ConsoleOsc() {
        delete console;
        delete close;
        delete tbrowser;
        delete layout1;
    }
    QTextEdit* Console() {return console;}
    void SetButtonEnabled() {close->setEnabled(true);}

signals:
    void OpenTBrowser();

private:
    QTextEdit *console;
    QPushButton *close;
    QCommandLinkButton *tbrowser;
    QGridLayout *layout1;
private slots:
    void on_tbrowser_clicked() {
        emit OpenTBrowser();
        close->click();
    }

};

}

#endif // GUI_CONSOLEOUTPUT_H
