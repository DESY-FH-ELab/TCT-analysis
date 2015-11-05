#ifndef GUI_CONSOLEOUTPUT_H
#define GUI_CONSOLEOUTPUT_H

#include <QProgressDialog>
#include <QLayout>
#include <QTextEdit>
#include <QPushButton>

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

        layout1 = new QGridLayout();
        layout1->addWidget(console);
        layout1->setContentsMargins(10,10,10,100);
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
    };
    ~ConsoleOutput() {
        delete console;
        delete layout1;
    }
    QTextEdit* Console() {return console;}

private:
    QTextEdit *console;
    QGridLayout *layout1;
private slots:
    void rename() {
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
        layout1->addWidget(console,0,0,1,4);
        layout1->addWidget(close,1,3,1,1);
        layout1->setContentsMargins(10,10,10,10);
        setLayout(layout1);

        setWindowTitle("Analysis Progress");
        setMinimumWidth(500);
        setMinimumHeight(400);
        setWindowModality(Qt::WindowModal);
        close->setEnabled(false);
        connect(close,SIGNAL(clicked()),this,SLOT(close()));

    };
    ~ConsoleOsc() {
        delete console;
        delete close;
        delete layout1;
    }
    QTextEdit* Console() {return console;}
    void SetButtonEnabled() {close->setEnabled(true);}

private:
    QTextEdit *console;
    QPushButton *close;
    QGridLayout *layout1;

};

}

#endif // GUI_CONSOLEOUTPUT_H
