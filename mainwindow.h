#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include"vypocetni.hpp"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_pushButton_2_clicked();

    void on_pushButton_clicked();

    void on_lineEdit_returnPressed();

    void on_lineEdit_editingFinished();

    void on_lineEdit_3_editingFinished();

    void on_lineEdit_5_editingFinished();

    void on_lineEdit_7_editingFinished();

    void on_lineEdit_9_cursorPositionChanged(int arg1, int arg2);

    void on_lineEdit_9_editingFinished();

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

    void on_pushButton_6_clicked();

    void on_pushButton_7_clicked();

    void on_pushButton_8_clicked();

    void on_pushButton_9_clicked();

    void on_pushButton_10_clicked();

private:
    Ui::MainWindow *ui;
    void readData();
};

#endif // MAINWINDOW_H
