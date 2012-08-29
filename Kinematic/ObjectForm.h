/********************************************************************************
** Form generated from reading UI file 'ObjectForm2.ui'
**
** Created: Tue May 1 18:46:41 2012
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OBJECTFORM2_H
#define UI_OBJECTFORM2_H

#include <Qt3Support/Q3Frame>
#include <Qt3Support/Q3ListBox>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionOpen;
    QAction *actionSave;
    QAction *actionSave_As;
    QAction *actionQuit;
    QAction *actionDraw_Colliding_Tris;
    QWidget *centralwidget;
    Q3ListBox *_lb;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuTools;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(188, 386);
        actionOpen = new QAction(MainWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        actionSave = new QAction(MainWindow);
        actionSave->setObjectName(QString::fromUtf8("actionSave"));
        actionSave_As = new QAction(MainWindow);
        actionSave_As->setObjectName(QString::fromUtf8("actionSave_As"));
        actionQuit = new QAction(MainWindow);
        actionQuit->setObjectName(QString::fromUtf8("actionQuit"));
        actionDraw_Colliding_Tris = new QAction(MainWindow);
        actionDraw_Colliding_Tris->setObjectName(QString::fromUtf8("actionDraw_Colliding_Tris"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        _lb = new Q3ListBox(centralwidget);
        _lb->setObjectName(QString::fromUtf8("_lb"));
        _lb->setGeometry(QRect(0, 0, 191, 361));
        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 188, 25));
        menubar->setMinimumSize(QSize(188, 25));
        menubar->setBaseSize(QSize(188, 25));
        menubar->setDefaultUp(false);
        menubar->setNativeMenuBar(false);
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuTools = new QMenu(menubar);
        menuTools->setObjectName(QString::fromUtf8("menuTools"));
        MainWindow->setMenuBar(menubar);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuTools->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addAction(actionSave);
        menuFile->addAction(actionSave_As);
        menuFile->addAction(actionQuit);
        menuTools->addAction(actionDraw_Colliding_Tris);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Objects", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        actionSave->setText(QApplication::translate("MainWindow", "Save", 0, QApplication::UnicodeUTF8));
        actionSave_As->setText(QApplication::translate("MainWindow", "Save As", 0, QApplication::UnicodeUTF8));
        actionQuit->setText(QApplication::translate("MainWindow", "Exit", 0, QApplication::UnicodeUTF8));
        actionDraw_Colliding_Tris->setText(QApplication::translate("MainWindow", "Draw Colliding Triangles", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuTools->setTitle(QApplication::translate("MainWindow", "Tools", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OBJECTFORM2_H
