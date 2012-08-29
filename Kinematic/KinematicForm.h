/********************************************************************************
** Form generated from reading UI file 'KinematicForm.ui'
**
** Created: Sun Apr 29 21:38:15 2012
**      by: Qt User Interface Compiler version 4.7.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_KINEMATICFORM_H
#define UI_KINEMATICFORM_H

#include <Qt3Support/Q3MimeSourceFactory>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>

QT_BEGIN_NAMESPACE

class Ui_KinematicForm
{
public:

    void setupUi(QDialog *KinematicForm)
    {
        if (KinematicForm->objectName().isEmpty())
            KinematicForm->setObjectName(QString::fromUtf8("KinematicForm"));
        KinematicForm->resize(772, 527);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(KinematicForm->sizePolicy().hasHeightForWidth());
        KinematicForm->setSizePolicy(sizePolicy);
        KinematicForm->setBaseSize(QSize(680, 600));

        retranslateUi(KinematicForm);

        QMetaObject::connectSlotsByName(KinematicForm);
    } // setupUi

    void retranslateUi(QDialog *KinematicForm)
    {
        KinematicForm->setWindowTitle(QApplication::translate("KinematicForm", "Kinematic editor", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class KinematicForm: public Ui_KinematicForm {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_KINEMATICFORM_H
