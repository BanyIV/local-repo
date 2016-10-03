#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QMessageBox>
#include <QFileDialog>

#include <qmath.h>
#include<iostream>

#include "QWT-6.1.3/include/qwt_plot.h"
#include "QWT-6.1.3/include/qwt_plot_curve.h"
#include "QWT-6.1.3/include/qwt_plot_grid.h"
#include "QWT-6.1.3/include/qwt_symbol.h"
#include "QWT-6.1.3/include/qwt_legend.h"
#include "QWT-6.1.3/include/qwt_plot_spectrogram.h"
#include "QWT-6.1.3/include/qwt_color_map.h"

QwtPlot *rez, *mapa;

class SpectrogramData: public QwtRasterData
{
public:
    SpectrogramData(double xmin = 0.0, double xmax = 1.1, double ymin = 0.0, double ymax = 1.1, double zmin = 0.0, double zmax = 1.1)
    {
        setInterval( Qt::XAxis, QwtInterval( xmin, xmax ) );
        setInterval( Qt::YAxis, QwtInterval( ymin, ymax ) );
        setInterval( Qt::ZAxis, QwtInterval( zmin, zmax ) );
    }

    virtual double value( double x, double y ) const
    {
        return hydraulic_head(x,y);
    }
};

void printout(double x[], int dim);
void ExportPlot(QwtPlot *arg, const char *filename);

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->lineEdit->setValidator(new QDoubleValidator(this));
    ui->lineEdit_2->setValidator(new QDoubleValidator(this));
    ui->lineEdit_3->setValidator(new QDoubleValidator(this));
    ui->lineEdit_4->setValidator(new QDoubleValidator(this));
    ui->lineEdit_5->setValidator(new QDoubleValidator(this));
    ui->lineEdit_6->setValidator(new QDoubleValidator(this));
    ui->lineEdit_7->setValidator(new QDoubleValidator(this));
    ui->lineEdit_8->setValidator(new QDoubleValidator(this));
    ui->lineEdit_9->setValidator(new QDoubleValidator(this));
    ui->lineEdit_10->setValidator(new QDoubleValidator(this));
    ui->lineEdit_11->setValidator(new QDoubleValidator(this));
    ui->lineEdit_12->setValidator(new QDoubleValidator(this));
    ui->lineEdit_13->setValidator(new QDoubleValidator(this));
    ui->lineEdit_14->setValidator(new QDoubleValidator(this));
    ui->lineEdit_15->setValidator(new QDoubleValidator(this));
    ui->lineEdit_16->setValidator(new QDoubleValidator(this));
    ui->lineEdit_17->setValidator(new QDoubleValidator(this));
    ui->lineEdit_18->setValidator(new QDoubleValidator(this));
    ui->lineEdit_19->setValidator(new QDoubleValidator(this));
    ui->lineEdit_20->setValidator(new QDoubleValidator(this));
    ui->lineEdit_21->setValidator(new QDoubleValidator(this));
}

MainWindow::~MainWindow()
{
    delete rez;
    delete mapa;
    delete ui;
}


void MainWindow::on_pushButton_2_clicked()
{
QMessageBox::critical(NULL,"problem",QString::number(N));
}

void MainWindow::on_pushButton_clicked()
{
   // readData();
testovaci_input();
    QGraphicsScene * schema = new QGraphicsScene();

    for (int i=0;i<6;i++)
    {
    QGraphicsLineItem * hloubka = new QGraphicsLineItem();
    hloubka->setLine(0,z[i],100,z[i]);
    schema->addItem(hloubka);
}
    ui->graphicsView->setScene(schema);
    ui->graphicsView->fitInView(schema->sceneRect());

    rez = new QwtPlot(ui->FrameHH);
    rez->setTitle("Řez hydraulické výšky v linii vrtů");
    rez->setCanvasBackground(Qt::white);
    rez->setAxisScale(QwtPlot::yLeft, z[0], z[N]);
    rez->setAxisScale(QwtPlot::xBottom, -L/2 , L*1.5);
    rez->insertLegend(new QwtLegend);

    QwtPlotCurve *vrstva = new QwtPlotCurve("vrstva");
    //krivka->setSamples( x, rezh, xDim);
    krivka->setPen( Qt::red, 1 );
    krivka->attach(rez);

    QwtPlotGrid *mrizka = new QwtPlotGrid();
    mrizka->attach(rez);

    rez->resize(ui->FrameHH->width(), ui->FrameHH->height());
    rez->show();
    
    
}

void MainWindow::on_lineEdit_returnPressed()
{
 //   ui->lineEdit->
}

void MainWindow::readData()
{
            // nacist data z gui do poli

            // hloubky rozhrani

                    z[0] = 0;
                    z[1] = ui->lineEdit->text().toDouble();
                    z[2] = ui->lineEdit_3->text().toDouble();
                    z[3] = ui->lineEdit_5->text().toDouble();
                    z[4] = ui->lineEdit_7->text().toDouble();
                    z[5] = ui->lineEdit_9->text().toDouble();

                    d[0] = z[1];
                    d[1] = z[2] - z[1];
                    d[2] = z[3] - z[2];
                    d[3] = z[4] - z[3];
                    d[4] = z[5] - z[4];


            // hydraulicke vodivosti
            K[0] = 0;
            K[1] = 0;
            K[2] = 0;
            K[3] = 0;
            K[4] = 0;

            // vydatnosti
            Q[0] = 0;
            Q[1] = 0;

            // polomery studni
            r[0] = 0;
            r[1] = 0;

            // polomery dosahu (odhad, nebo vnucena hodnota
            R[0] = 0;//ui->lineEdit_20->text().toDouble();
            R[1] = 0;//ui->lineEdit_21->text().toDouble();;

            // snizeni
            s[0] = 0;
            s[1] = 0;

            // vzdalenost mezi studnami
            L = 0;

            // puvodni hydraulicka vyska (pred cerpanim - nebo hloubka hladiny podzemni vody? bylo by praktictejsi
            H = 0;

            // veci odvozene:

            // transmisivita kolektoru
            T = 0;
            for(int i=0; i < N; i++)
            {
                    T+= K[i]*d[i];
            }

            // tohle je moje interni vec, tuhle promennou budu kontrolovat, jestli je true, kdybych si nebyl jist, jestli uz je vsechno nacteno
            // tenhle kus kodu pripravi z-souradnice vrstevnich rozhrani - obrati osu, polozi z=0 na podlozi
            double Z_base = z[N];

            for(int i = 0; i < N+1; i++)
            {
                z[i] = Z_base - z[i];
                //cerr << "z[" << i << "] = " << z[i] << endl;
            }

            for(int i = 0; i < (N+1)/2; i++)
            {
                double pom;
                pom = z[i];
                z[i] = z[N-i];
                z[N-i] = pom;
                //cerr << "z[" << i << "] = " << z[i] << endl;
            }

            H = Z_base - H;
}

void MainWindow::on_lineEdit_editingFinished()
{
    if( !ui->lineEdit->text().isEmpty() && !ui->lineEdit_2->text().isEmpty())
    {
        ui->lineEdit_3->setEnabled(true);
        ui->lineEdit_4->setEnabled(true);
        if(N==0)
            N++;
    }
    else
    {
        ui->lineEdit_3->setEnabled(false);
        ui->lineEdit_4->setEnabled(false);
        if(N==1)
            N--;
    }
}

void MainWindow::on_lineEdit_3_editingFinished()
{
    if( !ui->lineEdit_3->text().isEmpty() && !ui->lineEdit_4->text().isEmpty())
    {
        ui->lineEdit_5->setEnabled(true);
        ui->lineEdit_6->setEnabled(true);
        if(N==1)
            N++;
    }
    else
    {
        ui->lineEdit_5->setEnabled(false);
        ui->lineEdit_6->setEnabled(false);
        if(N==2)
            N--;
    }
}

void MainWindow::on_lineEdit_5_editingFinished()
{
    if( !ui->lineEdit_5->text().isEmpty() && !ui->lineEdit_6->text().isEmpty())
    {
        ui->lineEdit_7->setEnabled(true);
        ui->lineEdit_8->setEnabled(true);
        if(N==2)
            N++;
    }
    else
    {
        ui->lineEdit_7->setEnabled(false);
        ui->lineEdit_8->setEnabled(false);
        if(N==3)
            N--;
    }
}

void MainWindow::on_lineEdit_7_editingFinished()
{
    if( !ui->lineEdit_7->text().isEmpty() && !ui->lineEdit_8->text().isEmpty())
    {
        ui->lineEdit_9->setEnabled(true);
        ui->lineEdit_10->setEnabled(true);
        if(N==3)
            N++;
    }
    else
    {
        ui->lineEdit_9->setEnabled(false);
        ui->lineEdit_10->setEnabled(false);
        if(N==4)
            N--;
    }
}

void MainWindow::on_lineEdit_9_editingFinished()
{
    if( !ui->lineEdit_9->text().isEmpty() && !ui->lineEdit_10->text().isEmpty())
    {
        N=5;
    }
    else
    {
        N=4;
    }
}

void MainWindow::on_pushButton_3_clicked() // chceme spocitat a zobrazit hydraulickou vysku
{
    testovaci_input();

    if(rez != NULL)
        delete rez;

    // rozmery pole: x: 2*L, y: whatever
    int xDim = int(2*L);
    double x[xDim], rezh[xDim];

    // mrizka souradnic a funkcni hodnoty
    for(int i = 0; i < xDim; i++) {
        x[i] = i - L/2;
        if(isnan(x[i]))
        {
            cerr << "i    = " << i << endl;
            cerr << "xDim = " << xDim << endl;
            exit(0);
        }
        rezh[i] = hydraulic_head(x[i],0.0);
    }

    // a ted to cele vynest:
    rez = new QwtPlot(ui->FrameHH);
    rez->setTitle("Řez hydraulické výšky v linii vrtů");
    rez->setCanvasBackground(Qt::white);
    rez->setAxisScale(QwtPlot::yLeft, z[0], z[N]);
    rez->setAxisScale(QwtPlot::xBottom, -L/2 , L*1.5);
    rez->insertLegend(new QwtLegend);

    QwtPlotCurve *krivka = new QwtPlotCurve("hydraulická výška");
    krivka->setSamples( x, rezh, xDim);
    krivka->setPen( Qt::blue, 1 );
    krivka->attach(rez);

    QwtPlotGrid *mrizka = new QwtPlotGrid();
    mrizka->attach(rez);
    
    rez->resize(ui->FrameHH->width(), ui->FrameHH->height());
    rez->show();

}

void MainWindow::on_lineEdit_9_cursorPositionChanged(int arg1, int arg2)
{
    int k = arg1+arg2; // just to shut up that compile time warning 'bout unused params
    k++;
    // what?
}

void printout(double x[], int dim) // vypis pole pro manualni debug
{
    for(int i = 0; i < dim; i++)
        cerr << "["<< i << "] = " << x[i] << endl;
}

void MainWindow::on_pushButton_4_clicked() // TAB: REZ HYDRAULICKE VYSKY, zmena mezi na svisle ose
{

    double ymin = ui->hhymin->text().toDouble();
    double ymax = ui->hhymax->text().toDouble();

    if(ui->hhymin->text().isEmpty())
        ymin = z[0];

    if(ui->hhymax->text().isEmpty())
        ymax = z[N];

    rez->setAxisScale(QwtPlot::yLeft, ymin, ymax);
    rez->replot();
}

void MainWindow::on_pushButton_5_clicked() // TAB: REZ HYDRAUL. V., export
{
    ExportPlot(rez, "rez_hydraul_vyska.png");
}

void ExportPlot(QwtPlot *arg, const char *filename)
{
    //QPixmap qPix = QPixmap::grabWidget(arg);
    QPixmap qPix = arg->grab();

    if(qPix.isNull()){
        qDebug("Export: Failed to capture the plot for saving.");
        return;
    }

    //QFileDialog *saving = new QFileDialog(NULL,"Kepsn", "c:/");


    //QFileDialog(ui,"Export řezu hydraulické výšky", new QString("C:/"), new QString(".png"));
    qPix.save(filename,"PNG",85);
}

/* jen kus kodu - ukazkovy zpusob, jak vynest jednoduchy graf:
    QwtPlot *plot = new QwtPlot(ui->ramecek);
    plot->setTitle( "Plot Demo" );
    plot->setCanvasBackground( Qt::white );
    plot->setAxisScale( QwtPlot::yLeft, 0.0, 10.0 );
    plot->insertLegend( new QwtLegend() );

    QwtPlotGrid *grid = new QwtPlotGrid();
    grid->attach( plot );

    QwtPlotCurve *curve = new QwtPlotCurve();
    curve->setTitle( "Some Points" );
    curve->setPen( Qt::blue, 4 ),
    curve->setRenderHint( QwtPlotItem::RenderAntialiased, true );

    QwtSymbol *symbol = new QwtSymbol( QwtSymbol::Ellipse,
        QBrush( Qt::yellow ), QPen( Qt::red, 2 ), QSize( 4, 8 ) );
    curve->setSymbol( symbol );

    QPolygonF points;
    points << QPointF( 0.0, 4.4 ) << QPointF( 1.0, 3.0 )
        << QPointF( 2.0, 4.5 ) << QPointF( 3.0, 6.8 )
        << QPointF( 4.0, 7.9 ) << QPointF( 5.0, 7.1 );
    curve->setSamples( points );

    curve->attach( plot );

    plot->resize( ui->ramecek->width(), ui->ramecek->height() );
    plot->show();
*/

void MainWindow::on_pushButton_6_clicked() // TAB: MAPA: vypocet
{
    testovaci_input();

    if(mapa !=NULL)
        delete mapa;

    mapa = new QwtPlot(ui->widgetMapa);

    double xmin = -.5*L; //interval delky 2L
    double xmax = 1.5*L;

    double ymin =-1.0*ui->widgetMapa->width()/ui->widgetMapa->height() * L;
    double ymax = 1.0*ui->widgetMapa->width()/ui->widgetMapa->height() * L;

    //cerr << ymin << endl;
    //cerr << ymax << endl;

    //polointeligentni urceni mezi z:
    double arr[3] = { H, hydraulic_head(r[0],0), hydraulic_head(L+r[1],0) };

    for(int pruchod = 0; pruchod < 3; pruchod++)
        for(int k = 0; k < 2; k++)
            if(arr[k] > arr[k+1])
            {
                double p = arr[k];
                arr[k] = arr[k+1];
                arr[k+1] = p;
            }

    double zmin = arr[0];
    double zmax = arr[2];

    xmin = -25;
    xmax = 75;
    ymin = -50;
    ymax = 50;


    SpectrogramData *obsah = new SpectrogramData();
    obsah->setInterval( Qt::XAxis, QwtInterval( xmin, xmax ) );
    obsah->setInterval( Qt::YAxis, QwtInterval( ymin, ymax ) );
    obsah->setInterval( Qt::ZAxis, QwtInterval( zmin, zmax ) );

    QwtPlotSpectrogram *graf = new QwtPlotSpectrogram("Hydraulická výška");
    graf->setRenderThreadCount(4);
    graf->setCachePolicy(QwtPlotRasterItem::PaintCache);

    //kontury
    QList<double> contourLevels;
    for ( double level = zmin; level < zmax + .1; level += (zmax-zmin)/10 )
        contourLevels += level;
    graf->setContourLevels( contourLevels );

    graf->setData(obsah);
    graf->attach(mapa);

    //barvy
    QwtLinearColorMap *barvy = new QwtLinearColorMap(Qt::blue, Qt::cyan, QwtColorMap::RGB);
    barvy->addColorStop(H,Qt::green);
    graf->setColorMap(barvy);
    graf->setDisplayMode( QwtPlotSpectrogram::ContourMode);

    //_
    mapa->setAxisScale( QwtPlot::yRight, zmin, zmax);
    mapa->enableAxis(QwtPlot::yRight);

    mapa->setAxisScale(QwtPlot::xBottom, xmin, xmax);
    mapa->setAxisScale(QwtPlot::yLeft, ymin, ymax);
    //mapa->setAxisAutoScale(QwtPlot::yLeft);
    //mapa->setAxisAutoScale(QwtPlot::xBottom);

    mapa->resize(ui->widgetMapa->width(), ui->widgetMapa->height());
    mapa->show();


}

void MainWindow::on_pushButton_7_clicked() // TAB: MAPA: export
{
    ExportPlot(mapa, "mapa.png");
}

void MainWindow::on_pushButton_8_clicked() //dopočítáme snížení
{
  //readData();
testovaci_input();
  if ( ui->lineEdit_11->text().isEmpty() || ui->lineEdit_14->text().isEmpty())
  {
      QMessageBox::warning(NULL,"problem","Zadejte vydatnost!");
      return;
  }

    s[0]=wellDrawdown(0);
    s[1]=wellDrawdown(1);
  ui->lineEdit_13->setText(QString::number(s[0]));
  ui->lineEdit_16->setText(QString::number(s[1]));

}

void MainWindow::on_pushButton_9_clicked() //dopočítáme vydatnost
{
 //readData();
    testovaci_input();

  if (ui->lineEdit_20->text().isEmpty() || ui->lineEdit_21->text().isEmpty())
  {
      QMessageBox::warning(NULL,"problem","Zadejte oba dosahy!");
              return;
  }

  if ( ui->lineEdit_11->text().isEmpty() && ui->lineEdit_14->text().isEmpty())
  {
      QMessageBox::warning(NULL,"problem","Zadejte jednu vydatnost!");
      return;
  }

  if (ui->lineEdit_14->text().isEmpty() && ui->lineEdit_16->text().isEmpty())
  {
      QMessageBox::warning(NULL,"problem","Zadejte snizeni ve studni,ve ktere chcete znad vydatnost!");
      return;
  }

  if (ui->lineEdit_11->text().isEmpty() && ui->lineEdit_13->text().isEmpty())
  {
      QMessageBox::warning(NULL,"problem","Zadejte snizeni ve studni,ve ktere chcete znad vydatnost!");
      return;
  }

  if (ui->lineEdit_11->text().isEmpty())
  {
      Q[0]=wellYield(0);
      ui->lineEdit_11->setText(QString::number(Q[0]*1000));
  }

  if (ui->lineEdit_14->text().isEmpty())
  {
      Q[1]=wellYield(1);
      ui->lineEdit_14->setText(QString::number(Q[1]*1000));
  }
}

void MainWindow::on_pushButton_10_clicked() // odhad dosahu depr kuzele
{
    if (ui->lineEdit_13->text().isEmpty() || ui->lineEdit_16->text().isEmpty())
    {
        QMessageBox::warning(NULL,"problem","Zadejte snizeni!");
                return;
    }
   // readData();
    testovaci_input();
    R[0]=fabs(575*s[0]*sqrt(T));
    ui->lineEdit_20->setText(QString::number(R[0]));
    R[1]=fabs(575*s[1]*sqrt(T));
    ui->lineEdit_21->setText(QString::number(R[1]));
}
