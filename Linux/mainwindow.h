#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <vector>
#include <QThread>
#include <qelapsedtimer.h>
#include <QWaitCondition>
#include <time.h>
#include <QVector>

namespace Ui {
class MainWindow;
}


//typedef and structures

// memory element where store info related to the single person in database
typedef struct memory_element_s
{
    std::string name_t;
    quint8 qgramsNum;              //Number of Q-grams name + surname
    quint8 age;
    double latitude;
    double longitude;
    quint32 attributeClass;
} element;

// store in the calcs results, include the corresponded element ID
typedef struct element_result_s
{
    quint64 element_id;
    double euclideanDistance;
    double geoDistance;
    double stringDistance;
    double ageDistance;
}results;

//base element for tree node
typedef struct treeNode_s
{
   quint64 indexPos;
   treeNode_s* majorBranch;
   treeNode_s* minorBranch;
   quint32 classNumber;
   quint32 numElements;
   treeNode_s* parent;
   qint8 param;
   bool searched;
}treeNode;


//Classes declaration

class MainWindow;
//threads
class CalcStringDistThread : public QThread
{
    Q_OBJECT
protected:
    void run();

public:
    element* db;
    element test;
    int* qgIndex;
    quint32 qgramsNum;
    results* res;
    quint64 first;
    quint64 last;
    MainWindow* parent;
    int* strDistances = NULL;

public:
    void CalcQGramsDist();
    double CalcSingleStringDist(quint64 Index);
    void resetCalcs();
};

class CalcGeoDistThread : public QThread
{
    Q_OBJECT
protected:
    void run()
    {
        for (;first < last; first++)
             res[first].geoDistance = CalculateLocalizationDistance(&test,&db[res[first].element_id]);
    }
public:
    element* db;
    element test;
    results* res;
    quint64 first;
    quint64 last;

public:
    double CalculateLocalizationDistance(element* peopleA, element* peopleB);
};


class CalcAgeDistThread : public QThread
{
    Q_OBJECT
protected:
    void run()
    {
        for (;first < last; first++)
             res[first].ageDistance = CalculateAgeDistance(&test,&db[res[first].element_id]);


    }
public:
    element* db;
    element test;
    results* res;
    quint64 first;
    quint64 last;

public:
    double CalculateAgeDistance(element* peopleA, element* peopleB);
};


class CalcEuclideanDistThread : public QThread
{
    Q_OBJECT
protected:
    void run()
    {
        for (;_first < _last; _first++)
            CalculateEuclideanDistance(&_res[_first]);
    }
private:
    results* _res;
    quint64 _first;
    quint64 _last;
    double _kAge;
    double _kString;
    double _kGeo;
public:
    void SetParameter(results* resultsArray,
                      quint64 firstElem,
                      quint64 lastElem,
                      double KAge,
                      double KGeo,
                      double KString);
    void SetLastElement(quint64 last){_last = last;}
    void CalculateEuclideanDistance(results* res);
};

class Q_GramsProcessingThread : public QThread
{
    Q_OBJECT
protected:
    void run()
    {
        Q_GramsProcessing();
    }
private:
    QMutex* _mutex;
    quint64 _first;
    quint64 _last;
    MainWindow* _parent;
public:
    void SetParams(QMutex* mutex,
                   MainWindow* parent,
                   quint64 first,
                   quint64 last);
    void SetLastElement(quint64 last){_last = last;}
    void SetFirstElement(quint64 first){_first = first;}
    void Q_GramsProcessing();
};



class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    void on_horizontalSlider_valueChanged(int value);
    void on_pbCreateDB_clicked();
    void on_pbSearch_clicked();
    void on_pbFastSearch_clicked();
    void on_pbClear_clicked();
    void on_leKString_editingFinished();
    void on_leKGeo_editingFinished();
    void on_leKAge_editingFinished();
    void on_leClusterThr_editingFinished();

private:
    void closeEvent(QCloseEvent *event);
    void SaveSettings();
    void LoadSettings();

    Ui::MainWindow *ui;
    //thread dist management
    QElapsedTimer  taskTimer;

    CalcAgeDistThread ageThread;
    CalcEuclideanDistThread euclidThread;
    CalcGeoDistThread geoThread;
    CalcStringDistThread strThread;
    double* stringDistances;
    results* distResults;
    double kAge;
    double kString;
    double kGeo;
    quint64 maxClusterSize;
    QString m_sSettingsFile;

public:
    element* people;
    quint64 peopleCount;
    quint32 totClassCreated;
    treeNode* bestMatch;
    treeNode* root1;
    element userSearch;


    void ThreadHandler();
    void CreateMemoryBase();
    void ClusterizeMemBase(treeNode* branch);
    void ClusterizeMemBase2(treeNode* branch);
    void Q_GramsProcessing();
    void QuickSort(results* rs,  quint64 numResults, quint64 theshold);
    void QuickSearch(results* distances, quint32 resElements);
    void StraightTreeCheck(treeNode* currNode);
    void LinearSearch(quint64 numElements, results* distance);
    void DeleteNode(treeNode* currentNode);
    void Clusterize();
    void QGramsProcessing();
    void StartTreeSort();
    void TreeSort(treeNode* branch, quint64 *startPoint, element* peopleCopy);
    void ReverseTreeSearch(results* distances, quint32 resElements, treeNode* currNode,quint64* position);
    void NodeResetSearch(treeNode* currentNode);

    QVector<QByteArray> q_gramsIdx_v;
    QVector<QVector<quint64> > invertedTable;
};

#endif // MAINWINDOW_H
