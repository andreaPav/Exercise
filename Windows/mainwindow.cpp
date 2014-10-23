#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>
#include <QStringList>
#include <QFile>
#include <QTextStream>
#include <qsettings.h>
#include <qmutex.h>


#define K_STRING_WEIGHT 1                   //weight of the String distance - default
#define K_AGE_WEIGHT    1.2                 //weight of the Age distance - default
#define K_GEO_WEIGHT    0.04                //weight of the Geographical distance - default
#define CLUSTER_MAX_ELEMENTS 1000           //max cluster size
#define ELEMENT_NUMBER    1000000           //element of people multiplier (x slider value)
#define NUMBER_OF_CORE  4                   //number of core of the pc used
#define K_NEARER_NEIGH  10                  //number of nearest neighbour requested during the search
#define NEARER_PERCENTAGE 0.95  //(90%)     //precision requested (parameter not really valid..)



/**************************/
//Thread methods
/**************************/
/* fn: run()
 * brief: calculate the distance between a test string and an array of strings putting the results in another results array
 *         0 means identical strings, 100 + qgrams number difference means maximum distance
 */
void CalcStringDistThread::run()
{
    if(strDistances == NULL)
        CalcQGramsDist();
    //int threshold = 10;
    //int minLevelArray[20]={0};
    //int minLevel = 0;
    for (quint64 count= first; count < last; count++)
    {
        if(strDistances[res[count].element_id] > 0/*minLevel*/)
        {
            if(qgramsNum < db[res[count].element_id].qgramsNum)
                res[count].stringDistance = 100.0 - ((strDistances[res[count].element_id]*100.0) /(double)qgramsNum) + ((double)(db[res[count].element_id].qgramsNum - qgramsNum) );
            else
                res[count].stringDistance = 100.0 - ((strDistances[res[count].element_id]*100) /(double)db[res[count].element_id].qgramsNum) + ((double)(qgramsNum - db[res[count].element_id].qgramsNum));

            if(res[count].stringDistance < 0)
                res[count].stringDistance = 0;

            /*
            if(strDistances[res[count].element_id] >= 20)
            {
                if(minLevelArray[19]++ > threshold)
                    minLevel = 17;
            }
            else
            {
                if(minLevelArray[strDistances[res[count].element_id]]++ > threshold)
                    minLevel = strDistances[res[count].element_id]-1;
            }
            */
        }
        else
            res[count].stringDistance = 100 + (double)abs((qgramsNum - db[res[count].element_id].qgramsNum));
    }
}

//fn:CalcStringDist
//brief: this function calculate the number of q-grams identical (to the test string q-grams) for every string in the db
void CalcStringDistThread::CalcQGramsDist()
{
    //exit reason - not fair...
    if(test.name_t.length()<= 1)
    {
        //too short.. exit
        for (quint64 count= first; count < last; count++)
            res[count].stringDistance=100;
    }
    else
    {
        QString name = QString::fromUtf8(test.name_t.c_str());
        QByteArray *nameQgrams;
        qgramsNum = (test.name_t.length())-1;
        int qgIndex[qgramsNum];
        strDistances = new int[parent->peopleCount];
        std::fill_n(strDistances, parent->peopleCount, 0);
        nameQgrams = new QByteArray[qgramsNum];
        quint32 tempQgramsNum = 0;
        //split the string in q-grams (avoid qgrams that include spaces)
        for (quint32 c=0; c < qgramsNum; c++)
        {
            if(!name.mid(c,2).contains(" "))
            {
                nameQgrams[tempQgramsNum].append(name.mid(c,2).toLower().toLatin1());   // and filling with the two chars
                tempQgramsNum++;
            }
        }
        qgramsNum = tempQgramsNum;
        //store the indexes of all q-grams in an array
        for (quint32 num=0; num < qgramsNum; num++)
        {
            int c;
            qgIndex[num] = -1;
            for(c = 0; c < parent->q_gramsIdx_v.size(); c++)
            {
                if(parent->q_gramsIdx_v[c] == nameQgrams[num])
                {
                    qgIndex[num]= c;
                }
            }
        }
        delete[] nameQgrams;

        //fill the distances array (one each element of the db) with the number of identical q-grams
        for (quint64 row= 0; row < qgramsNum; row++)
        {
            if(qgIndex[row]!= -1)
            {
                for (int col=0; col < parent->invertedTable[qgIndex[row]].size(); col++)
                {
                    strDistances[parent->invertedTable[qgIndex[row]][col]]++;
                }
            }
        }
    }
    return;
}


//fn: resetCalcs
void CalcStringDistThread::resetCalcs()
{
    if(strDistances != NULL)
        delete[] strDistances;
    strDistances=NULL;
}

//fn: CalcSingleStringDist
//param quint64 Index  test between user and db[index]
//brief single distance calculation between the test string and N element of db
double CalcStringDistThread::CalcSingleStringDist(quint64 Index)
{ 
    double distance;

    if(strDistances[Index] > 0)
    {
        if(qgramsNum < db[Index].qgramsNum)
            distance = 100.0 - ((strDistances[Index]*100.0) /(double)qgramsNum) + ((double)(db[Index].qgramsNum - qgramsNum) );
        else
            distance = 100.0 - ((strDistances[Index]*100) /(double)db[Index].qgramsNum) + ((double)(qgramsNum - db[Index].qgramsNum));

        if(distance < 0)
            distance = 0;
    }
    else
        distance=100 + (double)abs((qgramsNum - db[Index].qgramsNum));

    return distance;
}


double CalcAgeDistThread::CalculateAgeDistance(element *peopleA, element *peopleB)
{
  return (double)abs(peopleA->age - peopleB->age);
}


//fn: CalculateLocalizationDistance
//param: two people element pointers
//brief: Distance between two gps point (stored as degree. decimal parts of degree)in km
double CalcGeoDistThread::CalculateLocalizationDistance(element *peopleA, element *peopleB)
{
    double _eQuatorialEarthRadius = 6378.1370L;
    double _d2r = (M_PI / 180.0L);
    double dlong = (peopleA->longitude - peopleB->longitude) * _d2r;
    double dlat = (peopleA->latitude - peopleB->latitude) * _d2r;
    double a = pow(sin(dlat/2.0), 2) + cos(peopleA->latitude * _d2r) * cos(peopleB->latitude * _d2r) * pow(sin(dlong/2.0), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    double distance = (double)(_eQuatorialEarthRadius * c);
    return distance;
}


void CalcEuclideanDistThread::SetParameter(results* resultsArray,
                                           quint64 firstElem,
                                           quint64 lastElem,
                                           double KAge,
                                           double KGeo,
                                           double KString)
{
    _res =  resultsArray;
    _first = firstElem;
    _last = lastElem;
    _kAge = KAge;
    _kGeo = KGeo;
    _kString = KString;
}


//fn: CalculateEuclideanDistance
//param results* res  results array pointer
void CalcEuclideanDistThread::CalculateEuclideanDistance(results* res)
{
    res->euclideanDistance = pow((res->stringDistance*_kString),2)+ pow((res->ageDistance*_kAge),2)+ pow((res->geoDistance*_kGeo),2);
    return;
}

/* fn: Q_GramsProcessing
 * brief:  this function parse all the element of the database, split it into qgrams and fill an inverted table to speed up
 *          the next searches
 */
void Q_GramsProcessingThread::Q_GramsProcessing()
{
    QByteArray *nameQgrams;
    for (; _first < _last; _first++)
    {
        //for every people name I update the QGrams inverted arrays
        _parent->people[_first].qgramsNum = (_parent->people[_first].name_t.length())-1;
        QString name = QString::fromUtf8(_parent->people[_first].name_t.c_str());
        nameQgrams = new QByteArray[_parent->people[_first].qgramsNum];      //two chr token array creation
        quint32 tempQgramsNum = 0;
        //split the string in qgrams, avoid qgrams that include spaces
        for (quint32 c=0; c < _parent->people[_first].qgramsNum; c++)
        {
            if(!name.mid(c,2).contains(" "))
            {
                nameQgrams[tempQgramsNum].append(name.mid(c,2).toLower().toLatin1());   // and filling with the two chars
                tempQgramsNum++;
            }
        }
        _parent->people[_first].qgramsNum = tempQgramsNum;

        //check if the QGrams are alredy present: otherwise I insert it into the vector
        bool checkPresence;
        for (int num=0; num < _parent->people[_first].qgramsNum; num++)
        {
            checkPresence = false;
            _mutex->lock();
            for(int c = 0; c < _parent->q_gramsIdx_v.size(); c++)
            {
                if(_parent->q_gramsIdx_v[c] == nameQgrams[num])
                {
                    _parent->invertedTable[c].push_back(_first);
                    checkPresence = true;
                }
            }
            if (checkPresence == false)
            {
                //insert qgram as vector into the inverted table
                if(nameQgrams[num].length() == 2)
                {
                    //new QGram found
                    QByteArray newQgram = nameQgrams[num];
                    //add new index to the vector
                    _parent->q_gramsIdx_v.push_back(newQgram);
                    QVector<quint64> tempVector;
                    //add a new column to the inverted table matrix
                    tempVector.push_back(_first);
                    _parent->invertedTable.push_back(tempVector);
                }
            }
            _mutex->unlock();
        }
        delete[] nameQgrams;
    }
}

void Q_GramsProcessingThread::SetParams(QMutex* mutex,
                                        MainWindow* parent,
                                        quint64 first,
                                        quint64 last)
{
    _mutex = mutex;
    _first = first;
    _last = last;
    _parent = parent;
}



/*************************************/
//Main window Methods
/*************************************/
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    people = NULL;
    root1 = NULL;
    //validator init
    QDoubleValidator* latValid = new QDoubleValidator(36.0L, 58.0L, 6, this);
    ui->leLat->setValidator(latValid);
    QDoubleValidator* lonValid = new QDoubleValidator(-10.0L, 40.0L, 6, this);
    ui->leLong->setValidator(lonValid);
    QDoubleValidator* kValid = new QDoubleValidator(this);
    ui->leKAge->setValidator(kValid);
    ui->leKString->setValidator(kValid);
    ui->leKGeo->setValidator(kValid);
    QIntValidator* ageValid = new QIntValidator(0,100,this);
    ui->leAge->setValidator(ageValid);
    QIntValidator* clusterValid = new QIntValidator(10,1000000,this);
    ui->leClusterThr->setValidator(clusterValid);
    //settings init and load
    m_sSettingsFile = "settings.ini";
    LoadSettings();
}


MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::closeEvent (QCloseEvent *event)
{
    if(people!= NULL)
        delete[] people;
    if(root1 != NULL)
        DeleteNode(root1);
    SaveSettings();
}


//fn: DeleteNode
//param: treeNode* start node to delete. all the children are deleted recursively starting from that point
//brief: delete all the children and the node itself
void MainWindow::DeleteNode(treeNode* currentNode)
{
    if (currentNode->minorBranch != NULL)
        DeleteNode(currentNode->minorBranch);
    if(currentNode->majorBranch != NULL)
        DeleteNode(currentNode->majorBranch);
    delete currentNode;
}


//fn: NodeResetSearch
//param: treeNode* start node to reset. all the children are called recursively starting from that point
//brief: Set the "searched" flag to false in all nodes of the tree
void MainWindow::NodeResetSearch(treeNode* currentNode)
{
    if (currentNode->minorBranch != NULL)
        NodeResetSearch(currentNode->minorBranch);
    if(currentNode->majorBranch != NULL)
        NodeResetSearch(currentNode->majorBranch);
    currentNode->searched=false;
}


//fn: LoadSettings
//Brief: load the K weight parameters and the cluster size from ini file
void MainWindow::LoadSettings()
{
    QSettings settings(m_sSettingsFile, QSettings::IniFormat);
    QString sText = settings.value("K_AGE", "").toString();
    if(sText.length() != 0)
        ui->leKAge->setText(sText);
    else
        ui->leKAge->setText(QString::number(K_AGE_WEIGHT));
    kAge=ui->leKAge->text().toDouble();
    sText = settings.value("K_GEO", "").toString();
    if(sText.length() != 0)
        ui->leKGeo->setText(sText);
    else
        ui->leKGeo->setText(QString::number(K_GEO_WEIGHT));
    kGeo = ui->leKGeo->text().toDouble();
    sText = settings.value("K_STR", "").toString();
    if(sText.length() != 0)
        ui->leKString->setText(sText);
    else
        ui->leKString->setText(QString::number(K_STRING_WEIGHT));
    kString=ui->leKString->text().toDouble();
    sText = settings.value("ClusterSize", "").toString();
    if(sText.length() != 0)
        ui->leClusterThr->setText(sText);
    else
        ui->leClusterThr->setText(QString::number(CLUSTER_MAX_ELEMENTS));
    maxClusterSize = ui->leClusterThr->text().toDouble();
}


//fn: SaveSettings
//Brief: save the K weight parameters and the cluster size from ini file
void MainWindow::SaveSettings()
{
    QSettings settings(m_sSettingsFile, QSettings::IniFormat);
    settings.setValue("K_AGE", ui->leKAge->text());
    settings.setValue("K_GEO", ui->leKGeo->text());
    settings.setValue("K_STR", ui->leKString->text());
    settings.setValue("ClusterSize", ui->leClusterThr->text());
}

//fn: Clusterize
//Brief: Start the clusterization procedure
void MainWindow::Clusterize()
{
    ui->teOutput->append("Start clusterization");
    taskTimer.start();
    totClassCreated = 0;
    //create the root of the kdtree
    root1 = new treeNode();
    root1->classNumber = 0;
    root1->numElements = peopleCount;
    root1->parent = NULL;
    root1->searched = false;
    ClusterizeMemBase(root1);
    ui->teOutput->append("Clusterization end in "+QString::number(taskTimer.elapsed())+ " ms");
    ui->teOutput->append(QString::number(totClassCreated)+ " classes created");
}


//fn: on_pbCreateDB_clicked
//Brief: Create a new in memory db
void MainWindow::on_pbCreateDB_clicked()
{
    CreateMemoryBase();
    //now is possible to do a simple search or clusterize the data
    ui->pbSearch->setEnabled(true);
    //now is possible to do a fast search..
    ui->pbFastSearch->setEnabled(true);
}


//fn: on_pbSearch_clicked
//Brief: Standard search procedure (linear or brute force)
void MainWindow::on_pbSearch_clicked()
{
    //take the paramentes of the new user that I have to compare
    //timer for calculate the search time
    taskTimer.start();
    QString tempName = ui->leName->text()+" "+ui->leSurname->text();
    ui->pbSearch->setEnabled(false);    //disable search buttons
    ui->pbFastSearch->setEnabled(false);
    userSearch.name_t = tempName.toStdString();
    userSearch.longitude = ui->leLong->text().toDouble();
    userSearch.latitude = ui->leLat->text().toDouble();
    userSearch.age=ui->leAge->text().toInt();
    //create the results array and store into every single result the correspondent People index
    distResults = new results[peopleCount];
    for (quint64 counter=0; counter < peopleCount; counter++)
        distResults[counter].element_id = counter;

    LinearSearch(peopleCount,distResults);

    ui->teOutput->append(QString::number(peopleCount)+"euclidean distance completed in "+QString::number(taskTimer.elapsed())+" ms");
    taskTimer.start();
    QuickSort(distResults,peopleCount, K_NEARER_NEIGH);
    ui->teOutput->append(QString::number(peopleCount)+" quick sorted completed in "+QString::number(taskTimer.elapsed())+" ms");

    //print the first N results
    for(int c = 0; c < K_NEARER_NEIGH; c++)
    {
        if(distResults[c].element_id < peopleCount)
        {
            tempName.clear();
            tempName = QString::fromUtf8(people[distResults[c].element_id].name_t.c_str());
            ui->teOutput->append("dis from "+QString::number(distResults[c].element_id)+" "+tempName+", age "
                                 +QString::number(people[distResults[c].element_id].age)
                                 +", lat"+QString::number(people[distResults[c].element_id].latitude)
                                 +", lon"+QString::number(people[distResults[c].element_id].longitude));
            ui->teOutput->append(QString::number(sqrt(distResults[c].euclideanDistance)));
        }
    }
    ui->pbSearch->setEnabled(true); //enable search button
    ui->pbFastSearch->setEnabled(true);
    strThread.resetCalcs(); //needed to reset previous search calcs and prepare for a new search

    if (distResults!= NULL)
        delete[] distResults;
    distResults = NULL;
}


void MainWindow::on_horizontalSlider_valueChanged(int value)
{
    ui->lblRecordsNum->setText(QString::number(value)+" million of records");
}

/*
 * fn: QuickSort
 * param: results* rs           pointer to the results array
 *        quint64 numResults    size of results's array (number of elements)
 *        quint64 threshold     threshold (number of elements to sort)
 * brief: partial sort of teh results array: the results from 0 to threshold are smaller than
 *          the results from threshold to numResults
 */
void MainWindow::QuickSort(results* rs, quint64 numResults, quint64 threshold)
{
    if (rs != NULL && numResults > threshold)
    {
        quint64 tempThs;
        float thsValue;
        do{
            tempThs = threshold;
            thsValue = rs[tempThs].euclideanDistance;
            for (quint64 count=0; count < numResults; count++)
            {
                if (count < tempThs)
                {
                    //under threshold
                    if (rs[count].euclideanDistance > thsValue)
                    {
                        std::swap(rs[tempThs],rs[count]);
                        //update the threshold
                        tempThs = count;
                    }
                }
                if (count > tempThs)
                {
                    //under threshold
                    if (rs[count].euclideanDistance < thsValue)
                    {
                        std::swap(rs[tempThs],rs[count]);
                        //update the threshold
                        tempThs = count;
                    }
                }
            }
        }while (threshold != tempThs);
    }
}


/* fn: CreateMemoryDatabase
 * brief:  this function create the inmemory database and fill it with randomic elements
 */
void MainWindow::CreateMemoryBase()
{

  //create the in memory database
  QStringList  names;
  QStringList surnames;
  taskTimer.start();
  //load the names/surnames database

  //load the population in memory database
  peopleCount = ui->horizontalSlider->value()* ELEMENT_NUMBER;
  people = new element[peopleCount];

  ui->teOutput->append("Start in memory database creation of "+QString::number(peopleCount)+" elements");
  //load CSV file of names
  QString fileName = "C:\\Personal\\exercise\\test_v2\\CSV_Database_of_First_Names.txt";
  QFile file(fileName);
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    while (!file.atEnd())
    {
        names << file.readLine();
        QCoreApplication::processEvents();
    }
    file.close();
  }
  //remove all "\n" at the end of any name (if present)
  for (int nCount=0; nCount < names.size();nCount++)
      if(names[nCount].endsWith("\n"))
          names[nCount].chop(1);

  //load CSV file of last names
  fileName = "C:\\Personal\\exercise\\test_v2\\CSV_Database_of_Last_Names.txt";
  file.setFileName(fileName);
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    while (!file.atEnd())
    {
      surnames << file.readLine();
      QCoreApplication::processEvents();
    }
    file.close();
  }
  //remove all "\n" at the end of any name (if present)
  for (int nCount=0; nCount < surnames.size();nCount++)
      if(surnames[nCount].endsWith("\n"))
          surnames[nCount].chop(1);

  for (quint64 count= 0; count < peopleCount; count++)
  {
      people[count].age = (quint8)(rand() % 100);
      QString tempName = names[rand() % names.size()] +" "+ surnames[rand() % surnames.size()];
      people[count].name_t = tempName.toStdString();
      people[count].latitude = ((((float) rand()) / (float) RAND_MAX)*22.0) +36.0;
      people[count].longitude = ((((float) rand()) / (float) RAND_MAX)*50.0) -10.0;
      people[count].attributeClass = 0;
      QCoreApplication::processEvents();
  }
  ui->teOutput->append("database creation completed in "+QString::number(taskTimer.elapsed())+" ms");

  Clusterize();
  //tree order
  StartTreeSort();

  QGramsProcessing();
}


/* fn: QGramsProcessing
 * brief:  this function start the thread taht calcs all the q-grams occurrency filling an inverted table
 *         Best performance using a single thread
 */
void MainWindow::QGramsProcessing()
{
    QMutex vectorMutex;
    Q_GramsProcessingThread* fastQ_GramsProc = new Q_GramsProcessingThread;
    taskTimer.start();
    ui->teOutput->append("Start Q-GRams processing");
    fastQ_GramsProc->SetParams(&vectorMutex,
                                 this,
                                 0,
                                 peopleCount);
    fastQ_GramsProc->start(QThread::HighPriority);

    bool finished=false;
    while(finished == false)     //wait untile the thread run is finished
    {
        finished=true;
        if (fastQ_GramsProc->isRunning())
            finished=false;
        QCoreApplication::processEvents();
    }
    delete fastQ_GramsProc;
    ui->teOutput->append("Q-GRams Processing completed in "+QString::number(taskTimer.elapsed())+" ms");
}


void MainWindow::StartTreeSort()
{
    ui->teOutput->append("Start tree sorting");
    taskTimer.start();
    quint64 startPoint = 0;
    element* peopleCopy = new element[peopleCount];
    //create the root of the kdtree
    if(root1 != NULL)
        TreeSort(root1,&startPoint,peopleCopy);
    element* swapper = people;
    people = peopleCopy;
    delete[] swapper;
    ui->teOutput->append("Tree sorting end in "+QString::number(taskTimer.elapsed())+ " ms");
    ui->teOutput->append(QString::number(totClassCreated)+ " classes created");
}

/* fn: TreeSort
 * param:  treeNode* branch     tree node to process
 *         quint64* startPoint  position of the first free element of peopleCopy array
 *         element* peopleCopy  pointer to the ordinated element array
 * brief:  this function reorder all the elements of the people array filling a separate copy of that array
 *         the tree node leaf value are updated with the index of the first element of that class
 */
void MainWindow::TreeSort(treeNode* branch, quint64* startPoint, element* peopleCopy)
{
    if(branch!= NULL)
    {
        //only leaf node has elements, so I need to reach the leaf
        if(branch->minorBranch == NULL && branch->majorBranch == NULL)
        {
            // index is the first element
            branch->indexPos = *startPoint;
            //swap the indexes
            for (quint64 counter=0; counter < peopleCount; counter++)
            {
                if (people[counter].attributeClass == branch->classNumber)
                {
                    peopleCopy[*startPoint] = people[counter];
                    //memcpy(&peopleCopy[*startPoint],&people[counter], sizeof(element));
                    (*startPoint)++;
                }
            }
        }
        else
        {
            if(branch->minorBranch != NULL)
                TreeSort(branch->minorBranch, startPoint,peopleCopy);
            if(branch->majorBranch != NULL)
                TreeSort(branch->majorBranch, startPoint,peopleCopy);
        }
    }
}


/*  ClusterizeMemBase(treeNode* branch)
 * recursive function that divide the memory in a KD-tree clusterized memory
 * the parameters used for the classification are latitude- longitude and age
 * name string is not used, because it's not an euclidean measure
 */
void MainWindow::ClusterizeMemBase(treeNode* branch)
{
    quint64 minorElem = 0;
    quint64 majorElem = 0;
    quint64 classElement = 0;

    double midLat = 0;
    double maxLat = -100;
    double minLat = 100;
    double midLon = 0;
    double maxLon = -100;
    double minLon = 100;
    double midAge = 0;
    double maxAge = 0;
    double minAge = 100;

    //calc the min /max /mid values for the subclass
    for(quint64 counter = 0; counter < peopleCount; counter ++)
    {
        if(people[counter].attributeClass == branch->classNumber)
        {
            classElement++;
            midLat += people[counter].latitude;
            midLon += people[counter].longitude;
            midAge += people[counter].age;
            if(people[counter].latitude > maxAge)
                maxLat = people[counter].latitude;
            if(people[counter].latitude < minLat)
                minLat = people[counter].latitude;
            if(people[counter].longitude > maxLon)
                maxLon = people[counter].longitude;
            if(people[counter].longitude < minLon)
                minLon = people[counter].longitude;
            if(people[counter].age > maxAge)
                maxAge = people[counter].age;
            if(people[counter].age < minAge)
                minAge = people[counter].age;
        }
    }

    midLat = midLat / classElement;
    midLon = midLon / classElement;
    midAge = midAge / classElement;
    //select the parameter with the highest variance.
    int param = 0;
    if(((maxAge-minAge)*0.1) > ((maxLat-minLat)) && ((maxAge-minAge)*0.1) > ((maxLon-minLon)))
        param = 2; //age selected
    if(((maxLon-minLon)) > ((maxAge-minAge)*0.1) && ((maxLon-minLon)) > ((maxLat-minLat)))
        param = 1;
    //otherwise Lat selected param=0

    QCoreApplication::processEvents();

    switch(param)
    {
    case 0: //Latitude
    {
        double closerLat = 100;
        //search the element closest to the mid value
        for(quint64 counter = 0; counter < peopleCount; counter ++)
        {
            if(people[counter].attributeClass == branch->classNumber)
            {
                if(closerLat > fabs(people[counter].latitude - midLat))
                {
                    branch->indexPos = counter;
                    closerLat = fabs(people[counter].latitude - midLat);
                }
            }
        }
        //split the classes in two parts
        for(quint64 counter = 0; counter < peopleCount; counter ++)
        {
            if(people[counter].attributeClass == branch->classNumber)
            {
                if(people[counter].latitude < people[branch->indexPos].latitude)
                {
                    people[counter].attributeClass = totClassCreated+1;
                    minorElem++;
                }
                else
                {
                   people[counter].attributeClass = totClassCreated+2;
                   majorElem++;
                }
            }
        }
        break;
    }
    case 1: //Longitude
    {
        double closerLon = 100;
        //search the element closest to the mid value
        for(quint64 counter = 0; counter < peopleCount; counter ++)
        {
            if(people[counter].attributeClass == branch->classNumber)
            {
                if(closerLon > fabs(people[counter].longitude - midLon))
                {
                    branch->indexPos = counter;
                    closerLon = fabs(people[counter].longitude - midLon);
                }
            }
        }
        //split the classes in two parts
        for(quint64 counter = 0; counter < peopleCount; counter ++)
        {
            if(people[counter].attributeClass == branch->classNumber)
            {
                if(people[counter].longitude < people[branch->indexPos].longitude)
                {
                    people[counter].attributeClass = totClassCreated+1;
                    minorElem++;
                }
                else
                {
                   people[counter].attributeClass = totClassCreated+2;
                   majorElem++;
                }
            }
        }
        break;
    }
    case 2: //Age
    {
        double closerAge = 100;
        //search the element closest to the mid value
        for(quint64 counter = 0; counter < peopleCount; counter ++)
        {
            if(people[counter].attributeClass == branch->classNumber)
            {
                if(closerAge > fabs(people[counter].age - midAge))
                {
                    branch->indexPos = counter;
                    closerAge = fabs(people[counter].age - midAge);
                }
            }
        }
        //split the classes in two parts
        for(quint64 counter = 0; counter < peopleCount; counter ++)
        {
            if(people[counter].attributeClass == branch->classNumber)
            {
                if(people[counter].age < people[branch->indexPos].age)
                {
                    people[counter].attributeClass = totClassCreated+1;
                    minorElem++;
                }
                else
                {
                   people[counter].attributeClass = totClassCreated+2;
                   majorElem++;
                }
            }
        }
        break;
    }
    default:
        //error
        break;
    }

    branch->param = param;
    branch->minorBranch = new treeNode();
    branch->majorBranch = new treeNode();
    branch->minorBranch->numElements = minorElem;
    branch->majorBranch->numElements = majorElem;
    branch->minorBranch->parent = branch;
    branch->majorBranch->parent = branch;
    branch->minorBranch->classNumber = totClassCreated+1;
    branch->majorBranch->classNumber = totClassCreated+2;
    totClassCreated +=2;

    if(majorElem > maxClusterSize)
        ClusterizeMemBase(branch->majorBranch);

    if(minorElem > maxClusterSize)
        ClusterizeMemBase(branch->minorBranch);
}


void MainWindow::on_pbFastSearch_clicked()
{
    //fast search algorithm
    //timer for calculate the search time
    taskTimer.start();
    //take the paramentes of the new user that I have to compare
    QString tempName = ui->leName->text()+" "+ui->leSurname->text();
    ui->pbFastSearch->setEnabled(false);
    userSearch.name_t = tempName.toStdString();
    userSearch.longitude = ui->leLong->text().toDouble();
    userSearch.latitude = ui->leLat->text().toDouble();
    userSearch.age=ui->leAge->text().toInt();

    strThread.db = people;
    strThread.test = userSearch;
    strThread.first = 0;
    strThread.last = peopleCount;
    strThread.parent = this;
    //precalculate all the string distances
    //strThread.CalcStringDist();

    treeNode* currNode = root1;
    //straight tree check /until I find the leaf
    StraightTreeCheck(currNode);
    //then quick search inside that leaf and quick sort of the results. Now the 10th result is the result to beat
    //distResults = new results[bestMatch->numElements];
    distResults = new results[peopleCount];

    quint64 resCount=0;
    for (quint64 count = bestMatch->indexPos ; count < (bestMatch->indexPos+bestMatch->numElements); count++)
    {
        //if (people[count].attributeClass == bestMatch->classNumber)
        //{
            distResults[resCount].element_id = count;
            resCount++;
        //}
    }
    //first Best Result found. Check the 10th distance in this leaf
    LinearSearch(bestMatch->numElements,distResults);
    QuickSort(distResults,bestMatch->numElements,K_NEARER_NEIGH);
    //now the distance[10].euclidean distance is the distance to beat. so I check
    //all the other element in reverse searching for a closer value
    bestMatch->searched = true;
    currNode = bestMatch->parent;
    quint64 position=K_NEARER_NEIGH+1;
    ReverseTreeSearch(distResults, bestMatch->numElements, currNode, &position);
    QuickSort(distResults,bestMatch->numElements, K_NEARER_NEIGH);
    QCoreApplication::processEvents();

    ui->teOutput->append(QString::number(peopleCount)+"euclidean distance completed in "+QString::number(taskTimer.elapsed())+" ms");

    for(int c = 0; c < K_NEARER_NEIGH; c++)
    {
        if(distResults[c].element_id < peopleCount)
        {
            QString tempName = QString::fromUtf8(people[distResults[c].element_id].name_t.c_str());
            ui->teOutput->append("dis from "+QString::number(distResults[c].element_id)+" "+tempName+", age "
                                 +QString::number(people[distResults[c].element_id].age)
                                 +", lat"+QString::number(people[distResults[c].element_id].latitude)
                                 +", lon"+QString::number(people[distResults[c].element_id].longitude));
            ui->teOutput->append(QString::number(sqrt(distResults[c].euclideanDistance)));
        }
    }
    QCoreApplication::processEvents();
    ui->pbFastSearch->setEnabled(true);

    NodeResetSearch(root1);
    strThread.resetCalcs();
    if(distResults!= NULL)
        delete[] distResults;
    distResults= NULL;
}

/* StraightTreeCheck
 * param: treeNode starting point
 * brief: recursive function that search inside the KD-tree the index approximately nearer to the user (best result)
 *        in that case is calulated the gps distance and the age distance only
 */
void MainWindow::StraightTreeCheck(treeNode* currNode)
{
    if(currNode!=NULL)
    {
        if(currNode->minorBranch != NULL && currNode->majorBranch != NULL)
        {
            switch(currNode->param)
            {
            case 0: // latitude
                if(userSearch.latitude > people[currNode->indexPos].latitude)
                    StraightTreeCheck(currNode->majorBranch);
                else
                    StraightTreeCheck(currNode->minorBranch);
                break;
            case 1: // longitude
                if(userSearch.longitude > people[currNode->indexPos].longitude)
                    StraightTreeCheck(currNode->majorBranch);
                else
                    StraightTreeCheck(currNode->minorBranch);
                break;
                break;
            case 2: // age
                if(userSearch.age > people[currNode->indexPos].age)
                    StraightTreeCheck(currNode->majorBranch);
                else
                    StraightTreeCheck(currNode->minorBranch);
                break;
                break;
            }
        }
        else
            bestMatch = currNode;
    }
}


/* LinearSearch
 * param: quint64 numElements  number of elements of the rasults array
 *        results* distance  results array pointer
 * brief: calculate the euclidean distance for all the elements of the distance array
 *        to the userSearch element (test element)
 */
void MainWindow::LinearSearch(quint64 numElements, results* distance)
{
    //thread Creation
    //every thread calculate 1/4 (4 is the number of cores) of the database
    CalcAgeDistThread* cAgeTh = new CalcAgeDistThread[NUMBER_OF_CORE];
    CalcGeoDistThread* cGeoTh = new CalcGeoDistThread[NUMBER_OF_CORE];
    CalcStringDistThread* cStrTh = &strThread;//new CalcStringDistThread;
    CalcEuclideanDistThread* cEucTh = new CalcEuclideanDistThread[NUMBER_OF_CORE];

    //thread init
    for (quint8 count=0; count < NUMBER_OF_CORE; count++)
    {
        cAgeTh[count].db = people;
        cAgeTh[count].test = userSearch;
        cAgeTh[count].res = distance;
        cAgeTh[count].first = count*(numElements/NUMBER_OF_CORE);
        cAgeTh[count].last = (count+1)*(numElements/NUMBER_OF_CORE);

        cGeoTh[count].db = people;
        cGeoTh[count].test = userSearch;
        cGeoTh[count].res = distance;
        cGeoTh[count].first = count*(numElements/NUMBER_OF_CORE);
        cGeoTh[count].last = (count+1)*(numElements/NUMBER_OF_CORE);

        cEucTh[count].SetParameter(distance,count*(numElements/NUMBER_OF_CORE),(count+1)*(numElements/NUMBER_OF_CORE),kAge,kGeo,kString);
    }

    //the string thread cannot be easily splitted, so a single thread is used
    cStrTh->db = people;
    cStrTh->test = userSearch;
    cStrTh->res = distance;
    cStrTh->first = 0;
    cStrTh->last = numElements;
    cStrTh->parent = this;

    // last element of the N-th thread is set to the maximum numElement value
    cAgeTh[NUMBER_OF_CORE-1].last = numElements;
    cGeoTh[NUMBER_OF_CORE-1].last = numElements;
    cEucTh[NUMBER_OF_CORE-1].SetLastElement(numElements);

    //check the connection between threads and Main thread

    for (quint8 count=0; count < NUMBER_OF_CORE; count++)
    {
        cAgeTh[count].start(QThread::HighPriority);
        cGeoTh[count].start(QThread::HighPriority);
    }
    cStrTh->start(QThread::HighPriority);

    bool finished=false;
    //test del cactus..
    while(finished == false)
    {
        finished=true;
        for (quint8 count=0; count < NUMBER_OF_CORE; count++)
        {
            if (cAgeTh[count].isRunning())
                finished=false;
            if (cGeoTh[count].isRunning())
                finished=false;
        }
        if (cStrTh->isRunning())
            finished=false;
        QCoreApplication::processEvents();
    }

    for (quint8 count=0; count < NUMBER_OF_CORE; count++)
        cEucTh[count].start(QThread::HighPriority);

    finished=false;
    while(finished == false)
    {
        finished=true;
        for (quint8 count=0; count < NUMBER_OF_CORE; count++)
            if (cEucTh[count].isRunning())
                finished=false;
        QCoreApplication::processEvents();
    }

    delete[] cAgeTh;
    delete[] cGeoTh;
    delete[] cEucTh;
}


/* ReverseTreeSearch
 * param:  results* distances  array where to store the results
 *         quint32 reElements  number of elements in result array
 *         treeNode* currNode  pointer to the node to search in
 *         quint64* position   index of the next value that is possible to modify into the results array
 * brief:  This recursive function search among the element of every single leaf node, starting from the closer
 *         to the farer. the first node given to the function couldn't be a leaf, otherwise the tree is completely not
 *         searched
 */
void MainWindow::ReverseTreeSearch(results* distances, quint32 resElements, treeNode* currNode, quint64* position)
{
    //select the node
    //is that a leaf?
    if(currNode->minorBranch == NULL && currNode->majorBranch == NULL)
    {
        if(currNode->searched == false) //not yet searched in
        {
            quint32 resIndex = *position;
            double threshold = (distances[K_NEARER_NEIGH].euclideanDistance*NEARER_PERCENTAGE);
            distances[resIndex].ageDistance = 0;
            distances[resIndex].stringDistance = 0;
            distances[resIndex].geoDistance = 0;
            distances[resIndex].euclideanDistance = 0;
            distances[resIndex].element_id = 0;


            quint64 endOfNode = currNode->indexPos+currNode->numElements;
            for (quint64 counter=currNode->indexPos; counter < endOfNode; counter ++)
            {
                if(people[counter].attributeClass != bestMatch->classNumber)
                {
                    distances[resIndex].ageDistance = ageThread.CalculateAgeDistance(&userSearch,&people[counter]);
                    distances[resIndex].euclideanDistance = (distances[resIndex].ageDistance*kAge)*(distances[resIndex].ageDistance*kAge) ;
                    if(distances[resIndex].euclideanDistance < threshold)
                    {
                        //calculate the string distance
                        distances[resIndex].stringDistance = strThread.CalcSingleStringDist(counter);
                        distances[resIndex].euclideanDistance+=(distances[resIndex].stringDistance*kString)*(distances[resIndex].stringDistance*kString);
                        if(distances[resIndex].euclideanDistance < threshold)
                        {
                            //calculate the geo distance and the euclidean too
                            distances[resIndex].geoDistance = geoThread.CalculateLocalizationDistance(&userSearch,&people[counter]);
                            distances[resIndex].euclideanDistance+= (distances[resIndex].geoDistance*kGeo)*(distances[resIndex].geoDistance*kGeo);
                            distances[resIndex].element_id = counter;
                            //if maximum not reached, go next
                            if(resIndex < resElements)
                                resIndex++;
                            else
                            {
                               //result array completely filled, I sort again and restart from the 10th element
                               QuickSort(distances,resElements, K_NEARER_NEIGH);
                               resIndex= K_NEARER_NEIGH+1;
                               threshold = distances[K_NEARER_NEIGH].euclideanDistance*NEARER_PERCENTAGE;
                            }
                            distances[resIndex].ageDistance = 0;
                            distances[resIndex].stringDistance = 0;
                            distances[resIndex].geoDistance = 0;
                            distances[resIndex].euclideanDistance = 0;
                        }
                    }
                }
            }
            *position = resIndex;
            currNode->searched = true;

        }
    }
    else
    {
         currNode->searched = true;
        //not a leaf.. do I have branch not searched yet?
        if(currNode->minorBranch->searched == false)
           ReverseTreeSearch(distances, resElements, currNode->minorBranch, position);
        if(currNode->majorBranch->searched == false)
           ReverseTreeSearch(distances, resElements, currNode->majorBranch, position);

        //now select the parent if present
        if(currNode->parent!= NULL)
        {
            if(currNode->parent->searched == false)
                ReverseTreeSearch(distances, resElements, currNode->parent, position);
        }
    }
}


void MainWindow::QuickSearch(results* distances, quint32 resElements)
{
    quint32 resIndex = K_NEARER_NEIGH+1;
    double threshold = (distances[K_NEARER_NEIGH].euclideanDistance*NEARER_PERCENTAGE);
    distances[resIndex].ageDistance = 0;
    distances[resIndex].stringDistance = 0;
    distances[resIndex].geoDistance = 0;
    distances[resIndex].euclideanDistance = 0;
    distances[resIndex].element_id = 0;

    for (quint64 counter=0; counter < peopleCount; counter ++)
    {
        if(people[counter].attributeClass != bestMatch->classNumber)
        {
            distances[resIndex].ageDistance = ageThread.CalculateAgeDistance(&userSearch,&people[counter]);
            distances[resIndex].euclideanDistance = (distances[resIndex].ageDistance*kAge)*(distances[resIndex].ageDistance*kAge) ;
            if(distances[resIndex].euclideanDistance < threshold)
            {
                //calculate the string distance
                distances[resIndex].stringDistance = strThread.CalcSingleStringDist(counter);
                distances[resIndex].euclideanDistance+=(distances[resIndex].stringDistance*kString)*(distances[resIndex].stringDistance*kString);
                if(distances[resIndex].euclideanDistance < threshold)
                {
                    //calculate the geo distance and the euclidean too
                    distances[resIndex].geoDistance = geoThread.CalculateLocalizationDistance(&userSearch,&people[counter]);
                    distances[resIndex].euclideanDistance+= (distances[resIndex].geoDistance*kGeo)*(distances[resIndex].geoDistance*kGeo);
                    distances[resIndex].element_id = counter;
                    //if maximum not reached, go next
                    if(resIndex < resElements)
                        resIndex++;
                    else
                    {
                       //result array completely filled, I sort again and restart from the 10th element
                       QuickSort(distances,resElements, K_NEARER_NEIGH);
                       resIndex= K_NEARER_NEIGH+1;
                       threshold = distances[K_NEARER_NEIGH].euclideanDistance*NEARER_PERCENTAGE;
                    }
                    distances[resIndex].ageDistance = 0;
                    distances[resIndex].stringDistance = 0;
                    distances[resIndex].geoDistance = 0;
                    distances[resIndex].euclideanDistance = 0;
                }
            }
        }
    }
    QuickSort(distances,resElements, K_NEARER_NEIGH);
    QCoreApplication::processEvents();
}




// Update value functions
void MainWindow::on_pbClear_clicked()
{
    ui->teOutput->clear();
}

void MainWindow::on_leKString_editingFinished()
{
    kString = ui->leKString->text().toDouble();
}

void MainWindow::on_leKGeo_editingFinished()
{
    kGeo = ui->leKGeo->text().toDouble();
}

void MainWindow::on_leKAge_editingFinished()
{
    kAge = ui->leKAge->text().toDouble();
}

void MainWindow::on_leClusterThr_editingFinished()
{
    maxClusterSize = ui->leClusterThr->text().toDouble();
}
