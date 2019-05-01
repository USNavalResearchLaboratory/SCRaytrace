
#include "config.h"

//#if defined (HAVE_BOOST_THREAD) && defined (HAVE_BOOST)

#if defined Boost_FOUND
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>
#endif


#include "scene.h"
#include "rtmiscfunc.h"
#include "constant.h"

#include <iostream>


Scene::Scene() 
{
  abs = Cbasis(Cvec(0,0,0),0,0,0);
  modelposition = ModelPosition();
  neonly = false;
  quiet = false;
  pphy = physicsSelect(THOMSON); // -- set the physics to Thomson scattering by default
  pphy->setParentScene(this);
  pmod = 0;
  runfracmax = false;
  runDumpInteg = false;
  
// #if !defined (HAVE_BOOST_THREAD) || !defined (HAVE_BOOST)
 
#if !defined Boost_FOUND
  std::cout << "WARNING: library not compiled with threading: normal raytracing will be performed." << std::endl;
#endif

}


void Scene::setPhysics(PhysicsType phytype)
{
  delete pphy;
  pphy=physicsSelect(phytype);
  pphy->setParentScene(this);
}

void Scene::setFracMax(float fracmax)
{
      if (fracmax <= 0.) {
        this->fracmax = 0.;
        this->runfracmax = false;
      } else if (fracmax > 1.) {
        this->fracmax = 1.;
        this->runfracmax = true;
      } else {
        this->fracmax = fracmax;
        this->runfracmax = true;
      }
}


void Scene::setDumpIntegrand(bool flagrunDumpInteg, float *pintegrand)
{
  if (flagrunDumpInteg) {
    this->runDumpInteg = true;
    this->pintegrand = pintegrand;
  } else {
    this->runDumpInteg = false;
  }
  
  printvar(this->runDumpInteg);
  
}


// define methods if compilation with boost
//#if defined (HAVE_BOOST_THREAD) && defined (HAVE_BOOST)
#if defined Boost_FOUND

boost::mutex io_mutex;

void Scene::losintegray(const unsigned int &i,const unsigned int &j,const unsigned int &threadid)
{
    pisrunning[threadid]=1;

    losinteg(i,j);
    if (runfracmax) losintegtofrac(i,j);

    pisrunning[threadid]=0;

}



void Scene::losintegchunk(const unsigned int &chunkid,const unsigned int &threadid)
{
	{

    boost::this_thread::disable_interruption di;

    pisrunning[threadid]=1;

    unsigned int adjustedchunksize=chunksize;
    if ((chunkid+1)==nbchunk) adjustedchunksize+=lastchunkremain;


   // -- progression
    float progresspercent=0.2;
    int progressflag=(int) ((float) adjustedchunksize * progresspercent);
    float progresspass=progresspercent;



    for (unsigned int k=0;k<adjustedchunksize;k++)
    {

        if ((quiet != 1) && (k > progressflag))
        {
            {
                boost::mutex::scoped_lock lock(io_mutex);
                cout << "Chunk "<< chunkid << " : " << progresspass*100 << "% " << endl;
            }
            progresspass+=progresspercent;
            progressflag=(int) ((float) adjustedchunksize * progresspass);
        }

        unsigned int pos=chunkid*chunksize+k;
        unsigned int i=pos % camera.detector.sxpix;
        unsigned int j=pos / camera.detector.sxpix;

        if (runDumpInteg) losintegDumpIntegrand(i,j); else losinteg(i,j);
        if (runfracmax) losintegtofrac(i,j);

    }

    cout << "Chunk "<< chunkid << " : 100% " << endl;
	
    pisrunning[threadid]=0;

  }
}




void Scene::computeImagebyRay(float *btot,float *bpol,float *netot,const unsigned int nbthread)
{

    this->btot=btot;
    this->bpol=bpol;
    this->netot=netot;

    // -- progression
    float progresspercent=0.2;
    int progressflag=(int) ((float) camera.detector.sypix * progresspercent);
    float progresspass=progresspercent;

    boost::thread *pthread=new boost::thread [nbthread];
    pisrunning=new bool [nbthread];

    // ---- init isrunning
    for (int i=0;i<nbthread;i++) pisrunning[i]=0;

    if (quiet != 1) std::cout << "NB Core : " << pthread[0].hardware_concurrency() << std::endl;

    unsigned int threadid=0;

    for (unsigned int j=0;j<camera.detector.sypix;j++)
    {

        // -- print progression
        if ((quiet != 1) && (j > progressflag))
        {
            cout << progresspass*100 << "% " << endl;
            progresspass+=progresspercent;
            progressflag=(int) ((float) camera.detector.sypix * progresspass);
        }

        for (unsigned int i=0;i<camera.detector.sxpix;i++)
        {

            // ---- wait
            while (pisrunning[threadid])
            {
                threadid++;
                if (threadid >= nbthread) threadid=0;
            }

            pthread[threadid].join();
            pisrunning[threadid]=1;
            boost::thread t(boost::bind(&Scene::losintegray,this,i,j,threadid));
            pthread[threadid]=boost::move(t);

        }
    }


// ---- wait for completion
    for (int i=0;i<nbthread;i++) pthread[i].join();

    delete[] pisrunning;
    delete[] pthread;

}


void Scene::computeImagebyChunk(float *btot,float *bpol,float *netot,const unsigned int nbthread,const unsigned int nbchunk)
{

    this->btot=btot;
    this->bpol=bpol;
    this->netot=netot;

    cout << "Obs dist to sun [Rsun] "<< obs.o.mag() << endl;
    cout << "1 AU [Rsun] "<< ONEAU_RSUN << endl;

    boost::thread *pthread=new boost::thread[nbthread];
    pisrunning=new bool[nbthread];

    // ---- init isrunning
    for (int i=0;i<nbthread;i++) pisrunning[i]=0;

    if (quiet != 1) std::cout << "NB Core : " << pthread[0].hardware_concurrency() << std::endl;

    unsigned int threadid=0;
    nbpix=camera.detector.sxpix*camera.detector.sypix;
    chunksize=nbpix / nbchunk;
    lastchunkremain=nbpix % nbchunk;
    this->nbchunk=nbchunk;

    for (unsigned int i=0;i<nbchunk;i++)
    {

        // ---- wait
        while (pisrunning[threadid])
        {
            threadid++;
            if (threadid >= nbthread) threadid=0;
        }

        pthread[threadid].join();
        pisrunning[threadid]=1;
        boost::thread t(boost::bind(&Scene::losintegchunk,this,i,threadid));
        pthread[threadid]=boost::move(t);


    }

    // ---- wait for completion
    for (int i=0;i<nbthread;i++) pthread[i].join();

    delete[] pisrunning;
    delete[] pthread;


}



// ---- provide alternative code if no boost
#else 



void Scene::losintegray(const unsigned int &i,const unsigned int &j,const unsigned int &threadid)
{

}



void Scene::losintegchunk(const unsigned int &chunkid,const unsigned int &threadid)
{

}




void Scene::computeImagebyRay(float *btot,float *bpol,float *netot,const unsigned int nbthread)
{

    this->btot=btot;
    this->bpol=bpol;
    this->netot=netot;

    // -- progression
    float progresspercent=0.2;
    int progressflag=(int) ((float) camera.detector.sypix * progresspercent);
    float progresspass=progresspercent;

    for (unsigned int j=0;j<camera.detector.sypix;j++)
    {

        // -- print progression
        if (((quiet) <= 1) && (j > progressflag))
        {
            cout << progresspass*100 << "% " << endl;
            progresspass+=progresspercent;
            progressflag=(int) ((float) camera.detector.sypix * progresspass);
        }

        for (unsigned int i=0;i<camera.detector.sxpix;i++) losinteg(i,j);
        if (runfracmax) for (unsigned int i=0;i<camera.detector.sxpix;i++) losintegtofrac(i,j);
    }

}


void Scene::computeImagebyChunk(float *btot,float *bpol,float *netot,const unsigned int nbthread,const unsigned int nbchunk)
{

computeImagebyRay(btot,bpol,netot,nbthread);

}

#endif

