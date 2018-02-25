// $Id: main.cpp,v 1.12 2010-09-17 15:24:43 thernis Exp $

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "constant.h"
#include "Cvec.h"
#include "Cmat.h"
#include "Cbasis.h"
#include "Cdvoxel.h"
#include "Clos.h"
#include "CModelBase.h"
#include "CControlPoint.h"
#include "Cnurbs.h"
#include "raytrace.h"
#include "rtmiscfunc.h"
#include "rtcloud.h"
#include <algorithm>

using namespace std;

#define NBCALL 500000000

void callbyvalue(Cvec v)
{
	float a=v.mag();
}

void callbyreference(Cvec &v)
{
	float a=v.mag();
}

void callbypointer(Cvec *pv)
{
	float a=pv->mag();
}

void callbyconstreference(const Cvec &v)
{
	float a=v.mag();
}




//! Main program, mostly for testing
int main(int argc,char **argv) {

    // ---- select the density model requested by the user
    int modelid=1;

    if (argc == 1) {
        cout << "Using default modelid: 1" << endl;
    } else {
        argv++;
        modelid=atoi(*argv);
    }

		Cvec vcall(10,20,30);

		if (modelid == -3) {
				cout << "Call by value...:" << NBCALL <<endl;
				for(unsigned int i=0;i<NBCALL;i++) callbyvalue(2*vcall);
				return 0;
		}

		if (modelid == -4) {
				cout << "Call by reference...:" << NBCALL <<endl;
				for(unsigned int i=0;i<NBCALL;i++) callbyreference(vcall);
				return 0;
		}

		if (modelid == -5) {
				cout << "Call by pointer...:" << NBCALL <<endl;
				for(unsigned int i=0;i<NBCALL;i++) callbypointer(&vcall);
				return 0;
		}

		if (modelid == -6) {
				cout << "Call by const reference...:" << NBCALL <<endl;
				for(unsigned int i=0;i<NBCALL;i++) callbyconstreference(2*vcall);
				return 0;
		}

    // -- development test is modelid == -1
    if (modelid == -1) {

        float xxx=0,yyy=1,zzz=1.,rrr,beta,ttt;

        if (argc == 5) {
            argv++;
            xxx=atof(*argv);
            argv++;
            yyy=atof(*argv);
            argv++;
            zzz=atof(*argv);
        }

        cvcoord(xxx,yyy,0,&rrr,&beta,&ttt);

        Cvec v(xxx,yyy,zzz);

        cout << "xxx : " << xxx << endl;
        cout << "yyy : " << yyy << endl;
        cout << "rrr : " << rrr << endl;
        cout << "beta : " << beta*RADEG << endl;
        cout << "ttt : " << ttt << endl;
        /*
                float rs=9.,rc=7.;
                cout << "rs : " << rs << endl;
                cout << "rc : " << rc << endl;
         
                Cvec vsk(rc*cos(beta),rs+rc*sin(beta),0);       
                Cvec vpc=v-vsk;
         
                cout << "v   : " << v << endl;
                cout << "Mag v   : " << v.norm() << endl;
          
                cout << "vsk : " << vsk << endl;
                cout << "Mag vsk : " << vsk.norm() << endl;
                cout << "vpc : " << vpc << endl;
                cout << "Mag vpc : " << vpc.norm() << endl;
                
                cout << "vpc[0] : " << vpc[0] << endl;
                cout << "vpc[1] : " << vpc[1] << endl;
                cout << "vpc[2] : " << vpc[2] << endl;
                
                
                Cmat mat(PI/4,1);
                cout << "mat : " << mat << endl;
                
                cout << "mat[0][0] : " << mat[0][0] << endl;
                cout << "mat[1][1] : " << mat[1][1] << endl;
                cout << "mat[1][2] : " << mat[1][2] << endl;
                cout << "mat[0][2] : " << mat[0][2] << endl;
                cout << "mat[2][2] : " << mat[2][2] << endl;
         
                cout << "mat.determinant() : " << mat.determinant() << endl;
                cout << "mat.inverse() : " << mat.inverse() << endl;
                cout << "mat.transpose() : " << mat.tranpose() << endl;
         
                Cmat matoto=Cmat(4,5,3,4,5,6,7,8,9);
                cout << "matoto : " << matoto << endl;
                cout << "matoto.subdet(1,0) : " << matoto.subdet(1,0) << endl;
                cout << "matoto.det() : " << matoto.determinant() << endl;
         
                cout << "matoto.inverse() : " << matoto.inverse() << endl;
         
                cout << "A * A^-1 : " << matoto*matoto.inverse() << endl;
        */
        Cbasis obs(Cvec(0,0,-100),0,0,0);
        Cdvoxel voxel(Cvec(0,0,0),Cmat(0,1),1.,1e5,obs);

        cout << "voxel : " << voxel << endl;

        /*
        vector<Cvec> vert=voxel.vertobs;

        float r,phi,theta;

        for (int ii=0;ii<8;ii++) {
            cout << "v"<<ii<<" : "<< vert[ii] << endl;
            cvcoord(vert[ii],&r,&phi,&theta);
            cout << "r : " << r
            << ", phi : " << phi*RADEG
            << ", theta : " << theta*RADEG << endl;
        }
        */
        // test losobs2xxyy function
        float xx,yy;
        Cvec vvv(1,1,100);
        losobs2xxyy(vvv,xx,yy);
        cout << "xx : " << xx << endl;
        cout << "yy : " << yy << endl;

        vector<float> xy(2,0);
        cout << "xy.size() : " << xy.size() << endl;

        // compute LOS param for the 8 vertices
        //vector<vector<float> > vxy(8, vector<float>(2,0));

        /*
        vector<vector<float> > vij=voxelvert2iijj(vert,5.45415e-05,5.45415e-05,512,512);

        for (int i=0;i<8;i++) {
            cout << "vxy["<<i<<"] : " << vij[0][i] << " , " << vij[1][i] << endl;
        }
        */
        
        vector<float> L;
        L.push_back(2);
        L.push_back(2.6);
        L.push_back(4.6);
        L.push_back(3.2);

        vector<float>::const_iterator it = max_element(L.begin(),L.end());
        cout << "The largest element is " << *it << endl;

        Cvec va=ChangeBase(Cvec(1.,0.,0.),Cbasis(Cvec(0,0,0),0.,0.,0.),
                           Cbasis(Cvec(10,0,0),0.,0.,PI));
        
        cout << "va : " << va << endl;
        
        va+=Cvec(1,1,1);
        cout << "va : " << va << endl;

        vector<Cvec> vt;
        vt.push_back(Cvec(1,2,3));
        
        cout << "vt[0] : " << vt[0] << endl;
        
        vt[0]+=Cvec(1,2,3);
        
        cout << "vt[0] : " << vt[0] << endl;
        
        cout << "atan2(0,1) : " << atan2(double(0),double(1))*RADEG << endl;
        cout << "atan2(1,0) : " << atan2(double(1),double(0))*RADEG << endl;
        cout << "atan2(0,-1) : " << atan2(double(0),double(-1))*RADEG << endl;
        cout << "atan2(-1,0) : " << atan2(double(-1),double(0))*RADEG << endl;
        cout << "atan2(1,-1) : " << atan2(double(1),double(-1))*RADEG << endl;
        cout << "atan2(-1,-1) : " << atan2(double(-1),double(-1))*RADEG << endl;
        cout << "atan2(-1,1) : " << atan2(double(-1),double(1))*RADEG << endl;
        
        
        // ---- try write binary file
        /*
        fstream binFile;
        binFile.open("binfiletest.dat",ios::out | ios::binary);
        if (!binFile.is_open()) {
          cout << "Cannot open file for writting !" << endl;
          return 0;
        }
        
        unsigned int data=54;
        binFile.write((char *)(&data),sizeof(unsigned int)); 
        
        binFile.close();
        */
        
        Cbasis abs(Cvec(0,0,0),0,0,0);

        // ---- observer position definition
        Cbasis obs2;
        obs2=Cbasis(0.,0.,214.,0.,0.,0.);

        Cvec vlosabs;
        float rho;
        float crpix[2]={1023.5,1023.5};
        float pc[4]={1.,0.,0.,1.};
        float *fovpix;
        
        *fovpix=7.02980e-05; // rad/pix: cor2 A
        
        /*
        getimpactandlos(1023.5,0.,obs2,abs,2,0.,vlosabs,rho,crpix,fovpix);

        
        cout << endl << "vlosabs : " << vlosabs << endl;
        cout << "rho : " << rho << endl;
        */
        getimpactandlos(0.,1023.5,obs2,abs,2,0.,vlosabs,rho,crpix,fovpix,pc);

        
        cout << endl << "vlosabs : " << vlosabs << endl;
        cout << "rho : " << rho << endl;

        float xxxx=-(*fovpix)*1023.;
        getimpactandlos(-0.072,0.,obs2,abs,2,0.,vlosabs,rho);
        
        cout << endl << "vlosabs : " << vlosabs << endl;
        cout << "rho : " << rho << endl<<endl;

        float llr[3]={2.,-1.,10.};
        
        
        Cvec c=carrington2cart(PI-llr[0],llr[1],llr[2]);
        cout << "lon : " << llr[0] << ", lat : " << llr[1] << ", rad : " << llr[2] << endl;
        
        cout << "c : " << c << endl<<endl;
        
        float lonlatrad[3];
        cart2carrington(c,lonlatrad);
        lonlatrad[0]=PI-lonlatrad[0];
        
        cout << "lon : " << lonlatrad[0] << ", lat : " << lonlatrad[1] << ", rad : " << lonlatrad[2] << endl;
        
        cout << endl;
        
        // -- test rtcloud
        cout << "Test rtcloud..." << endl;
        Cvec Vlosobs(1.,1.,1.);
        float crpx[2];
        crpx[0]=255.5;crpx[1]=255.5;
        float fovpx=15./3600.;
        float pci[4]={1,0,0,1};
        int projtypecode=2;
        float pv2_1=5.;
        float xx0,yy0;
        
        losobs2xxyy(Vlosobs,crpx,&fovpx,pci,projtypecode,pv2_1,xx0,yy0);
        
        cout << "Vlosobs : " << Vlosobs << endl;
        cout << "atan2(0,0) : " << atan2(0.,0.) << endl;
        cout << "x,y : " << xx0 << " , " << yy0 << endl;
        
        
        cout << endl << "test roots of polynomial 2nd deg..." << endl;
        float root1,root2;
        float determin=get2ndDegPolyRoots(-1.,2.,3.,root1,root2);
        
        cout << "det : " << determin << endl;
        cout << "root1 : " << root1 << endl;
        cout << "root2 : " << root2 << endl;
        
        cout << "Chalut ! " << endl;
        
	cout << atan2(1.,0.5)*RADEG << endl;
	cout << atan2(0.,-1.)*RADEG << endl;
	cout << atan2(-1.,-1.)*RADEG << endl;



        dumpBuildInfo();
        
        
#ifdef HAVE_LIBCCFITS
cout << "CCfits enabled" << endl;
#else
cout << "No CCfits" << endl;
#endif
        
        
        return 1;
    }


    cout << "modelid : " << modelid << endl;

    cout << "Bujur Christophe !" << endl;

    Cvec v1(1.,2,3),v2(1.,1.,1),vr;
    Cmat m(1.,1);
    Cbasis abs(Cvec(0,0,0),0,0,0);
    Clos los(512,-20,20);

    float scalprod = pscal(v1,v2);
    vr = pvect(v1,v2);

    cout << "v1 : " << v1 << endl;
    cout << "v2 : " << v2 << endl;
    cout << "sp : " << scalprod << endl;
    cout << "vr : " << vr << endl;

    CModelBase *pmod;
    pmod=modelselect(modelid);

    float *pparam;
    pparam=NULL;
    pmod->initParam(pparam);


    float temperature=0;
    Cvec v(10,10,10);

    float ner=pmod->Density(v);

    cout << "Density is " << ner << " at " << v << endl;


    // ---- test getParam
    //int nbparam;
    //nbparam=pmod->getNbparam();
    //cout << "nb param : " << nbparam << endl;


    delete pmod;


    // ---- POI basis
    Cbasis poi(Cvec(0,0,0),10.*DTOR,0,0);

    cout << "u  : " << poi.u << endl;
    cout << "ui : " << poi.ui << endl;


    // test of Cnurbs

    vector<float> KnotVec;
    KnotVec.push_back(0);
    KnotVec.push_back(0);
    KnotVec.push_back(0);
    KnotVec.push_back(1);
    KnotVec.push_back(1);
    KnotVec.push_back(2);
    KnotVec.push_back(2);
    KnotVec.push_back(3);
    KnotVec.push_back(3);
    KnotVec.push_back(3);

    cout << "KnotVec : " ;
    for (unsigned i=0;i<KnotVec.size();i++)
        cout << KnotVec[i] << ",";
    cout << endl;

    vector<CControlPoint> cp(7);
    cp[0]=CControlPoint(4,0,0,1);
    cp[1]=CControlPoint(8,0,1,0.5);
    cp[2]=CControlPoint(6,4,2,1);
    cp[3]=CControlPoint(4,8,3,0.5);
    cp[4]=CControlPoint(2,4,2,1);
    cp[5]=CControlPoint(0,0,1,0.5);
    cp[6]=CControlPoint(4,0,0,1);

    Cnurbs nurbs(7,KnotVec,10,cp);

    //nurbs.printBasis(2);

    nurbs.printCurve();


    //vector<Cvec> curve;
    //curve.push_back(Cvec(0,0,0));
    //curve.push_back(Cvec(1,1,1));
    //curve.push_back(Cvec(2,2,2));


    //cout << "curve[0] : " << curve[0] << endl;
    //cout << "curve[1].v[0] : " << curve[1].v[0] << endl;

    float rrr,phi,theta;
    cvcoord(0,-1,0,&rrr,&phi,&theta);
    cout << "rrr : " << rrr << endl;
    cout << "phi : " << phi*RADEG << endl;
    cout << "theta : " << theta*RADEG << endl;


    Cbasis nps(Cvec(0,0,0),
               0,0,0,
               0,0,0.5);

    cout << "nps.u : " << nps.u << endl;
    cout << "nps.ui : " << nps.ui << endl;
    cout << "nps.o : " << nps.o << endl;


    // -- test commutativity
    Cvec vout;
    vout=v1*3;
    cout << "v1 * 3 : " << vout << endl;
    vout=3*v1;
    cout << "3 * v1 : " << vout << endl;

    // -- test call with array
    float tott[3],polt[3];
    getThomsonGeomFactor(v1.norm(),3.,tott[0],polt[0]);
    cout << "tott[0] : " << tott[0] << ", polt[0] : " << polt[0] << endl;
    getThomsonGeomFactor(v1.norm(),3.,tott[1],polt[1]);
    cout << "tott[1] : " << tott[1] << ", polt[1] : " << polt[1] << endl;

    cout << "Perfect !" << endl;

    
    
    // -- test quick raytrace
    rtparam fp;
    
    
    
    
    return 0;
}


/*
* $Log: main.cpp,v $
* Revision 1.12  2010-09-17 15:24:43  thernis
* Cast values in atan2 function to avoid compilation crash in Solaris
*
* Revision 1.11  2010-08-16 20:43:45  colaninn
* Added cstdlib and algorithm
*
* Revision 1.10  2009/02/09 20:51:11  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.9  2008/09/23 14:08:12  thernis
* - add new testing models, - implement integration in front of the instrument
*
* Revision 1.8  2007/07/20 21:45:38  thernis
* Show if CCfits was used or not
*
* Revision 1.7  2007/07/19 19:52:13  thernis
* Implement test for second degree polynomial root calculation testing
*
* Revision 1.6  2007/07/10 21:15:57  thernis
* Implement raytracing of a cloud of points.
*
* Revision 1.5  2007/05/14 17:19:40  thernis
* Add CVS id and log in all files
*
*/
