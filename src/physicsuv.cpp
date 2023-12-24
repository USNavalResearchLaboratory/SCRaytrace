
#include "physicsuv.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <rpc/xdr.h>
#include <string>
#include <sstream>
#include <cmath>

#include "config.h"
#include "rtmiscfunc.h"
#include "scene.h"


const unsigned int PhysicsUV::NBSAMP;
const unsigned int PhysicsUV::NBLINES;


PhysicsUV::PhysicsUV()
{
    physicsName="UV emission";
    wavelengthId=0;

    isgood=false;
//     char *ppathtoxdr;
    string pathtoxdr;
    // ---- get the path where emission files are located
//     ppathtoxdr=getenv("RT_UVEMISSIONPATH");
    pathtoxdr.assign(SCRAYTRACE_DATA_DIR);

//     if (ppathtoxdr!=NULL) pathtoxdr.assign(ppathtoxdr);
// 
//     if (pathtoxdr.empty())
//     {
//         cout << "The RT_UVEMISSIONPATH environment variable is not defined\nI just use the current directory" << endl;
//         pathtoxdr="~/work/cpp/scraytrace/data";
//     }

//pathtoxdr="/home/arnaud/work/cpp";

    // ---- load the emission vs temperature tables
    // -- yisel
    FILE *pFile;
    string filename=pathtoxdr+"/yisel.xdr";
    pFile = fopen (filename.c_str() , "r");
    if (pFile == NULL)
    {
        cout << "Error opening file "<< filename << endl;
        exit(1);
    }

    XDR xdrs;
    xdrstdio_create(&xdrs, pFile, XDR_DECODE);

    float a;

    for (unsigned int j=0;j<NBLINES;j++) for (unsigned int i=0;i<NBSAMP;i++)
        {
            if (!xdr_float(&xdrs, &a))
            {
                printf("Error parsing yisel.xdr file!\n");
                fclose(pFile);
                exit(1);
            }
            yisel[i+j*NBSAMP]=log10(a+1e-30);
        }

    xdr_destroy(&xdrs);
    fclose (pFile);


    // -- emis
    filename.clear();
    filename=pathtoxdr+"/emis.xdr";
    pFile = fopen (filename.c_str() , "r");
    if (pFile == NULL)
    {
        cout << "Error opening file " << filename << endl;
        exit(1);
    }

    xdrstdio_create(&xdrs, pFile, XDR_DECODE);

    for (unsigned int j=0;j<NBLINES;j++) for (unsigned int i=0;i<NBSAMP;i++)
        {
            if (!xdr_float(&xdrs, &a))
            {
                printf("Error parsing emis.xdr file!\n");
                fclose(pFile);
                exit(1);
            }
            emis[i+j*NBSAMP]=log10(a);
        }

    xdr_destroy(&xdrs);
    fclose (pFile);


    // -- ti
    filename.clear();
    filename=pathtoxdr+"/ti.xdr";
    pFile = fopen (filename.c_str() , "r");
    if (pFile == NULL)
    {
        cout << "Error opening file "<< filename << endl;
        exit(1);
    }

    xdrstdio_create(&xdrs, pFile, XDR_DECODE);

    for (unsigned int i=0;i<NBSAMP;i++)
    {
        if (!xdr_float(&xdrs, &a))
        {
            printf("Error parsing ti.xdr file!\n");
            fclose(pFile);
            exit(1);
        }
        ti[i]=a;
    }

    xdr_destroy(&xdrs);
    fclose (pFile);


    // -- teg
    filename.clear();
    filename=pathtoxdr+"/teg.xdr";
    pFile = fopen (filename.c_str() , "r");
    if (pFile == NULL)
    {
        cout << "Error opening file "<< filename << endl;
        exit(1);
    }

    xdrstdio_create(&xdrs, pFile, XDR_DECODE);

    for (unsigned int i=0;i<NBSAMP;i++)
    {
        if (!xdr_float(&xdrs, &a))
        {
            printf("Error parsing teg.xdr file!\n");
            fclose(pFile);
            exit(1);
        }
        teg[i]=a;
    }

    xdr_destroy(&xdrs);

    // -- set isgood flag to true since it looks like everything was alright
    isgood=true;

}


//! Return True if the tables have been loaded properly
bool PhysicsUV::IsGood()
{
    return isgood;
}


//! Return yisel
float PhysicsUV::getyisel(const unsigned int &i,const unsigned int &j)
{
    return yisel[i+j*NBSAMP];
}

//! Compute the emissivity function of the line and the temperature
float PhysicsUV::calcEmissivity(const unsigned int &waveid,const float &te)
{
    if (waveid > 3) return 0.;
    //unsigned int klineoffset=kline-1;
    float logte=log10(te);
    float yif=nearestneighbor1dinterp(logte,ti,NBSAMP,&yisel[waveid*NBSAMP]);
    float emiss=nearestneighbor1dinterp(logte,teg,NBSAMP,&emis[waveid*NBSAMP]);

    return(pow(float(10.),yif)*pow(float(10.),emiss));
}




bool PhysicsUV::computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout)
{
    float temperature;
    float ner=pparentscene->pmod->Density(ChangetoDensityCoord(pparentscene->modelposition,vs),temperature);
    float uv=calcEmissivity(wavelengthId,temperature);

    neout=ner;
    bpout=uv;
    btout=uv*ner*ner;

    return 0;
}


void PhysicsUV::getConstFactors(float &btf,float &bpf,float &nef, float rho)
{
    btf=pparentscene->los.ds;
    bpf=btf;
    nef=RSUN_CM*pparentscene->los.ds;
}
