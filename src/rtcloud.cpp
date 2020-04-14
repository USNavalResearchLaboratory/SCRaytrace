

#include <iostream>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "Cvec.h"
#include "Clos.h"
#include "CModelBase.h"
#include "Cbasis.h"
#include "raytrace.h"
#include "rtcloud.h"
using namespace std;

extern "C" int rtcloud ( int argc, void **argv )
{
    // ---- extract input parameters
    int *is= ( int* ) argv[0]; // -- x size of the CCD
    int *js= ( int* ) argv[1]; // -- y size of the CCD
    float *fovpix= ( float* ) argv[2]; // -- resolution of the pix, in rad
    float *obspos= ( float* ) argv[3]; // -- position of the observer
    float *obsang= ( float* ) argv[4]; // -- angle position of the observer
    float *nepos= ( float* ) argv[5]; // -- Position of the reference point of the Ne
    float *neang= ( float* ) argv[6]; // -- angle position of the Ne
    float *btot= ( float* ) argv[7]; // -- output image
    float *crpix= ( float* ) argv[8]; // -- position of the optical axis, in pix
    int *quiet= ( int* ) argv[9]; // -- flag for quiet mode
    float *hlonlat= ( float* ) argv[10]; // -- car lon lat height of the Ne
    float *protmat= ( float* ) argv[11]; // -- (out) rotation matrix for the Ne [2,2]
    float *obslonlat= ( float* ) argv[12]; // -- car lon lat height of the obs [3]
    int obslonlatflag=* ( ( int* ) argv[13] ); // -- flag if obslonlat must be used
    int projtypecode=* ( ( int* ) argv[14] ); // -- projection type
    float pv2_1=* ( ( float* ) argv[15] ); // -- pv2_1 mu parameter for the azp projection
    float *pc= ( float* ) argv[16]; // -- pc matrix [2,2]
    float *plist= ( float* ) argv[17]; // -- list of points [3,N]
    int nbp=* ( ( int* ) argv[18] ); // -- number of points in the list of points
    float *pci= ( float* ) argv[19]; // -- inverse of pc matrix [2,2]
    int flaglistout=* ( ( int* ) argv[20] ); // -- set to 1 is the user wants a the list of points
    float *plistout= ( float* ) argv[21]; // -- (out) list of pixel position given the list of points plist
    float flagclip=* ( ( int* ) argv[22] ); // -- set to 1 if you don't want to compute the points masked by the sun
    int *plistbehind= ( int* ) argv[23] ; // -- (out) 1 if the point is behind disk, else 0
    float *nerotcntr= ( float* ) argv[24]; // -- Position of the rotation center within Ne
    float *nerotang= ( float* ) argv[25]; // -- rotation angles within Ne
    float *netranslation= ( float* ) argv[26]; // -- Shift of the density model, within the model coordinate system
    int *nerotaxis= (int*) argv[27]; // -- Axis order for the rotation angles within the Ne

    // ---- get the params from IDL
    if ( *quiet == 0 )
    {
        dumpBuildInfo();
        cout << "In rtcloud..." << endl;
        cout << "is : " << *is << " , js : " << *js << endl;
        cout << "fovpix : " << *fovpix << endl;
        cout << "obspos : " << *obspos;
        for ( unsigned int i=1;i<3;i++ )
            cout << " , " << * ( obspos + i );
        cout << endl << "obsang : " << *obsang;
        for ( unsigned int i=1;i<3;i++ )
            cout << " , "  << * ( obsang + i );
        cout << endl << "nepos : " << *nepos;
        for ( unsigned int i=1;i<3;i++ )
            cout << " , "  << * ( nepos + i );
        cout << endl << "neang : " << *neang;
        for ( unsigned int i=1;i<3;i++ )
            cout << " , "  << * ( neang + i );
        cout << endl << "nerotcntr : " << *nerotcntr;
        for ( unsigned int i=1;i<3;i++ )
            cout << " , "  << * ( nerotcntr + i );
        cout << endl << "nerotang : " << *nerotang;
        for ( unsigned int i=1;i<3;i++ )
            cout << " , "  << * ( nerotang + i );
        cout << endl << "nerotaxis : " << *nerotaxis;
        for ( unsigned int i=1;i<3;i++ )
            cout << " , "  << * ( nerotaxis + i );
        cout << endl << "netranslation : " << *netranslation;
        for ( unsigned int i=1;i<3;i++ )
            cout << " , "  << * ( netranslation + i );
        cout << endl;
        cout << "Hlonlat : " << hlonlat[0] << " , " << hlonlat[1] << " , " << hlonlat[2] << endl;
        cout << "obslonlatflag : " << obslonlatflag << endl;
        cout << "obslonlat : " << *obslonlat;
        for ( unsigned int i=1;i<3;i++ )
            cout << " , "  << * ( obslonlat + i );
        cout << endl;
        cout << "projtype : " << projtypecode << endl;
        cout << "crpix : " << crpix[0] << " , " << crpix[1] << endl;
        cout << "pc [0,1,2,3] : " << pc[0];
        for ( unsigned int i=1;i<=3;i++ ) cout << " , " << pc[i];
        cout << endl;
        cout << "pv2_1 : " << pv2_1 << endl;
        
        
    }

    // -- array position ptr
    float *posbtot;
    posbtot=btot;

    float rho,r,s;

    // ---- absolute basis definition
    Cbasis abs ( Cvec ( 0,0,0 ),0,0,0 );

    // ---- observer position definition
    Cbasis obs;
    if ( obslonlatflag == 0 )
        obs=Cbasis ( Cvec ( obspos[0],obspos[1],obspos[2] ),obsang[0],obsang[1],obsang[2] );
    else
        obs=Cbasis ( obslonlat[0],obslonlat[1],obslonlat[2],obsang[0],obsang[1],obsang[2] );

    // ---- Ne position
    ModelPosition modpos;
    modpos.setBasis(Cbasis(Cvec(nepos[0],nepos[1],nepos[2]),
               neang[0], neang[1],neang[2],hlonlat[0],hlonlat[1],hlonlat[2]));
    modpos.rotation.setCenter(Cvec(nerotcntr[0],nerotcntr[1],nerotcntr[2]));
    modpos.rotation.setRotationPerAxis(nerotang[0],nerotaxis[0],nerotang[1],nerotaxis[1],nerotang[2],nerotaxis[2]);
    modpos.setTranslation(Cvec(netranslation[0],netranslation[1],netranslation[2]));

    // -- dump the rotation matrix in rotmat output variable
    for ( unsigned int i=0,k=0;i<3;i++ )
        for ( unsigned int j=0;j<3;j++ )
            protmat[k++]=modpos.modelbasis.u.m[j][i];
//            protmat[k++]=nps.u.m[j][i];

    // -- progression
    float progresspercent=0.2;
    int progressflag= ( int ) ( ( float ) nbp * progresspercent );
    float progresspass=progresspercent;

    // -- loop for all the points of the point list
    float ii,jj;
    for ( unsigned int i=0;i<nbp;i++ )
    {
        // -- print progression
        if ( ( ( *quiet ) == 0 ) && ( i > progressflag ) )
        {
            cout << progresspass*100 << "% " << endl;
            progresspass+=progresspercent;
            progressflag= ( int ) ( ( float ) nbp * progresspass );
        }

        // ---- compute position of the current point depending on the positioning of the Ne
        // -- compute pos of the point in the abs basis
        int pointpos=i*3;
        Cvec OMne=Cvec ( plist[pointpos],
                         plist[pointpos+1],
                         plist[pointpos+2] );
        // -- compute the LOS passing by that point

        Cvec Vlosabs=ChangeDensityCoordtoAbs(modpos,OMne)-obs.o;

        Cvec Vlosobs=obs.u * Vlosabs;

        // compute impact parameters
        float rho=psldist ( obs.o,Vlosabs,abs.o );
        if ( rho < 1. )
        {
            Cvec Qlos=orthoproj ( obs.o,Vlosabs,abs.o );
            if ( Vlosabs.norm() > ( ( Qlos-obs.o ).norm()-sqrt ( fabs(0.99-Qlos.magsqr()) ) ) )
            {
                plistbehind[i]=1;
                if ( flagclip != 0 ) continue;
            }
        }


        // -- compute pixel pos corresponding to that LOS
        losobs2xxyy ( Vlosobs,crpix,fovpix,pci,projtypecode,pv2_1,ii,jj );

        if ( flaglistout == 0 )
        {
            // -- place the calculated point in the output image
            unsigned int iii= ( unsigned int ) ( ii+0.5 );
            unsigned int jjj= ( unsigned int ) ( jj+0.5 );

            // -- set the image point to 1 if within the boundaries of the image
            if ( iii >= 0 && iii < ( *is ) && jjj>=0 && jjj < ( *js ) ) posbtot[iii+ ( *is ) *jjj]=1.;
        }
        else
        {
            // -- put the calculated points in the list
            plistout[0+2*i]=ii;plistout[1+2*i]=jj;
        }

    }

    return 1;

}
