
#include <fstream>
#include <sstream>
#include <algorithm>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_SUNMATH_H
#include <sunmath.h>
#endif
#include "raytrace.h"
#include "rtmiscfunc.h"
#include "rtthread.h"

using namespace std;


//! Main raytracing routine in cartesian coordinates
int raytracemain(rtparam fp) {
    // extract parameters from the structure
    int *is;
    is=fp.pis;
    int *js;
    js=fp.pjs;
    float *fovpix;
    fovpix=fp.pfovpix;
    float *obspos;
    obspos=fp.pobspos;
    float *obsang;
    obsang=fp.pobsang;
    float *nepos;
    nepos=fp.pnepos;
    float *neang;
    neang=fp.pneang;
    int *losnbp;
    losnbp=fp.plosnbp;
    float *losrange;
    losrange=fp.plosrange;
    int *pmodelid;
    pmodelid=fp.pmodelid;
    float *btot;
    btot=fp.pbtot;
    float *bpol;
    bpol=fp.pbpol;
    float *netot;
    netot=fp.pnetot;
    float *pmodparam;
    pmodparam=fp.pmparam;
    float *crpix;
    crpix=fp.pcrpix;
    float *rhoim;
    rhoim=fp.prhoim;
    float *mmlon;
    mmlon=fp.pmmlon;
    float *mmlat;
    mmlat=fp.pmmlat;
    float *rrr;
    rrr=fp.prrr;
    int *ppofinteg;
    ppofinteg=fp.ppofinteg;
    int *pfrontinteg;
    pfrontinteg=fp.pfrontinteg;
    int *quiet;
    quiet=fp.pquiet;
    int flagneonly;
    flagneonly=*(fp.pneonly);
    int *proi;
    proi=fp.proi;
    float *ppoiang;
    ppoiang=fp.ppoiang;
    float occrad;
    occrad=*(fp.poccrad);
    float adapthres;
    adapthres=*(fp.padapthres);
    int maxsubdiv;
    maxsubdiv=*(fp.pmaxsubdiv);
    float limbdark;
    limbdark=*(fp.plimbdark);
    float *protmat;
    protmat=fp.protmat;
    int i,j,k;
    int obslonlatflag;
    obslonlatflag=*(fp.obslonlatflag);
    int projtypecode=*(fp.projtypecode);
    float pv2_1=*(fp.pv2_1);
    unsigned int uvinteg=fp.uvinteg;    // kept for backward compatibility but deprecated
    float disttofracmax=fp.disttofracmax;
    float *nerotcntr;
    nerotcntr=fp.pnerotcntr;
    float *nerotang;
    nerotang=fp.pnerotang;
    float *netranslation;
    netranslation=fp.pnetranslation;
    int *nerotaxis;
    nerotaxis=fp.pnerotaxis;
    float *losDepthIn;
    losDepthIn = fp.losDepthIn;
    float *losDepthOut;
    losDepthOut = fp.losDepthOut;
    int evalDepth;
    evalDepth = fp.evalDepth;
    
    // ---- get the params from IDL
    if (*quiet == 0) {
    dumpBuildInfo();
        cout << "Model : " << *pmodelid << endl;
        cout << "is : " << *is << " , js : " << *js << endl;
        cout << "fovpix : " << *fovpix << endl;
        cout << "obspos : " << *obspos;
        for(i=1;i<3;i++)
            cout << " , " << *(obspos + i);
        cout << endl << "obsang : " << *obsang;
        for(i=1;i<3;i++)
            cout << " , "  << *(obsang + i);
        cout << endl << "nepos : " << *nepos;
        for(i=1;i<3;i++)
            cout << " , "  << *(nepos + i);
        cout << endl << "neang : " << *neang;
        for(i=1;i<3;i++)
            cout << " , "  << *(neang + i);
        cout << endl << "nerotcntr : " << *nerotcntr;
        for(i=1;i<3;i++)
            cout << " , "  << *(nerotcntr + i);
        cout << endl << "nerotang : " << *nerotang;
        for(i=1;i<3;i++)
            cout << " , "  << *(nerotang + i);
        cout << endl << "nerotaxis : " << *nerotaxis;
        for(i=1;i<3;i++) cout << " , "  << *(nerotaxis + i);
        cout << endl << "netranslation : " << *netranslation;
        for(i=1;i<3;i++)
        cout << " , "  << *(netranslation + i);
        cout << endl << "losnbp : " << *losnbp << endl;
        cout << "losrange : " << *losrange;
        for(i=1;i<2;i++)
            cout << " , "  << *(losrange + i);
        cout << endl;
        cout << "Hlonlat : " << fp.Hlonlat[0] << " , " << fp.Hlonlat[1] << " , " << fp.Hlonlat[2] << endl;
        cout << "occrad : " << occrad << endl;
        cout << "adapthres : " << adapthres << endl;
        cout << "maxsubdiv : " << maxsubdiv << endl;
        cout << "limbdark : " << limbdark << endl;
        cout << "obslonlatflag : " << obslonlatflag << endl;
        cout << "projtype : " << projtypecode << endl;
//         cout << "uvinteg : " << uvinteg << endl;
        cout << "disttofracmax : " << disttofracmax << endl;
    }

    // ---- select the density model requested by the user
    CModelBase *pmod;
    pmod=modelselect(*pmodelid);
    pmod->initParam(pmodparam);
    
    // -- array position ptr
    float *posbtot,*posbpol,*posne,*posrhoim,*posrrr,*plosDepthIn,*plosDepthOut;
    posbtot=btot;
    posbpol=bpol;
    posne=netot;
    posrhoim=rhoim;
    posrrr=rrr;
    int *posroi;
    posroi=proi;
    plosDepthIn=losDepthIn;
    plosDepthOut=losDepthOut;

    *plosDepthIn = 58;
    
    int pofinteg=*ppofinteg;
    int frontinteg=*pfrontinteg;

    // -- Thomson scattering constants
    float u=limbdark; // -- limb darkening: in prevision that user can change it
    float constfactor=constfactor(u);

    float xstep=*fovpix;
    float ystep=*fovpix;
    float xstart=-xstep*(*js-1)/2.; // is and js are the size of the output image in the image usual coordinates
    float ystart=-ystep*(*is-1)/2.;

    float xx,yy;

    float rho,r,s;
    float ner,cosomega,sinomega;
    float cossquareomega,sinsquareomega;
    float rhooverr;

    float a,b,c,d;
    float logterm;
    float polterm;

    // ---- absolute basis definition
    Cbasis abs(Cvec(0,0,0),0,0,0);

    // ---- observer position definition
    Cbasis obs;
    if (obslonlatflag == 0) 
    obs=Cbasis(Cvec(obspos[0],obspos[1],obspos[2]),obsang[0],obsang[1],obsang[2]);   
    else 
    obs=Cbasis(fp.obslonlat[0],fp.obslonlat[1],fp.obslonlat[2],obsang[0],obsang[1],obsang[2]);
    
    // ---- Ne position
    //    Cbasis nps(Cvec(nepos[0],nepos[1],nepos[2]),
    //               neang[0], neang[1], neang[2], fp.Hlonlat[0], fp.Hlonlat[1], fp.Hlonlat[2], Cvec(neshift[0], neshift[1], neshift[2]));
    //Cbasis nps(Cvec(nepos[0],nepos[1],nepos[2]),
    //           neang[0], neang[1],neang[2],fp.Hlonlat[0],fp.Hlonlat[1],fp.Hlonlat[2]);

    //Cbasis necoord2(Cvec(nerotcntr[0],nerotcntr[1],nerotcntr[2]),nerotang[0],nerotang[1],nerotang[2]);

    ModelPosition modpos;
    modpos.setBasis(Cbasis(Cvec(nepos[0],nepos[1],nepos[2]),
            neang[0], neang[1],neang[2],fp.Hlonlat[0],fp.Hlonlat[1],fp.Hlonlat[2]));
    modpos.rotation.setCenter(Cvec(nerotcntr[0],nerotcntr[1],nerotcntr[2]));
    modpos.rotation.setRotationPerAxis(nerotang[0],nerotaxis[0],nerotang[1],nerotaxis[1],nerotang[2],nerotaxis[2]);
    modpos.setTranslation(Cvec(netranslation[0],netranslation[1],netranslation[2]));

    // -- dump the rotation matrix in rotmat output variable
    for (unsigned int i=0,k=0;i<3;i++)
        for (unsigned int j=0;j<3;j++)
            protmat[k++]=modpos.modelbasis.u.m[j][i];

    // ---- POI basis (Plane Of Integration)
    Cbasis poi(Cvec(0,0,0),ppoiang[0],ppoiang[1],ppoiang[2]);//,fp.Hlonlat[0],fp.Hlonlat[1],fp.Hlonlat[2]);

    // ---- force integration around POI if orientation is present
    if (ppoiang[0]!=0 || ppoiang[1]!=0 || ppoiang[2]!=0)
        pofinteg=1;

    // -- LOS integration parameter definition
    Clos los(*losnbp,*losrange,*(losrange+1));

    // -- LOS vector in absolute basis
    Cvec vlosabs,vlosstep;

    // -- Sun cntr ortho proj on the LOS
    Cvec qlos;

    // -- current LOS point position
    Cvec vs,vslospos;

    // -- compute Sun center position on the image
    Cvec vobs2sun=obs.u * (-obs.o);
    if (*quiet == 0)
        cout << "vobs2sun : " << vobs2sun << endl;

    float phicntr=PI/2,thetacntr=0.;
    if (vobs2sun.v[2] != 0)
    thetacntr=applyprojection(atan(sqrt(sqr(vobs2sun.v[0])+sqr(vobs2sun.v[1]))/vobs2sun.v[2]),projtypecode,pv2_1);

    if (vobs2sun.v[1] !=0) {
        phicntr=atan(vobs2sun.v[0]/vobs2sun.v[1]);
        if (vobs2sun.v[1] < 0)
            phicntr+=PI;
    } else if (vobs2sun.v[0] < 0)
        phicntr=3*PI/2;

    // -- compute the Sun center in pix if not passed
    crpix[0]=thetacntr*cos(phicntr)/ystep+float(*is-1)/2;
    crpix[1]=thetacntr*sin(phicntr)/xstep+float(*js-1)/2;

    if (*quiet == 0) {
        cout << "phi   : " << phicntr << endl;
        cout << "theta : " << thetacntr << endl;

        cout << "crpix[0] : " << crpix[0] << ", ";
        cout << "crpix[1] : " << crpix[1] << endl;
    }

    // -- progression
    float progresspercent=0.2;
    int progressflag=(int) ((float) *js * progresspercent);
    float progresspass=progresspercent;

    // X axis : vertical, up: pix subscript j
    // Y axis : horizontal, right: pix subscript i

    // ---- compute min and max longitudes and latitudes
    *mmlat=xstart;
    *(mmlat+1)=xstart+(*js)*xstep;
    *mmlon=ystart;
    *(mmlon+1)=ystart+(*is)*ystep;

    xx=xstart;
    for(j=0;j<(*js);j++,xx+=xstep) {

        // -- print progression
        if (((*quiet) == 0) && (j > progressflag)) {
            cout << progresspass*100 << "% " << endl;
            progresspass+=progresspercent;
            progressflag=(int) ((float) *js * progresspass);
        }
        yy=ystart;
        for(i=0; i<(*is); i++, yy+=ystep, posbtot++, posbpol++, posne++, plosDepthIn++, plosDepthOut++) {
            // -- compute the LOS direction and impact parameter
            getimpactandlos(xx,yy,obs,abs,projtypecode,pv2_1,vlosabs,rho);
            //cout << "vlosabs : " << vlosabs << endl;

            // ---- fill out distance maps
            // -- impact parameter
            *(posrhoim++)=rho;
            // -- dist intersection LOS - plane of the sky containing the Sun cntr
            vslospos=obs.o+vlosabs*(-obs.o.v[2]/vlosabs.v[2]);
            *(posrrr++)=vslospos.norm();

            // -- skip integration if LOS within the occulter
            if (rho < occrad)
                continue;

            // -- skip integration if not in the roi
            if (*(posroi++)==0)
                continue;

    // ---- LOS integration
    if (adapthres <=0) {
        losinteg(pofinteg,frontinteg,obs,vlosabs,abs,poi,los,rho,pmod,pmodparam,modpos,u,constfactor,flagneonly,posbtot,posbpol,posne,plosDepthIn,plosDepthOut,evalDepth);
        } else {
        losintegadaptstep(pofinteg,frontinteg,obs,vlosabs,abs,poi,los,rho,pmod,pmodparam,modpos,u,constfactor,flagneonly,adapthres,maxsubdiv,posbtot,posbpol,posne);
        }
    // ---- find distance to half total b max
    if (disttofracmax > 0.) {
        losintegdisttofracmax(pofinteg,frontinteg,obs,vlosabs,abs,poi,los,rho,pmod,pmodparam,modpos,u,constfactor,flagneonly,posbtot,posbpol,posne,disttofracmax);
        }
    }
    }
    
    if (((*quiet) <= 1) ) {
        cout << "100% " << endl;
    }

    delete pmod;
    return EXIT_SUCCESS;
}



// ---- wrapper called by IDL
extern "C" int raytracewl(int argc, void **argv) {

    rtparam fp;

    fp.pis=(int*) argv[0];
    fp.pjs=(int*) argv[1];
    fp.pfovpix=(float*) argv[2];
    fp.pobspos=(float*) argv[3];
    fp.pobsang=(float*) argv[4];
    fp.pnepos=(float*) argv[5];
    fp.pneang=(float*) argv[6];
    fp.plosnbp=(int*) argv[7];
    fp.plosrange=(float*) argv[8];
    fp.pmodelid=(int*) argv[9];
    fp.pbtot=(float*) argv[10];
    fp.pbpol=(float*) argv[11];
    fp.pnetot=(float*) argv[12];
    fp.pmparam=(float*) argv[13];
    fp.pcrpix=(float*) argv[14];
    fp.prhoim=(float*) argv[15];
    fp.pmmlon=(float*) argv[16];
    fp.pmmlat=(float*) argv[17];
    fp.prrr=(float*) argv[18];
    fp.ppofinteg=(int*) argv[19];
    fp.pquiet=(int*) argv[20];
    fp.pneonly=(int*) argv[21];
    fp.proi=(int*) argv[22];
    fp.ppoiang=(float*) argv[23];
    fp.Hlonlat=(float*) argv[24];
    fp.poccrad=(float*) argv[25];
    fp.padapthres=(float*) argv[26];
    fp.pmaxsubdiv=(int*) argv[27];
    fp.plimbdark=(float*) argv[28];
    fp.protmat=(float*) argv[29];
    fp.obslonlat=(float*) argv[30];
    fp.obslonlatflag=(int*) argv[31];
    fp.projtypecode=(int*) argv[32];
    fp.pv2_1=(float*) argv[33];
    fp.pfrontinteg=(int*) argv[34];
    fp.uvinteg=*((unsigned int*) argv[35]); // kept for backward compatibility but deprecated
    fp.disttofracmax=*((float*) argv[36]);
    fp.pnerotcntr=(float*) argv[37];
    fp.pnerotang=(float*) argv[38];
    fp.pnetranslation=(float*) argv[39];
    fp.pnerotaxis=(int*) argv[40];
    fp.losDepthIn=(float*) argv[41];
    fp.losDepthOut=(float*) argv[42];
    fp.evalDepth=*((int*) argv[43]);
    
cout << "Yo!" << endl;

    // -- call the raytracing routine
    int rtxitstat=raytracemain(fp);
    return rtxitstat;
}


//! Raytracing using WCS coordinate system
extern "C" int rtraytracewcs(int argc, void **argv) {
// ---- extract input parameters

int *is=(int*) argv[0];
int *js=(int*) argv[1];
float *fovpix=(float*) argv[2];
float *obspos=(float*) argv[3];
float *obsang=(float*) argv[4];
float *nepos=(float*) argv[5];
float *neang=(float*) argv[6];
int *losnbp=(int*) argv[7];
float *losrange=(float*) argv[8];
int *pmodelid=(int*) argv[9];
float *btot=(float*) argv[10];
float *bpol=(float*) argv[11];
float *netot=(float*) argv[12];
float *pmodparam=(float*) argv[13];
float *crpix=(float*) argv[14];
float *rhoim=(float*) argv[15];
float *mmlon=(float*) argv[16];
float *mmlat=(float*) argv[17];
float *rrr=(float*) argv[18];
int *ppofinteg=(int*) argv[19];
int *quiet=(int*) argv[20];
int flagneonly=*((int*) argv[21]);
int *proi=(int*) argv[22];
float *ppoiang=(float*) argv[23];
float *hlonlat=(float*) argv[24];
float occrad=*((float*) argv[25]);
float adapthres=*((float*) argv[26]);
int maxsubdiv=*((int*) argv[27]);
float limbdark=*((float*) argv[28]);
float *protmat=(float*) argv[29];
float *obslonlat=(float*) argv[30];
int obslonlatflag=*((int*) argv[31]);
int projtypecode=*((int*) argv[32]);
float pv2_1=*((float*) argv[33]);
float *pc=(float*) argv[34];
int *pfrontinteg=(int*) argv[35];
unsigned int uvinteg=*((unsigned int*) argv[36]); // kept for backward compatibility but deprecated
float *nerotcntr=(float*) argv[37];
float *nerotang=(float*) argv[38];
float *netranslation=(float*) argv[39];
int *nerotaxis=(int*) argv[40];
float *losDepthIn=(float*) argv[41];
float *losDepthOut=(float*) argv[42];
int evalDepth=*((int*) argv[43]);

int i,j,k;

    // ---- get the params from IDL
if (*quiet == 0) {
    dumpBuildInfo();
    cout << "In rtraytracewcs..." << endl;
    cout << "Model : " << *pmodelid << endl;
    cout << "is : " << *is << " , js : " << *js << endl;
    cout << "fovpix : " << *fovpix << endl;
    cout << "obspos : " << *obspos;
    for(i=1;i<3;i++)
    cout << " , " << *(obspos + i);
    cout << endl << "obsang : " << *obsang;
    for(i=1;i<3;i++)
    cout << " , "  << *(obsang + i);
    cout << endl << "nepos : " << *nepos;
    for(i=1;i<3;i++)
    cout << " , "  << *(nepos + i);
    cout << endl << "neang : " << *neang;
    for(i=1;i<3;i++)
    cout << " , "  << *(neang + i);
    cout << endl << "nerotcntr : " << *nerotcntr;
    for(i=1;i<3;i++)
    cout << " , "  << *(nerotcntr + i);
    cout << endl << "nerotang : " << *nerotang;
    for(i=1;i<3;i++)
    cout << " , "  << *(nerotang + i);
    cout << endl << "nerotaxis : " << *nerotaxis;
    for(i=1;i<3;i++)
    cout << " , "  << *(nerotaxis + i);
    cout << endl << "netranslation : " << *netranslation;
    for(i=1;i<3;i++)
    cout << " , "  << *(netranslation + i);
    cout << endl << "losnbp : " << *losnbp << endl;
    cout << "losrange : " << *losrange;
    for(i=1;i<2;i++)
    cout << " , "  << *(losrange + i);
    cout << endl;
    cout << "Hlonlat : " << hlonlat[0] << " , " << hlonlat[1] << " , " << hlonlat[2] << endl;
    cout << "occrad : " << occrad << endl;
    cout << "adapthres : " << adapthres << endl;
    cout << "maxsubdiv : " << maxsubdiv << endl;
    cout << "limbdark : " << limbdark << endl;
    cout << "obslonlatflag : " << obslonlatflag << endl;
    cout << "projtype : " << projtypecode << endl;
//     cout << "uvinteg : " << uvinteg << endl;
    cout << "crpix : " << crpix[0] << " , " << crpix[1] << endl;
    cout << "pc [0,1,2,3] : " << pc[0];
    for (i=1;i<=3;i++) cout << " , " << pc[i];
    cout << endl;
    
    cout << "evalDepth : " << evalDepth << endl;
    
}

// ---- select the density model requested by the user
CModelBase *pmod;
pmod=modelselect(*pmodelid);
pmod->initParam(pmodparam);
    
// -- array position ptr
float *posbtot,*posbpol,*posne,*posrhoim,*posrrr,*plosDepthIn,*plosDepthOut;
posbtot=btot;
posbpol=bpol;
posne=netot;
posrhoim=rhoim;
posrrr=rrr;
int *posroi;
posroi=proi;
plosDepthIn=losDepthIn;
plosDepthOut=losDepthOut;


int pofinteg=*ppofinteg;
int frontinteg=*pfrontinteg;

// -- Thomson scattering constants
float u=limbdark; // -- limb darkening: in prevision that user can change it
float constfactor=constfactor(u);  

float rho,r,s;
float ner,cosomega,sinomega;
float cossquareomega,sinsquareomega;
float rhooverr;

float a,b,c,d;
float logterm;
float polterm;

    // ---- absolute basis definition
Cbasis abs(Cvec(0,0,0),0,0,0);

    // ---- observer position definition
Cbasis obs;
if (obslonlatflag == 0) 
    obs=Cbasis(Cvec(obspos[0],obspos[1],obspos[2]),obsang[0],obsang[1],obsang[2]);   
else 
    obs=Cbasis(obslonlat[0],obslonlat[1],obslonlat[2],obsang[0],obsang[1],obsang[2]);
    
    // ---- Ne position
    //    Cbasis nps(Cvec(nepos[0],nepos[1],nepos[2]),
    //               neang[0], neang[1],neang[2],hlonlat[0],hlonlat[1],hlonlat[2],Cvec(neshift[0],neshift[1],neshift[2]));
    //     Cbasis nps(Cvec(nepos[0],nepos[1],nepos[2]),
    //                neang[0], neang[1],neang[2],hlonlat[0],hlonlat[1],hlonlat[2]);
    // 
    //     Cbasis necoord2(Cvec(nerotcntr[0],nerotcntr[1],nerotcntr[2]),nerotang[0],nerotang[1],nerotang[2]);

    ModelPosition modpos;
    modpos.setBasis(Cbasis(Cvec(nepos[0],nepos[1],nepos[2]),
            neang[0], neang[1],neang[2],hlonlat[0],hlonlat[1],hlonlat[2]));
    modpos.rotation.setCenter(Cvec(nerotcntr[0],nerotcntr[1],nerotcntr[2]));
    modpos.rotation.setRotationPerAxis(nerotang[0],nerotaxis[0],nerotang[1],nerotaxis[1],nerotang[2],nerotaxis[2]);

    //    modpos.setRotation(Cbasis(Cvec(nerotcntr[0],nerotcntr[1],nerotcntr[2]),nerotang[0],nerotang[1],nerotang[2]));
    modpos.setTranslation(Cvec(netranslation[0],netranslation[1],netranslation[2]));

// -- dump the rotation matrix in rotmat output variable
for (unsigned int i=0,k=0;i<3;i++)
    for (unsigned int j=0;j<3;j++)
    protmat[k++]=modpos.modelbasis.u.m[j][i];
//      protmat[k++]=nps.u.m[j][i];

    // ---- POI basis
Cbasis poi(Cvec(0,0,0),ppoiang[0],ppoiang[1],ppoiang[2]);

    // ---- force integration around POI if orientation is present
if (ppoiang[0]!=0 || ppoiang[1]!=0 || ppoiang[2]!=0)
    pofinteg=1;


    // -- LOS integration parameter definition
Clos los(*losnbp,*losrange,*(losrange+1));

    // -- LOS vector in absolute basis
Cvec vlosabs,vlosstep;

    // -- Sun cntr ortho proj on the LOS
Cvec qlos;

    // -- current LOS point position
Cvec vs,vslospos;

    // -- progression
    float progresspercent=0.2;
    int progressflag=(int) ((float) *js * progresspercent);
    float progresspass=progresspercent;
        
    for(j=0;j<(*js);j++) {

        // -- print progression
    if (((*quiet) <= 1) && (j > progressflag)) {
        cout << progresspass*100 << "% " << endl;
        progresspass+=progresspercent;
        progressflag=(int) ((float) *js * progresspass);
    }

    for(i=0;i<(*is);i++, posbtot++, posbpol++, posne++, plosDepthIn++, plosDepthOut++) {
        // -- compute the LOS direction and impact parameter
        
        getimpactandlos(float(i),float(j),obs,abs,projtypecode,pv2_1,vlosabs,rho,crpix,fovpix,pc);
        
            // ---- fill out distance maps
            // -- impact parameter
        *(posrhoim++)=rho;
            // -- dist intersection LOS - plane of the sky containing the Sun cntr
        vslospos=obs.o+vlosabs*(-obs.o.v[2]/vlosabs.v[2]);
        *(posrrr++)=vslospos.norm();

            // -- skip integration if LOS within the occulter
        if (rho < occrad)
        continue;

            // -- skip integration if not in the roi
        if (*(posroi++)==0)
        continue;

    // ---- LOS integration
    if (adapthres <=0) {
        losinteg(pofinteg,frontinteg,obs,vlosabs,abs,poi,los,rho,pmod,pmodparam,modpos,u,constfactor,flagneonly,posbtot,posbpol,posne,plosDepthIn,plosDepthOut,evalDepth);
        } else {
        losintegadaptstep(pofinteg,frontinteg,obs,vlosabs,abs,poi,los,rho,pmod,pmodparam,modpos,u,constfactor,flagneonly,adapthres,maxsubdiv,posbtot,posbpol,posne);
        }
    }
    }

    if (((*quiet) <= 1) ) {
        cout << "100% " << endl;
    }
    
    delete pmod;
    return EXIT_SUCCESS;

}


//! raytracing program that integrates only a cirular profile
extern "C" int rtwlcirc(int argc, void **argv) {

    int i,j,k;

    // ---- get the params from IDL
    int *is,*js;
    float *fovpix;
    float *obspos,*obsang;
    float *nepos,*neang;
    int *losnbp;
    float *losrange;
    int *pmodelid;
    float *btot;
    float *bpol;
    float *netot;
    float *pmodparam;
    float *crpix; //!<sun center position (returned)
    int *nbang;
    float *radius;
    int *quiet;
    int flagneonly=0; // force brightness calculation for the moment
    float ppoiang[3]={0,0,0}; // POI orientation set to 0 for the moment
    float *Hlonlat;
    float *pcntr;
    float *plimbdark;
    float *pobslonlat;
    int *pobslonlatflag;

    is=(int*) argv[0];
    js=(int*) argv[1];
    fovpix=(float*) argv[2];
    obspos=(float*) argv[3];
    obsang=(float*) argv[4];
    nepos=(float*) argv[5];
    neang=(float*) argv[6];
    losnbp=(int*) argv[7];
    losrange=(float*) argv[8];
    pmodelid=(int*) argv[9];
    btot=(float*) argv[10];
    bpol=(float*) argv[11];
    netot=(float*) argv[12];
    pmodparam=(float*) argv[13];
    crpix=(float*) argv[14];
    nbang=(int*) argv[15];
    radius=(float*) argv[16];
    quiet=(int*) argv[17];
    Hlonlat=(float*) argv[18];
    pcntr=(float*) argv[19];
    plimbdark=(float*) argv[20];
    pobslonlat=(float*) argv[21];
    pobslonlatflag=(int*) argv[22];
        
    if (*quiet == 0) {
        cout << "Model : " << *pmodelid << endl;
        cout << "is : " << *is << " , js : " << *js << endl;
        cout << "fovpix : " << *fovpix << endl;
        cout << "obspos : " << *obspos;
        for(i=1;i<3;i++)
            cout << " , " << *(obspos + i);
        cout << endl << "obsang : " << *obsang;
        for(i=1;i<3;i++)
            cout << " , "  << *(obsang + i);
        cout << endl << "nepos : " << *nepos;
        for(i=1;i<3;i++)
            cout << " , "  << *(nepos + i);
        cout << endl << "neang : " << *neang;
        for(i=1;i<3;i++)
            cout << " , "  << *(neang + i);
        cout << endl << "losnbp : " << *losnbp << endl;
        cout << "losrange : " << *losrange;
        for(i=1;i<2;i++)
            cout << " , "  << *(losrange + i);
        cout << endl;
        cout << "Circular profile integration" << endl;
        cout << "Nb angules : " <<  *nbang << endl;
        cout << "Circular profile radius : " << *radius << endl;
        cout << "cntr : " << pcntr[0] << " , " << pcntr[1] << endl;
        cout << "Limb Darkening : " << *plimbdark << endl;
        cout << "obslonlatflag : " << *pobslonlatflag << endl;
    }

    // ---- select the density model requested by the user
    CModelBase *pmod;
    pmod=modelselect(*pmodelid);
    pmod->initParam(pmodparam);

    // -- array position ptr
    float *posbtot,*posbpol,*posne;
    posbtot=btot;
    posbpol=bpol;
    posne=netot;
    float *plosDepthIn,*plosDepthOut;
    
    int pofinteg=0;
    int frontinteg=0;

    // -- Thomson scattering constants
    float u=*plimbdark; // limb darkening
    float constfactor=constfactor(u);

    float radpix=(*radius) * (*fovpix);
    float astart=0.;
    float astep=DTOR*360./(*nbang);

    float ang;
    float xx,yy;

    float rho,r,s;
    float ner,cosomega,sinomega;
    float cossquareomega,sinsquareomega;
    float rhooverr;
    // -- temperature: not useful in Thomson scattering but declared
    //    for compatibility of the models with the radio raytracing
    //float temperature;

    float a,b,c,d;
    float logterm;
    float polterm;

    // ---- absolute basis definition
    Cbasis abs(Cvec(0,0,0),0,0,0);

    // ---- observer position definition
    Cbasis obs;
    if (*pobslonlatflag == 0) 
    obs=Cbasis(Cvec(obspos[0],obspos[1],obspos[2]),obsang[0],obsang[1],obsang[2]);   
    else 
    obs=Cbasis(pobslonlat[0],pobslonlat[1],pobslonlat[2],obsang[0],obsang[1],obsang[2]);

    // ---- Ne position
    Cbasis nps(Cvec(*nepos,*(nepos+1),*(nepos+2)),
            *neang, *(neang+1), *(neang+2),Hlonlat[0],Hlonlat[1],Hlonlat[2]);

    ModelPosition modpos;
    modpos.setBasis(nps);

    // ---- POI basis
    Cbasis poi(Cvec(0,0,0),ppoiang[0],ppoiang[1],ppoiang[2]);

    // -- LOS integration parameter definition
    Clos los(*losnbp,*losrange,*(losrange+1));

    // -- LOS vector in absolute basis
    Cvec vlosabs,vlosstep;

    // -- Sun cntr ortho proj on the LOS
    Cvec qlos;

    // -- current LOS point position
    Cvec vs;

    // -- compute Sun center position on the images
    Cvec suncntrobs;
    suncntrobs=orthoprojpointplane(Cvec(0,0,0),obs.u.row(3),-obs.o);
    Cvec vobs2sun=obs.u * (-obs.o);
    cout << "vobs2sun : " << vobs2sun << endl;

    float phicntr=PI/2,thetacntr=0.;
    if (vobs2sun.v[2] != 0)
        thetacntr=atan(sqrt(sqr(vobs2sun.v[0])+sqr(vobs2sun.v[1]))/vobs2sun.v[2]);

    if (vobs2sun.v[1] !=0) {
        phicntr=atan(vobs2sun.v[0]/vobs2sun.v[1]);
        if (vobs2sun.v[1] < 0)
            phicntr+=PI;
    } else if (vobs2sun.v[0] < 0)
        phicntr=3*PI/2;

    cout << "phi   : " << phicntr << endl;
    cout << "theta : " << thetacntr << endl;

    crpix[0]=thetacntr*cos(phicntr)/(*fovpix)+float(*is-1)/2;
    crpix[1]=thetacntr*sin(phicntr)/(*fovpix)+float(*js-1)/2;

    cout << "crpix[0] : " << crpix[0] << ", ";
    cout << "crpix[1] : " << crpix[1] << endl;

    // -- center shift of the circular profile
    //    to match the Sun center
    //  float solar_r=fabs(atan(1./obs.o.norm()));
    //  float xshft=suncntrobs.v[0]*solar_r; //-xstart;
    //  float yshft=suncntrobs.v[1]*solar_r; //-ystart;

    float yshft=*fovpix*(pcntr[0]-float(*is-1)/2);
    float xshft=*fovpix*(pcntr[1]-float(*js-1)/2);

    // -- progression
    float progresspercent=0.2;
    int progressflag=(int) ((float) *nbang * progresspercent);
    float progresspass=progresspercent;

    ang=astart;
    for(j=0;j<(*nbang);j++,ang+=astep,posbtot++,posbpol++,posne++) {
        // -- print progression
        if (((*quiet) == 0) && (j > progressflag)) {
            cout << progresspass*100 << "% " << endl;
            progresspass+=progresspercent;
            progressflag=(int) ((float) *nbang * progresspass);
        }
        xx=xshft+radpix*cos(ang);
        yy=yshft+radpix*sin(ang);

        // -- compute the LOS direction and impact parameter
        getimpactandlos(xx,yy,obs,abs,vlosabs,rho);

    // ---- LOS integration
    losinteg(pofinteg,frontinteg,obs,vlosabs,abs,poi,los,rho,pmod,pmodparam,modpos,u,constfactor,flagneonly,posbtot,posbpol,posne,plosDepthIn,plosDepthOut,0);

    }
    delete pmod;
    return 1;
}



//! raytracing program that integrates only along a straight line segment profile
extern "C" int rtwlseg(int argc, void **argv) {

    int i,j,k;

    // ---- get the params from IDL
    int *is,*js;
    float *fovpix;
    float *obspos,*obsang;
    float *nepos,*neang;
    int *losnbp;
    float *losrange;
    int *pmodelid;
    float *btot;
    float *bpol;
    float *netot;
    float *pmodparam;
    float *crpix; //!< sun center position (returned)
    int *nbpix;
    float *segment;
    int *quiet;
    int flagneonly=0; // force brightness calculation for the moment
    float ppoiang[3]={0,0,0}; // POI orientation set to 0 for the moment
    float *Hlonlat;
    float *plimbdark;
    float *pobslonlat;
    int *pobslonlatflag;
    
    is=(int*) argv[0];
    js=(int*) argv[1];
    fovpix=(float*) argv[2];
    obspos=(float*) argv[3];
    obsang=(float*) argv[4];
    nepos=(float*) argv[5];
    neang=(float*) argv[6];
    losnbp=(int*) argv[7];
    losrange=(float*) argv[8];
    pmodelid=(int*) argv[9];
    btot=(float*) argv[10];
    bpol=(float*) argv[11];
    netot=(float*) argv[12];
    pmodparam=(float*) argv[13];
    crpix=(float*) argv[14];
    nbpix=(int*) argv[15];
    segment=(float*) argv[16];
    quiet=(int*) argv[17];
    Hlonlat=(float*) argv[18];
    plimbdark=(float*) argv[19];
    pobslonlat=(float*) argv[20];
    pobslonlatflag=(int*) argv[21];

    if (*quiet == 0) {
        cout << "Model : " << *pmodelid << endl;
        cout << "is : " << *is << " , js : " << *js << endl;
        cout << "fovpix : " << *fovpix << endl;
        cout << "obspos : " << *obspos;
        for(i=1;i<3;i++)
            cout << " , " << *(obspos + i);
        cout << endl << "obsang : " << *obsang;
        for(i=1;i<3;i++)
            cout << " , "  << *(obsang + i);
        cout << endl << "nepos : " << *nepos;
        for(i=1;i<3;i++)
            cout << " , "  << *(nepos + i);
        cout << endl << "neang : " << *neang;
        for(i=1;i<3;i++)
            cout << " , "  << *(neang + i);
        cout << endl << "losnbp : " << *losnbp << endl;
        cout << "losrange : " << *losrange;
        for(i=1;i<2;i++)
            cout << " , "  << *(losrange + i);
        cout << endl;
        cout << "Straight line profile integration" << endl;
        cout << "Nb pix along the line segment : " <<  *nbpix << endl;
        cout << "point 1 : " << segment[0] <<","<< segment[1] << endl;
        cout << "point 2 : " << segment[2] <<","<< segment[3] << endl;
        cout << "Limb Darkening : " << *plimbdark << endl;
        cout << "obslonlatflag : " << *pobslonlatflag << endl;
    }

    // ---- select the density model requested by the user
    CModelBase *pmod;
    pmod=modelselect(*pmodelid);
    pmod->initParam(pmodparam);

    // -- array position ptr
    float *posbtot,*posbpol,*posne;
    posbtot=btot;
    posbpol=bpol;
    posne=netot;

    float *plosDepthIn,*plosDepthOut;
    
    int pofinteg=0;
    int frontinteg=0;

    // -- Thomson scattering constants
    float u=*plimbdark; // limb darkening
    float constfactor=constfactor(u);

    float xx,yy;

    float rho,r,s;
    float ner,cosomega,sinomega;
    float cossquareomega,sinsquareomega;
    float rhooverr;
    // -- temperature: not useful in Thomson scattering but declared
    //    for compatibility of the models with the radio raytracing
    //float temperature;

    float a,b,c,d;
    float logterm;
    float polterm;//,tanterm;

    // ---- absolute basis definition
    Cbasis abs(Cvec(0,0,0),0,0,0);

    // ---- observer position definition
    Cbasis obs;
    if (*pobslonlatflag == 0) 
    obs=Cbasis(Cvec(obspos[0],obspos[1],obspos[2]),obsang[0],obsang[1],obsang[2]);   
    else 
    obs=Cbasis(pobslonlat[0],pobslonlat[1],pobslonlat[2],obsang[0],obsang[1],obsang[2]);

    // ---- Ne position
    Cbasis nps(Cvec(*nepos,*(nepos+1),*(nepos+2)),
            *neang, *(neang+1), *(neang+2),Hlonlat[0],Hlonlat[1],Hlonlat[2]);

    ModelPosition modpos;
    modpos.setBasis(nps);

    // ---- POI basis
    Cbasis poi(Cvec(0,0,0),ppoiang[0],ppoiang[1],ppoiang[2]);

    // -- LOS integration parameter definition
    Clos los(*losnbp,*losrange,*(losrange+1));

    // -- LOS vector in absolute basis
    Cvec vlosabs,vlosstep;

    // -- Sun cntr ortho proj on the LOS
    Cvec qlos;

    // -- current LOS point position
    Cvec vs;

    // -- compute Sun center position on the images
    Cvec suncntrobs;
    suncntrobs=orthoprojpointplane(Cvec(0,0,0),obs.u.row(3),-obs.o);
    Cvec vobs2sun=obs.u * (-obs.o);

    float phicntr=PI/2,thetacntr=0.;
    if (vobs2sun.v[2] != 0)
        thetacntr=atan(sqrt(sqr(vobs2sun.v[0])+sqr(vobs2sun.v[1]))/vobs2sun.v[2]);

    if (vobs2sun.v[1] !=0) {
        phicntr=atan(vobs2sun.v[0]/vobs2sun.v[1]);
        if (vobs2sun.v[1] < 0)
            phicntr+=PI;
    } else if (vobs2sun.v[0] < 0)
        phicntr=3*PI/2;

    crpix[0]=thetacntr*cos(phicntr)/(*fovpix)+float(*is-1)/2;
    crpix[1]=thetacntr*sin(phicntr)/(*fovpix)+float(*js-1)/2;

    float ystart=*fovpix*(segment[0]-float(*is-1)/2);
    float xstart=*fovpix*(segment[1]-float(*js-1)/2);

    float yend=*fovpix*(segment[2]-float(*is-1)/2);
    float xend=*fovpix*(segment[3]-float(*js-1)/2);

    float ystep=(yend-ystart)/(*nbpix-1);
    float xstep=(xend-xstart)/(*nbpix-1);
    if (*quiet == 0) {
        cout << "crpix[0] : " << crpix[0] << ", ";
        cout << "crpix[1] : " << crpix[1] << endl;
        cout << "xstart : " << xstart << endl;
        cout << "ystart : " << ystart << endl;
        cout << "xend : " << xend << endl;
        cout << "yend : " << yend << endl;
        cout << "xstep : " << xstep << endl;
        cout << "ystep : " << ystep << endl;
    }

    // -- progression
    float progresspercent=0.2;
    int progressflag=(int) ((float) *nbpix * progresspercent);
    float progresspass=progresspercent;

    xx=xstart;
    yy=ystart;

    for(j=0;j<(*nbpix);j++,xx+=ystep,yy+=xstep,posbtot++,posbpol++,posne++) {
        // -- print progression
        if (((*quiet) == 0) && (j > progressflag)) {
            cout << progresspass*100 << "% " << endl;
            progresspass+=progresspercent;
            progressflag=(int) ((float) *nbpix * progresspass);
        }

        // -- compute the LOS direction and impact parameter
        getimpactandlos(xx,yy,obs,abs,vlosabs,rho);

    // ---- LOS integration
    losinteg(pofinteg,frontinteg,obs,vlosabs,abs,poi,los,rho,pmod,pmodparam,modpos,u,constfactor,flagneonly,posbtot,posbpol,posne,plosDepthIn,plosDepthOut,0);

    }
    delete pmod;
    return 1;
}



//! Raytracing for a density cube
extern "C" int rtdenscube(int argc, void **argv) {

    // ---- get the params from IDL
    int *is,*js;
    float *fovpix;
    float *obspos,*obsang;
    float *nepos,*neang;
    float *btot;
    float *bpol;
    float *netot;
    float *pmodparam;
    float *crpix; //!<sun center position (returned)
    float *rho;
    float *mmlon;
    float *mmlat;
    float *rrr;
    int *ppofinteg;
    int *quiet;
    int *flagneonly;
    int *proi;
    float *ppoiang;
    float *Hlonlat;
    float *poccrad;
    float *plimbdark;
    int *pcubesize; //!< [x,y,z] size of the density cube edge in pix
    float *pvoxsize; //!< size of the voxel edge in Rsun
    float *dccenter; //!< [x,y,z] density cube center in pix
    int *pflagdenscubetype; //!< size of the density cube edge in pix
    unsigned char *pdensitylistfilename; //!< filename of the density list, if requested
    float *pobslonlat;
    int *pobslonlatflag;
    
    
    is=(int*) argv[0];
    js=(int*) argv[1];
    fovpix=(float*) argv[2];
    obspos=(float*) argv[3];
    obsang=(float*) argv[4];
    nepos=(float*) argv[5];
    neang=(float*) argv[6];
    btot=(float*) argv[7];
    bpol=(float*) argv[8];
    netot=(float*) argv[9];
    pmodparam=(float*) argv[10];
    crpix=(float*) argv[11];
    rho=(float*) argv[12];
    mmlon=(float*) argv[13];
    mmlat=(float*) argv[14];
    rrr=(float*) argv[15];
    ppofinteg=(int*) argv[16];
    quiet=(int*) argv[17];
    flagneonly=(int*) argv[18];
    proi=(int*) argv[19];
    ppoiang=(float*) argv[20];
    Hlonlat=(float*) argv[21];
    poccrad=(float*) argv[22];
    plimbdark=(float*) argv[23];
    pcubesize=(int*) argv[24];
    pvoxsize=(float*) argv[25];
    dccenter=(float*) argv[26];
    pflagdenscubetype=(int*) argv[27];
    pdensitylistfilename=(unsigned char*) argv[28];
    pobslonlat=(float*) argv[29];
    pobslonlatflag=(int*) argv[30];

    // ---- copy the filename in a string
    string lstfn;
    // -- I couldn't find a other way to process an string from IDL
    //    In IDL: bytarr('name.dat'). In C++ I get a array of char with a 0 in between each char ! Don't know why !
    unsigned int i=0;
    while (pdensitylistfilename[i] != '\0' && i <= 255) {
        lstfn+= pdensitylistfilename[i];
        i+=2; // the step is 2 to jump over the inserted 0
    }

    // ---- get the params from IDL
    if (*quiet == 0) {
        cout << "is : " << *is << " , js : " << *js << endl;
        cout << "fovpix : " << *fovpix << endl;
        cout << "obspos : " << *obspos;
        for(int i=1;i<3;i++)
            cout << " , " << *(obspos + i);
        cout << endl << "obsang : " << *obsang;
        for(int i=1;i<3;i++)
            cout << " , "  << *(obsang + i);
        cout << endl << "nepos : " << *nepos;
        for(int i=1;i<3;i++)
            cout << " , "  << *(nepos + i);
        cout << endl << "neang : " << *neang;
        for(int i=1;i<3;i++)
            cout << " , "  << *(neang + i);

        cout << endl;
        cout << "Hlonlat : " << Hlonlat[0] << " , " << Hlonlat[1] << " , " << Hlonlat[2] << endl;
        cout << "occrad : " << *poccrad << endl;
        cout << "limbdark : " << *plimbdark << endl;

        cout << "pcubesize : " << pcubesize[0] << ", "<< pcubesize[1]<< ", " <<pcubesize[2]<< endl;
        cout << "pvoxsize : " << *pvoxsize << endl;
        cout << "dccenter : " << dccenter[0] << ", " << dccenter[1] << ", " << dccenter[2] << endl;
        cout << "denscube : " << pmodparam[0] << endl;
        cout << "pflagdenscubetype : " << *pflagdenscubetype << endl;

        cout << "lstfn : " << lstfn << endl;
        cout << "obslonlatflag : " << *pobslonlatflag << endl;

    }

    // ---- test the formating of the input
    unsigned long posoffset=0;
    switch (*pflagdenscubetype) {
    case 1 :
        cout << "Compacted density cube detected" <<endl;
        posoffset=7;
        break;

    case 2 :
        cout << "Not implemented yet" << endl;
        return 1;
        break; // not really useful here :-)

    case 3 :
        cout << "Density list file" << endl;
        break;
    }

    // -- array position ptr
    float *posbtot,*posbpol,*posne,*posrhoim,*posrrr;
    posbtot=btot;
    posbpol=bpol;
    posne=netot;
    posrhoim=rho;
    posrrr=rrr;
    int *posroi;
    posroi=proi;

    int pofinteg=*ppofinteg;

    // -- Thomson scattering constants
    float u=*plimbdark; // -- limb darkening: in prevision that user can change it
    float constfactor=constfactor(u);

    float xstep=*fovpix;
    float ystep=*fovpix;
    float jcntr=(float)(*js-1)/2.;
    float icntr=(float)(*is-1)/2.;

    float xstart=-xstep*jcntr; // is and js are the size of the output image in the image usual coordinates
    float ystart=-ystep*icntr;

    // ---- absolute basis definition
    Cbasis abs(Cvec(0,0,0),0,0,0);

    // ---- observer position definition
    Cbasis obs;
    if (*pobslonlatflag == 0) 
    obs=Cbasis(Cvec(obspos[0],obspos[1],obspos[2]),obsang[0],obsang[1],obsang[2]);   
    else 
    obs=Cbasis(pobslonlat[0],pobslonlat[1],pobslonlat[2],obsang[0],obsang[1],obsang[2]);
    
    // ---- Ne position
    Cbasis nps(Cvec(*nepos,*(nepos+1),*(nepos+2)),
            *neang, *(neang+1), *(neang+2),Hlonlat[0],Hlonlat[1],Hlonlat[2]);

    // ---- POI basis
    Cbasis poi(Cvec(0,0,0),ppoiang[0],ppoiang[1],ppoiang[2]);


    // ---- force integration around POI if orientation is present
    if (ppoiang[0]!=0 || ppoiang[1]!=0 || ppoiang[2]!=0)
        pofinteg=1;

    // -- compute Sun center position on the image
    Cvec vobs2sun=obs.u * (-obs.o);
    if (*quiet == 0)
        cout << "vobs2sun : " << vobs2sun << endl;

    float phicntr=PI/2,thetacntr=0.;
    if (vobs2sun.v[2] != 0)
        thetacntr=atan(sqrt(sqr(vobs2sun.v[0])+sqr(vobs2sun.v[1]))/vobs2sun.v[2]);

    if (vobs2sun.v[1] !=0) {
        phicntr=atan(vobs2sun.v[0]/vobs2sun.v[1]);
        if (vobs2sun.v[1] < 0)
            phicntr+=PI;
    } else if (vobs2sun.v[0] < 0)
        phicntr=3*PI/2;

    crpix[0]=thetacntr*cos(phicntr)/ystep+float(*is-1)/2;
    crpix[1]=thetacntr*sin(phicntr)/xstep+float(*js-1)/2;


    // -- test if it's a density list
    ifstream binFile;
    if (*pflagdenscubetype == 3) {
        binFile.open(lstfn.c_str(),ios::in | ios::binary);
        if (!binFile.is_open()) {
            cout << "Cannot open "<< lstfn<<"!" << endl;
            return 0;
        }

        // -- check if density list file
        int densityfiletype;
        binFile.read(reinterpret_cast < char * > (&densityfiletype), sizeof(densityfiletype));
        cout << "densityfiletype : " << densityfiletype << endl;

        // -- overread pvoxsize
        binFile.read(reinterpret_cast < char * > (pvoxsize), sizeof(float));
        cout << "*pvoxsize : " << *pvoxsize << endl;


        if (densityfiletype != 3) {
            cout << "The density list file "<<  lstfn<< " has the wrong format." << endl;
            cout << "Cannot process." << endl;
            binFile.close();
            return 0;
        }


    }


    // ---- loop for each of the density cube voxel
    // -- init
    Cbasis vox=nps; // init voxel coord sys: only the origin will change in the loop

    unsigned long voxpos=posoffset;
    unsigned long countunder=0,countunder2=0,countbig=0;
    Cvec cubecntr(dccenter[0],dccenter[1],dccenter[2]); // density cube center
    float hvoxsize=*pvoxsize/2; // half of the voxel size to shift to voxel center

    Cvec *vlosimobs,*vlosimabs; // save the values of the LOS vectors to avoid recalculating them
    vlosimobs=new Cvec[*is * *js];
    vlosimabs=new Cvec[*is * *js];

    // -- init of the first voxel
    vox.o=nps.o+nps.ui*((*pvoxsize * (Cvec(0,0,0)-cubecntr) )-hvoxsize);
    Cdvoxel voxel(vox,*pvoxsize,pmodparam[voxpos],obs);

    // -- init projection pos of the vertices on the image
    std::vector<float> vi(voxel.vertobs.size()),vj(voxel.vertobs.size());
    std::vector<std::vector<float> > vij;
    vij.push_back(vi);
    vij.push_back(vj);

    // -- progression
    unsigned int nbvox;
    if (*pflagdenscubetype == 3) {
    // -- get the file size to calculate the number of voxel
    int filesize=FileSize(lstfn.c_str());
    cout << "file size : " << filesize << endl;
    nbvox=((((unsigned int) filesize)-sizeof(int)-sizeof(float))/sizeof(float))/4;
    cout << "Nb Voxels : " << nbvox << endl;
    } else {
    nbvox=pcubesize[0]*pcubesize[1]*pcubesize[2];
    }
    float progresspercent=0.2;
    unsigned int progressflag=(unsigned int) ((float) nbvox * progresspercent);
    float progresspass=progresspercent;
    
    
    // ---- check if it's a density list file
    if (*pflagdenscubetype == 3) {
        // ---- process in case of density loop file
    unsigned int voxcount=0;
        while (!binFile.eof()) {
            // ---- read position and density of voxel
            float xv,yv,zv,dv; // position and density of the record
            binFile.read(reinterpret_cast<char *>(&xv),sizeof(float));
            binFile.read(reinterpret_cast<char *>(&yv),sizeof(float));
            binFile.read(reinterpret_cast<char *>(&zv),sizeof(float));
            binFile.read(reinterpret_cast<char *>(&dv),sizeof(float));
            
            // -- print progression
            if (((*quiet) == 0) && (voxcount > progressflag)) {
            cout << progresspass*100 << "% " << endl;
            progresspass+=progresspercent;
            progressflag=(int) ((float) nbvox * progresspass);
            }
            voxcount++;
            
            // -- loop if the density within the voxel is quasi null
            if (dv < 0.1)
                continue;

            // -- position of the voxel origin vertice
            vox.o=nps.o+nps.ui*(Cvec(xv,yv,zv)-hvoxsize);

            voxel.UpdateOriginandDensity(vox.o,dv);
            
            // ---- proceed to voxel projection on image plane
#include "voxinteg.h"

        }

        binFile.close();

    } else {

        for (int k=0;k < pcubesize[2];k++)
            for (int j=0;j < pcubesize[1];j++)
                for (int i=0;i < pcubesize[0];i++,voxpos++) {
        
                // -- print progression
        if (((*quiet) == 0) && (voxpos > progressflag)) {
            cout << progresspass*100 << "% " << endl;
            progresspass+=progresspercent;
            progressflag=(int) ((float) nbvox * progresspass);
        }
        
                    // -- loop if the density within the voxel is quasi null
                    if (pmodparam[voxpos] < 0.1)
                        continue;

                    // -- position of the voxel origin vertice
                    vox.o=nps.o+nps.ui*((*pvoxsize * (Cvec(i,j,k)-cubecntr) )-hvoxsize);

                    voxel.UpdateOriginandDensity(vox.o,pmodparam[voxpos]);
                    // ---- proceed to voxel projection on image plane
#include "voxinteg.h"

                }
    }

    delete [] vlosimobs;
    delete [] vlosimabs;
    cout << "countunder : " << countunder << endl;
    cout << "countunder2 : " << countunder2 << endl;
    cout << "countbig : " << countbig << endl;

    return 1;
}


//! return the carrington position on the solar limb for a given pixel: useful for EUVI
extern "C" int rtGetCarPosOnDisk(int argc, void **argv) {
// ---- extract input parameters
int *is=(int*) argv[0]; // -- x size of the CCD
int *js=(int*) argv[1]; // -- y size of the CCD
float *fovpix=(float*) argv[2]; // -- resolution of the pix, in rad
float *obspos=(float*) argv[3]; // -- position of the observer
float *obsang=(float*) argv[4]; // -- angle position of the observer
float *crpix=(float*) argv[5]; // -- position of the optical axis, in pix
int *quiet=(int*) argv[6]; // -- flag for quiet mode
float *obslonlat=(float*) argv[7]; // -- car lon lat height of the obs [3]
int obslonlatflag=*((int*) argv[8]); // -- flag if obslonlat must be used
int projtypecode=*((int*) argv[9]); // -- projection type
float pv2_1=*((float*) argv[10]); // -- pv2_1 mu parameter for the azp projection
float *pc=(float*) argv[11]; // -- pc matrix [2,2]
float *ppixpos=(float*) argv[12]; // -- pixel position where to compute the CAR coord
float *pci=(float*) argv[13]; // -- inverse of pc matrix [2,2]
float *pcarlonlat=(float*) argv[14]; // -- (out) carrington position
int flagintersection=*((int*) argv[15]); // -- (out) 1 if intersection, 0 if no
int i,j,k;

// ---- init the different coordinate systems
// -- absolute basis definition
Cbasis abs(Cvec(0,0,0),0,0,0);

// -- observer position definition
Cbasis obs;
if (obslonlatflag == 0) 
    obs=Cbasis(Cvec(obspos[0],obspos[1],obspos[2]),obsang[0],obsang[1],obsang[2]);
else
    obs=Cbasis(obslonlat[0],obslonlat[1],obslonlat[2],obsang[0],obsang[1],obsang[2]);

// ---- compute the Vlos for the given pix
float rrr,alpha;
Cvec Vlosobs = xxyy2losobs(ppixpos[0],ppixpos[1],rrr,alpha,crpix,fovpix,pc,projtypecode,pv2_1);

// ---- compute intersection of the los with the solar sphere of 1Rsun
Cvec Vlosabs=obs.ui*Vlosobs;

// -- intersection of the LOS with the sphere
float root1,root2;
float det=get2ndDegPolyRoots(obs.o.magsqr()-1.,2*pscal(obs.o,Vlosabs),Vlosabs.magsqr(),root1,root2);

// -- no intersection if no roots: return
if (det < 0.) {
    flagintersection=0;
    return 1;
}

// -- use first root (smallest root) for foreground intersection
Cvec OP=root1*Vlosabs+obs.o;

// -- convert in spherical coordinates
// - shift the axis to put the orientation of the ABS coord sys in the Carrington one
float rho;
OP.leftshift();
cvcoord(OP,&rho,pcarlonlat,pcarlonlat+1);

pcarlonlat[0]+=PI/2.;
if (pcarlonlat[0] > TWOPI) pcarlonlat[0]-=TWOPI;

return 1;
}


//! Compute ii,jj pix pos for the voxel vertices
std::vector<std::vector<float> > voxelvert2iijj(std::vector<Cvec> &vert,const float &xstep,const float &ystep,const float &icntr,const float &jcntr) {
    std::vector<std::vector<float> > vij;
    std::vector<float> vi(vert.size()),vj(vert.size());
    float xx,yy;
    for(unsigned int i=0;i<vert.size();i++) {
        losobs2xxyy(vert[i],xx,yy);
        vj[i]=xx/xstep+jcntr;
        vi[i]=yy/ystep+icntr;
    }
    vij.push_back(vi);
    vij.push_back(vj);
    return vij;
}

//! Compute ii,jj pix pos for the voxel vertices
std::vector<float> voxelvert2iijj(Cvec &vert,const float &xstep,const float &ystep,const float &icntr,const float &jcntr) {
    std::vector<float> vij(2);
    float xx,yy;

    losobs2xxyy(vert,xx,yy);

    vij[0]=yy/ystep+icntr;
    vij[1]=xx/xstep+jcntr;

    return vij;
}


//! Convert LOS orientation in obs to pix coordinates
void losobs2xxyy(const Cvec &Vlosobs,float &xx,float &yy) {
    if (Vlosobs[0] == 0 && Vlosobs[1] == 0) {
        xx=0;
        yy=0;
        return;
    }

    float alpha=atan2(Vlosobs[0],Vlosobs[1]);
    float ca=cos(alpha),sa=sin(alpha);
    float rrr;
    if (fabs(ca) < fabs(sa)) {
        rrr=asin(Vlosobs[0]/Vlosobs.mag()/sa);
    } else {
        rrr=asin(Vlosobs[1]/Vlosobs.mag()/ca);
    }

    xx=rrr*sin(alpha);
    yy=rrr*cos(alpha);
}

void losobs2xxyydebug(const Cvec &Vlosobs,float &xx,float &yy) {

    if (Vlosobs[0] == 0 && Vlosobs[1] == 0) {
        cout << "In losobs2xxyydebug() : Vlosobs[0] == 0 && Vlosobs[1] == 0" << endl;
        xx=0;
        yy=0;
        return;
    }

    float alpha=atan2(Vlosobs[0],Vlosobs[1]);
    float ca=cos(alpha),sa=sin(alpha);
    float rrr;
    if (fabs(ca) < fabs(sa)) {
        rrr=asin(Vlosobs[0]/Vlosobs.mag()/sa);
    } else {
        rrr=asin(Vlosobs[1]/Vlosobs.mag()/ca);
    }

    float rold=acos(Vlosobs[2]/Vlosobs.mag());

    cout << "rrr  : "<< rrr<<endl;
    cout << "rold : "<< rold<<endl;
    cout << "Vlosobs : "<< Vlosobs<<endl;
    cout << "Vlosobs.mag() : "<< Vlosobs.mag()<<endl;
    cout << "ca : "<< ca<< ", sa : "<< sa<<endl;

    cout << "alpha : "<< alpha<<endl;

    xx=rrr*sin(alpha);
    yy=rrr*cos(alpha);

}


// ---- model call wrapper for Dr Marque who likes Fortran
extern "C" float modelwrapper(int& modelid,
                                float& x,float& y,float& z,
                                float& pmparam,float& temperature) {

    // ---- select the density model requested by the user
    CModelBase *pmod;
    pmod=modelselect(modelid);
        pmod->initParam(&pmparam);

    // ---- call density model
    float ner=pmod->Density(Cvec(x,y,z),temperature);

    return ner;
}


// ---- build a density cloud file for the frontend display (outputtype=0)
// ---- or an integer grid for IDL work (outputtype=1)
// ---- note that outputtype=2 builds a cube and works quietly
extern "C" int buildcloud(int argc, void **argv) {
    int modelid=*((int*) argv[0]);
    unsigned int cubesidenbpix=(unsigned int) *((int*) argv[1]);
    float cubesidersun=*((float*) argv[2]);
    int outputtype=*((int*) argv[3]);
    float* pmparam=(float*) argv[4];
    float* cubeoriginrsun=(float*) argv[5]; //!< [x0,y0,z0] position of first voxel
    int quiet=*((int*) argv[6]);

    bool noisy=(quiet == 0);

//     noisy=1;
//      if (outputtype == 2)
//          noisy=0;

    // ---- display parameters
    if (noisy) {
        cout << "modelid : " << modelid << endl;
        cout << "Cube Side pix : " << cubesidenbpix << endl;
        cout << "Cube Side Rsun : " << cubesidersun << endl;
        cout << "cubeoriginrsun : " <<cubeoriginrsun[0]<<" "<<cubeoriginrsun[1]<<" "<<cubeoriginrsun[2]<< endl;
    }

    // ---- select the density model requested by the user
    CModelBase *pmod;
    pmod=modelselect(modelid);
    pmod->initParam(pmparam);
    float nel=0;

    // ---- open output file
    ostringstream ostr;
    switch (outputtype) {
    case 0 :
        if (noisy)
            cout << "Creating dragger format" << endl;
        ostr << "output" << modelid << ".txt" ;
        break;

    case 1 :
        if (noisy)
            cout << "Creating IDL cube" << endl;
        ostr << "cube" << modelid << ".txt" ;
        break;

    case 2 :
        if (noisy)
            cout << "Creating binary cube" << endl;
        ostr << "cube" << modelid << ".dat" ;
        break;

    case 3 :
        if (noisy)
            cout << "Creating binary density list" << endl;
        ostr << "cube" << modelid << ".dat" ;
        break;

    default :
        if (noisy)
            cout << "Creating IDL cube, by default" << endl;
        ostr << "cube" << modelid << ".txt" ;
        outputtype=1;

    }

    // ---- open file in text or bin mode
    if (noisy)
        cout << "Writing file : " << ostr.str() << endl;
    ofstream myFile;
    fstream binFile;
    if (outputtype < 2) {
        myFile.open(ostr.str().c_str(),ios::out);
        if (!myFile.is_open()) {
            cout << "Cannot open file for writting !" << endl;
            return 0;
        }
    } else {
        binFile.open(ostr.str().c_str(),ios::out | ios::binary);
        if (!binFile.is_open()) {
            cout << "Cannot open file for writting !" << endl;
            return 0;
        }
    }

    // -- progression
    float ProgressPercent=0.2;
    float ProgressTotalCount=(float) (cubesidenbpix*cubesidenbpix*cubesidenbpix);
    unsigned int ProgressFlag=(unsigned int) (ProgressTotalCount * ProgressPercent);
    float ProgressPass=ProgressPercent;
    unsigned int ProgressCount=0;

    // ---- main loop that builds the cube
    float step=cubesidersun/((float) (cubesidenbpix-1));
    float x0=cubeoriginrsun[0],y0=cubeoriginrsun[1],z0=cubeoriginrsun[2];
    float x,y,z;
    unsigned int i,j,k;

    if (outputtype < 2) {
        // ---- text formated output
        // -- the first line of the file contain the modelid, the size
        //      and the resolution of the cube
        myFile << modelid << " " << cubesidenbpix << " " << cubesidersun << endl;

        // - scan space and write file
        for (k=0,z=z0;k<cubesidenbpix;k++,z+=step)
            for (j=0,y=y0;j<cubesidenbpix;j++,y+=step)
                for (i=0,x=x0;i<cubesidenbpix;i++,x+=step) {
                    if (Cvec(x,y,z).mag() > 2.0)
                        nel=pmod->Density(Cvec(x,y,z));
                    else
                        nel=0;

                    if (outputtype >= 1) {
                        myFile << i << " " << j << " " << k << " " << nel << endl;
                    } else {
                        if (nel > 0)
                            myFile << x << " " << y << " " << z << " " << nel << endl;
                    }

                    // -- print progression
                    if (noisy) {
                    if (ProgressCount > ProgressFlag) {
                            cout << ProgressPass*100 << "% " << endl;
                        ProgressPass+=ProgressPercent;
                        ProgressFlag=(int) (ProgressTotalCount * ProgressPass);
                    }
                    }
                    // -- inc progression counter
                    ProgressCount++;
                }
        // -- close density file
        myFile.close();

    } else {
        // ---- binary file output
        binFile.write(reinterpret_cast<char *>(&outputtype),sizeof(int)); // 1 : the density cube binary file type: if 2: density cube, if 3: density point list
        if (outputtype==2) {
            binFile.write(reinterpret_cast<char *>(&cubesidenbpix),sizeof(unsigned int));  // 2: the size of the cube side in pix (unsigned int)
            binFile.write(reinterpret_cast<char *>(&cubesidersun),sizeof(float));  // 3: the cube side size in Rsun (float)
            // 4,5,6 : position of the first voxel of the density cube
            binFile.write(reinterpret_cast<char *>(&x0),sizeof(float));
            binFile.write(reinterpret_cast<char *>(&y0),sizeof(float));
            binFile.write(reinterpret_cast<char *>(&z0),sizeof(float));
        } else {
            // ---- else it's a density voxel list
            // -- save the voxel side size
            binFile.write(reinterpret_cast<char *>(&step),sizeof(float));
        }
        for (k=0,z=z0;k<cubesidenbpix;k++,z+=step)
            for (j=0,y=y0;j<cubesidenbpix;j++,y+=step)
                for (i=0,x=x0;i<cubesidenbpix;i++,x+=step) {
                    if (Cvec(x,y,z).mag() > 1.1)
                        nel=pmod->Density(Cvec(x,y,z));
                    else
                        nel=0;
                    if (outputtype==2)
                        binFile.write(reinterpret_cast<char *>(&nel),sizeof(float));
                    else if (nel > 0.1 && outputtype==3) {
                        binFile.write(reinterpret_cast<char *>(&x),sizeof(float));
                        binFile.write(reinterpret_cast<char *>(&y),sizeof(float));
                        binFile.write(reinterpret_cast<char *>(&z),sizeof(float));
                        binFile.write(reinterpret_cast<char *>(&nel),sizeof(float));
                    }

                    // -- print progression
                    if (noisy) {
                    if (ProgressCount > ProgressFlag) {
                        cout << ProgressPass*100 << "% " << endl;
                        ProgressPass+=ProgressPercent;
                        ProgressFlag=(int) (ProgressTotalCount * ProgressPass);
                    }
                    }
                    // -- inc progression counter
                    ProgressCount++;
                }
        // -- close density file
        binFile.close();
    }
    return 1;
}


// ---- Get the electron density at a point [x,y,z] in space
extern "C" int getdensity(int argc, void **argv) {
    int* pmodelid=(int*) argv[0];
    float* point=(float*) argv[1];  // -- [x,y,z]
    float* pmparam=(float*) argv[2];
    float* neout=(float*) argv[3];

    // ---- select the density model requested by the user
    CModelBase *pmod;
    pmod=modelselect(*pmodelid);
    pmod->initParam(pmparam);

    // ---- init param and call the density
    Cvec v(point[0],point[1],point[2]);
        float tempe=0;
    float ner=pmod->Density(v);

    *neout=ner;

    return 1;

}

int getprojection(int &is,int &js,float &fovpix,
                float *obspos,float *obsang,
                float *nepos,float *neang,
                float *posin,float *posout,
                float *Hlonlat,int &nbpoints,
                bool flagclip) {

    float xstep=fovpix;
    float ystep=fovpix;
    float xstart=-xstep*(js-1)/2.; // is and js are the size of the output image in the image usual coordinates
    float ystart=-ystep*(is-1)/2.;

    // ---- absolute basis definition
    Cbasis abs(Cvec(0.,0.,0.),0,0,0);

    // ---- observer position definition
    Cbasis obs(Cvec(obspos[0],obspos[1],obspos[2]),
            obsang[0], obsang[1], obsang[2]);

    // ---- Ne position
    Cbasis nps(Cvec(nepos[0],nepos[1],nepos[2]),
            neang[0], neang[1], neang[2],
            Hlonlat[0],Hlonlat[1],Hlonlat[2]);

    Cvec Vlosobs,Vlosabs,OMne,OOobs=obs.o,Qlos;
    float rk,phi,theta,r,alpha;
    int pointpos;

    // -- comput projection for each points
    for (int i=0;i<nbpoints;i++) {
        pointpos=i*3;
        OMne=Cvec(posin[pointpos],
                posin[pointpos+1],
                posin[pointpos+2]);
        // compute LOS vector
        // in Absolute coord sys
        Vlosabs=nps.ui*OMne-OOobs;
        // -- check if behind the sun if requested
        if (flagclip) {
            // compute impact parameters
            float rho=psldist(obs.o,Vlosabs,abs.o);
            if (rho < 1.) {
                Qlos=orthoproj(obs.o,Vlosabs,abs.o);
                if (Vlosabs.norm() > ((Qlos-obs.o).norm()-sqrt(1.-Qlos.norm()*Qlos.norm()))) {
                    // It's behind: don't display
                    pointpos=i*2;
                    
#ifdef HAVE_SUNMATH_H
                    posout[pointpos]=quiet_nan(0);
                    posout[pointpos+1]=quiet_nan(0);
#else
                    posout[pointpos]=NAN;
                    posout[pointpos+1]=NAN;
#endif
                    continue;
                }
            }
        }

        // in Observer coord sys: image sys
        Vlosobs=obs.u*(Vlosabs/Vlosabs.norm());

        // X image axis : rotation around Y axis
        // Y image axis : rotation around X axis
        pointpos=i*2;
        posout[pointpos]=(Vlosobs.v[1]-ystart)/ystep;
        posout[pointpos+1]=(Vlosobs.v[0]-xstart)/xstep;

    }

    return 1;
}

extern "C" int getprojectionidl(int argc,void **argv) {

    int *is,*js;
    float *fovpix;
    float *obspos,*obsang;
    float *nepos,*neang;
    float *posin;
    float *posout;
    float *Hlonlat;
    int *nbpoints;
    int *pclip;

    is=(int*) argv[0];
    js=(int*) argv[1];
    fovpix=(float*) argv[2];
    obspos=(float*) argv[3];
    obsang=(float*) argv[4];
    nepos=(float*) argv[5];
    neang=(float*) argv[6];
    posin=(float*) argv[7];
    posout=(float*) argv[8];
    Hlonlat=(float*) argv[9];
    nbpoints=(int*) argv[10];
    pclip=(int*) argv[11];

    int flagok=getprojection(*is,*js,*fovpix,
                            obspos,obsang,
                            nepos,neang,
                            posin,posout,Hlonlat,
                            *nbpoints,(*pclip != 0));

    return flagok;
}



extern "C" int getdefaultparameteridl(int argc,void **argv) {
    int* pmodelid=(int*) argv[0];
    int* pnodefaultflag=(int*) argv[1];

    // ---- select the density model requested by the user
    CModelBase *pmod;
    pmod=modelselect(*pmodelid);

    // ---- get the default parameter description list
    std::vector<moddefparam> vp;
    int flagcase=1; //!< 1 : parameters undefined.
    pmod->dumpDefaultParamForIDL(vp,flagcase);
    *pnodefaultflag=flagcase;


    // ---- dump the structure in a text file (for the moment it's the easiest way to go)
    std::ofstream myFile("defmodparam.txt",ios::out);
    for (unsigned int i=0;i<vp.size();i++) {
        myFile <<   //i             << endl <<
        vp[i].keyword << endl <<
        vp[i].def     << endl <<
        vp[i].desc    << endl <<
        vp[i].units   << endl;
    }

    // -- Close file
    myFile.close();

    return 1;
}



int FileSize(const char* sFileName)
{
std::ifstream f;
f.open(sFileName, std::ios_base::binary | std::ios_base::in);
if (!f.good() || f.eof() || !f.is_open()) { return 0; }
f.seekg(0, std::ios_base::beg);
std::ifstream::pos_type begin_pos = f.tellg();
f.seekg(0, std::ios_base::end);
return static_cast<int>(f.tellg() - begin_pos);
}

void dumpBuildInfo()
{

#ifdef PACKAGE
    cout << "Package : " << PACKAGE <<endl;
#endif
#ifdef VERSION
    cout << "Version : " << VERSION <<endl;
#endif
    cout << "Compilation date : " << __DATE__ << " " << __TIME__ << endl; 
}


extern "C" int dumpbuildinfo(int argc,void **argv){
dumpBuildInfo();
return 1;
}



extern "C" int testParamPass(int param) {
//     int* pparam=(int*) argv[0];
    
    cout << "Param : " << param * 2 << " Hello ! " << endl;
 
    return EXIT_SUCCESS;;   
}


