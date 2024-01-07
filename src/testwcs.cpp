/*! \file testwcs.cpp 
 * \brief Test functionality of wcslib.
 *
 */

#include <string.h>
#include <iostream>
#include <wcslib/wcslib.h>

#include "camera.h"

using namespace std;


int main(int argc,char **argv) {
    cout << "In testwcs." << endl;

    struct prjprm prj;
    const double tol = 1.0e-9;
    int  status;

    prjini(&prj);

    // AZP: zenithal/azimuthal perspective.
    static const char pcode[] = "ARC";
    // prj.pv[1] =   0.5;
    // prj.pv[2] =  30.0;

    strcpy(prj.code, pcode);
    prjset(&prj);

    double x, y;
    x = 100; y = 100;
    double phi, theta;
    int stat;

    if (prj.prjx2s(&prj, 1, 0, 1, 1, &x, &y, &phi, &theta, &stat) == 1) {
        printf("   %3s(X2S) ERROR 1: %s\n", pcode, prj_errmsg[1]);
    }

    cout << "x     : " << x << endl;
    cout << "y     : " << y << endl;
    cout << "phi   : " <<   phi << endl;
    cout << "theta : " << theta << endl;
    cout << "stat  : " << stat << endl;
    // prjprt(&prj);


    // ---- Now test the implementation in the camera.h
    // cout << endl << "Checking Camera class..." << endl;
    // Detector detA;
    // Camera camA(0.1, ARC, detA);

    // camA.printWCS();

    // ---- INitialize a WCS structure
    wcsprm* wcs = new wcsprm;
    wcs->flag = -1;
    const int NAXIS = 2;

    wcsini(1, NAXIS, wcs);
    float crpix[2] = {991.547, 1015.11};
    wcs->crpix[0] = 991.547;
    wcs->crpix[1] = 1015.11;
    double *pcij;
    pcij = wcs->pc;
    // *(pcij++) = 0.975380644444; // [0, 0] = [i, j]
    // *(pcij++) = 0.202545920733; // [0, 1]
    // *(pcij++) = -0.202348686778; // [1, 0]
    // *(pcij++) = 0.975099606459; // [1, 1]
    *(pcij++) = 1; // [0, 0] = [i, j]
    *(pcij++) = 0; // [0, 1]
    *(pcij++) = 0; // [1, 0]
    *(pcij++) = 1; // [1, 1]
    float fovpix_deg = 0.0211525;
    wcs->cdelt[0] = fovpix_deg;
    wcs->cdelt[1] = fovpix_deg;
    string CTYPE[2] = {"HPLN-ZPN", "HPLT-ZPN"};
    // string CTYPE[2] = {"HPLN-ARC", "HPLT-ARC"};
    strcpy(wcs->ctype[0], CTYPE[0].c_str());
    strcpy(wcs->ctype[1], CTYPE[1].c_str());
    wcs->crval[0] = 31.89;
    wcs->crval[1] = -7.57;
    wcs->lonpole = 180;
    wcs->latpole = 0;

    #define NUMB_PV 9
    int NPV = NUMB_PV;
    struct pvcard PV[NUMB_PV];

    // PV values based on WISPR inner header
    PV[0].i = 1;			// Longitude is on axis 1.
    PV[0].m = 1;			// Parameter number 1.
    PV[0].value =  0.0;		// Fiducial native longitude.

    PV[1].i = 1;			// Longitude is on axis 1.
    PV[1].m = 2;			// Parameter number 2.
    PV[1].value = 90.0;		// Fiducial native latitude.

    PV[2].i = 1;			// Longitude is on axis 1.
    PV[2].m = 3;			// Parameter number 3.
    PV[2].value = 180.0;	// 

    // Set the PVi_m key values for ZPN projection.
    PV[3].i = 2;			// Latitude is on axis 2.
    PV[3].m = 0;			// Parameter number 0.
    PV[3].value = -1.98789997796E-08;

    PV[4].i = 2;			// Latitude is on axis 2.
    PV[4].m = 1;			// Parameter number 1.
    PV[4].value = 1.00501000881;

    PV[5].i = 2;			// Latitude is on axis 2.
    PV[5].m = 2;			// Parameter number 2.
    PV[5].value = 0.0729582980275;

    PV[6].i = 2;			// Latitude is on axis 2.
    PV[6].m = 3;			// Parameter number 3.
    PV[6].value = 0.275292009115;

    PV[7].i = 2;			// Latitude is on axis 2.
    PV[7].m = 4;			// Parameter number 4.
    PV[7].value = -0.701880991459;

    PV[8].i = 2;			// Latitude is on axis 2.
    PV[8].m = 5;			// Parameter number 5.
    PV[8].value = 1.97518002987;

    wcs->npv = NPV;
    int i, j;
    for (i = 0; i < NPV; i++) {
        wcs->pv[i] = PV[i];
    }

    if (wcsset(wcs)) {
        wcsperr(wcs, "Error");
        std::cout << "WCS no good" << std::endl;
    } else {
        std::cout << "All good" << std::endl;
    }


    // cout << endl << "Print WCS structure:" << endl;
    // wcsprt(wcs);
    # define NCOORD 6
    double pixel[NCOORD][2] = {{1024, 1024}, 
                          {991.547, 1015.11},
                          {991.547, 1015.11 + 300},
                          {991.547, 1015.11 - 300},
                          {991.547 + 300, 1015.11},
                          {991.547 - 300, 1015.11}};
    double imgcrd[NCOORD][2];
    double phi0[NCOORD];
    double theta0[NCOORD];
    double world[NCOORD][2];
    int stat0[NCOORD];
    if (wcsp2s(wcs, NCOORD, 2, pixel[0], imgcrd[0], 
                phi0, theta0, world[0], stat0)) {
        wcsperr(wcs, "  ");
        cout << "Hoy!" << endl;
    }

    cout << "stat : " << stat0[0] << endl;

    cout << "Pixel: " << endl;
    for (i = 0; i < NCOORD; i++)
    {
        for(j = 0; j < 2; j++)
        {
            cout << pixel[i][j] << ", ";
        }
        cout << endl;
    }

    cout << "theta: " << endl;
    for (i = 0; i < NCOORD; i++)
    {
        cout << theta0[i] << ", ";
    }
    cout << endl;

    cout << "phi: " << endl;
    for (i = 0; i < NCOORD; i++)
    {
        cout << phi0[i] << ", ";
    }
    cout << endl;


    cout << "world: " << endl;
    for (i = 0; i < NCOORD; i++)
    {
        for(j = 0; j < 2; j++)
        {
            cout << world[i][j] << ", ";
        }
        cout << endl;
    }

    cout << "imgcrd: " << endl;
    for (i = 0; i < NCOORD; i++)
    {
        for(j = 0; j < 2; j++)
        {
            cout << imgcrd[i][j] << ", ";
        }
        cout << endl;
    }




    cout << endl << "Checking Camera class using full wcs..." << endl;
    int pv_i[NUMB_PV] = {1,1,1,2,2,2,2,2,2};
    int pv_m[NUMB_PV] = {1,2,3,0,1,2,3,4,5};
    float pv[NUMB_PV] = {0.0, 90., 180.,
                        -1.98789997796E-08,
                         1.00501000881,
                         0.0729582980275,
                         0.275292009115,
                         -0.701880991459,
                         1.97518002987};
    Detector detB(2048, 1920, 20.48, 19.20);
    Camera camB(fovpix_deg / RADEG,
                ZPN,
                detB,
                crpix[0], crpix[1],
                NPV,
                pv, pv_i, pv_m);

    Cvec v_wcs, v_rt;
    v_wcs = camB.ij2loswcs(1024, 1024);
    v_rt = camB.ij2los(1024, 1024);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;

    v_wcs = camB.ij2loswcs(991.547, 1015.11);
    v_rt = camB.ij2los(991.547, 1015.11);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;

    v_wcs = camB.ij2loswcs(991.547, 1015.11 + 300);
    v_rt = camB.ij2los(991.547, 1015.11 + 300);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;

    v_wcs = camB.ij2loswcs(991.547, 1015.11 - 300);
    v_rt = camB.ij2los(991.547, 1015.11 - 300);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;

    v_wcs = camB.ij2loswcs(991.547 + 300, 1015.11);
    v_rt = camB.ij2los(991.547 + 300, 1015.11);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;

    v_wcs = camB.ij2loswcs(991.547 - 300, 1015.11);
    v_rt = camB.ij2los(991.547 - 300, 1015.11);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;


    cout << endl << "Checking Camera class using full wcs, ARC case..." << endl;
    int pv_iC[3] = {1,1,1};
    int pv_mC[3] = {1,2,3};
    float pvC[3] = {0.0, 90., 180.};
    Detector detC(2048, 1920, 20.48, 19.20);
    Camera camC(fovpix_deg / RADEG,
                ARC,
                detC,
                crpix[0], crpix[1],
                3,
                pvC, pv_iC, pv_mC);
    v_wcs = camC.ij2loswcs(1024, 1024);
    v_rt = camC.ij2los(1024, 1024);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;

    v_wcs = camC.ij2loswcs(991.547, 1015.11);
    v_rt = camC.ij2los(991.547, 1015.11);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;

    v_wcs = camC.ij2loswcs(991.547, 1015.11 + 300);
    v_rt = camC.ij2los(991.547, 1015.11 + 300);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;

    v_wcs = camC.ij2loswcs(991.547, 1015.11 - 300);
    v_rt = camC.ij2los(991.547, 1015.11 - 300);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;

    v_wcs = camC.ij2loswcs(991.547 + 300, 1015.11);
    v_rt = camC.ij2los(991.547 + 300, 1015.11);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;


    v_wcs = camC.ij2loswcs(991.547 - 300, 1015.11);
    v_rt = camC.ij2los(991.547 - 300, 1015.11);
    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt  : " << v_rt << endl;


    cout << endl << "Checking Camera" << endl;
    Camera camD;

    // -- Test AZP projection
    cout << endl << "Checking AZP" << endl;
    Camera camF(fovpix_deg / RADEG,
                    AZP,
                    detC,
                    crpix[0], crpix[1],
                    0.3);  // pv2_1

    v_wcs = camF.ij2loswcs(1024, 1024);
    v_rt = camF.ij2losold(1024, 1024);

    cout << "v_wcs : " << v_wcs << endl;
    cout << "v_rt : " << v_rt << endl;
    cout << "crpix[0] : " << crpix[0] << endl;
    cout << "crpix[1] : " << crpix[1] << endl;
    cout << "fovpix [rad] : " << fovpix_deg / RADEG << endl;





    delete wcs;

    cout << "Done." << endl;

}


