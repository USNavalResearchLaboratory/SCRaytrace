

bool PhysicsMie::computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btot,float &bpol,float &density)
{

    density = pparentscene->pmod->Density(ChangetoDensityCoord(pparentscene->modelposition,vs));

//     if (neout <= 1e-1) return 1;
//    if (pparentscene->neonly) return 0;


    // -- compute distance observer - point of LOS
    Cvec vx = vs - pparentscene->obs.o;
    float x = vx.mag();

/*    printvar(x);
    printvar(vs);
    printvar(pparentscene->obs.o);*/
    
    
    // -- Scattering phase function is constant for now
    float distObs2Sun = (pparentscene->obs.o).mag();
//     float sinTheta = rho / distObs2Sun;
//     float cosTheta = pscal(pparentscene->obs.o,vx) / distObs2Sun / x;
//    float theta = atan2(sinTheta, cosTheta); //asin(sinTheta);
//    float cosTheta = cos(theta);
//     float phi = atan2(cosTheta - x / distObs2Sun, sinTheta);
    
    float cosPhi = 0.5 * (distObs2Sun * distObs2Sun - x * x - r * r) / (x * r);
    float phi = acos(cosPhi);
/*    printvar(rho);
    printvar(theta);
    printvar(phi);*/
    
    
    
    
    // -- compute phase function here (function of phi)
//     mie->calcS(abs(phi));
//     s1 = mie->getS1();
//     s2 = mie->getS2();

    unsigned int idx;
    idx = (unsigned int) (float(SIZELOOKUPTABLE-1) * phi / PI);
    
    if (idx >= SIZELOOKUPTABLE) idx = SIZELOOKUPTABLE-1;
    
    float integrand = density / ( r * r );

    // -- assume Rsun = 1, and Lsun = 1
    
    btot = (PI / (1. + r * r)) * (i1table[idx] + i2table[idx]) * (0.5 / (k * k)) * density;     
//     btot = (PI / (1. + r * r)) * (i1table[idx] + i2table[idx]) * (0.5 / (k * k)) * density; 
    //bpol = PI / (1. + r * r) * (i2table[idx] - i1table[idx]) * 0.5 / (k * k) * density; 

    bpol = density / ( r * r ) * (i2table[idx] + i1table[idx]) * 0.5 / (k * k);
    
    
    
    // ---- Formulas from Van de Hulst, see also Stokes parameters I Q U V
    /*
    btot = integrand * (i1table[idx] + i2table[idx]);      // I
    bpol = integrand * (i2table[idx] - i1table[idx]);      // Q

    density = phi;
   */
    
    
    
    
    return 0;
}







bool PhysicsMie::computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btot,float &bpol,float &density)
{

    density = pparentscene->pmod->Density(ChangetoDensityCoord(pparentscene->modelposition,vs));


    // -- compute distance observer - point of LOS
    Cvec vx = vs - pparentscene->obs.o;
    float x = vx.mag();

    
    
    // -- Scattering phase function is constant for now
    float distObs2Sun = (pparentscene->obs.o).mag();
    
    float cosPhi = (distObs2Sun * distObs2Sun - x * x - r * r) / (2. * x * r);
    float phi = acos(cosPhi);

    unsigned int idPhi;
    idPhi = (unsigned int) (float(SIZELOOKUPTABLE-1) * phi / PI);
    
    if (idPhi >= SIZELOOKUPTABLE) idPhi = SIZELOOKUPTABLE-1;
    

    
    btot = 0.5 * (i1table[idPhi] + i2table[idPhi])  / (r * r * r);
    bpol = 0.5 * (i2table[idPhi] - i1table[idPhi])  / (r * r * r);

    
    
    //     btot = (0.5 * (i1table[idPhi] + i2table[idPhi]) / (k * k))  / pow(r, 3.);
//     bpol = (0.5 * (i2table[idPhi] - i1table[idPhi]) / (k * k))  / pow(r, 3.);

    
    
//     btot = (PI / (1. + r * r)) * (i1table[idPhi] + i2table[idPhi]) * (0.5 / (k * k)) * density;     
//     btot = (PI / (1. + r * r)) * (i1table[idPhi] + i2table[idPhi]) * (0.5 / (k * k)) * density; 
    //bpol = PI / (1. + r * r) * (i2table[idPhi] - i1table[idPhi]) * 0.5 / (k * k) * density; 
//     bpol = density / ( r * r ) * (i2table[idPhi] + i1table[idPhi]) * 0.5 / (k * k);
    
    
    return 0;
}



