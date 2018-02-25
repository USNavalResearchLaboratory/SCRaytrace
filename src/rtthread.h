/*
 *  rtthread.h
 *  $Id: rtthread.h,v 1.3 2009/04/13 21:01:32 thernis Exp $
 *
 */

#ifndef RTTHREAD_H
#define RTTHREAD_H

/** \file rtthread.h
 * \brief Raytrace with threads to speed up on multi-core.
 */
 
int rtthread(int imsize0,int imsize1,float fovpix,float *obspos,float *obsang,float *nepos,float *neang,int losnbp,float *losrange,int modelid,float *btot,float *bpol,float *netot,float *pmodparam,float *crpix,int quiet,int neonly,float *hlonlat,float occrad,float limbdark,float *obslonlat,int obslonlatflag,int projtypecode,float pv2_1,float *pc,int frontinteg,unsigned int nbthreads,unsigned int nbchunk,float *nepos2,float *neang2);
int rtthreadtest();
//! Raytracing with threads to speed up on multi-core computers
extern "C" int rtthreadidl(int argc, void **argv);



#endif
