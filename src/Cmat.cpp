//
// File: Cmat.cc
// $Id: Cmat.cpp,v 1.4 2009/03/06 21:22:46 thernis Exp $
//

#include "Cmat.h"

Cmat::Cmat()
{
	// constructor code here
	for(int i=0;i<3;i++) {
		for(int j=0;j<3;j++) {
		m[i][j]=0;
		}
	}
}

Cmat::Cmat(const Cmat &a) {
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      m[i][j]=a.m[i][j];
    }
  }
}

Cmat::Cmat(const float a[3][3]) {
  int i,j;
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      m[i][j]=a[i][j];
    }
  }
}

Cmat::Cmat(const float a) {
  int i,j;
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      m[i][j]=a;
    }
  }
}
  
Cmat::Cmat(float a00,float a10,float a20,float a01,float a11,float a21,float a02,float a12,float a22) {
  m[0][0]=a00;m[1][0]=a10;m[2][0]=a20;
  m[0][1]=a01;m[1][1]=a11;m[2][1]=a21;
  m[0][2]=a02;m[1][2]=a12;m[2][2]=a22;
}

// -- column
Cvec Cmat::column(int col) {
  Cvec r;
  switch (col) {
  case 1 : r.v[0]=m[0][0];r.v[1]=m[0][1];r.v[2]=m[0][2];
    break;
  case 2 : r.v[0]=m[1][0];r.v[1]=m[1][1];r.v[2]=m[1][2];
    break;
  case 3 : r.v[0]=m[2][0];r.v[1]=m[2][1];r.v[2]=m[2][2];
    break;
  default: r.v[0]=m[0][0];r.v[1]=m[0][1];r.v[2]=m[0][2];
    break;
  }
  return r;
}

// -- row
Cvec Cmat::row(int row) {
  Cvec r;
  switch (row) {
  case 1 : r.v[0]=m[0][0];r.v[1]=m[1][0];r.v[2]=m[2][0];
    break;
  case 2 : r.v[0]=m[0][1];r.v[1]=m[1][1];r.v[2]=m[2][1];
    break;
  case 3 : r.v[0]=m[0][2];r.v[1]=m[1][2];r.v[2]=m[2][2];
    break;
  default: r.v[0]=m[0][0];r.v[1]=m[1][0];r.v[2]=m[2][0];
    break;
  }
  return r;
}


// --- op overload
std::ostream& operator << (std::ostream& os,const Cmat& m) {
  os << std::endl 
    << m.m[0][0] << " , " << m.m[1][0] << " , " << m.m[2][0] << std::endl
     << m.m[0][1] << " , " << m.m[1][1] << " , " << m.m[2][1] << std::endl 
     << m.m[0][2] << " , " << m.m[1][2] << " , " << m.m[2][2];

  return os;
}


Cmat& Cmat::operator = (const Cmat &a) {
  int i,j;
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      m[i][j]=a.m[i][j];
    }
  }
  return *this;
}


Cmat Cmat::operator * (Cmat b) {
  Cmat a;
  
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      a.m[i][j]=m[0][j]*b.m[i][0]+
	m[1][j]*b.m[i][1]+
	m[2][j]*b.m[i][2];
    }
  }
  return(a); 
}

Cmat operator * (const Cmat &a,const Cmat &b){
  Cmat c;
  
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      c.m[i][j]=a.m[0][j]*b.m[i][0]+
          a.m[1][j]*b.m[i][1]+
          a.m[2][j]*b.m[i][2];
    }
  }
  return(c); 
  
  
}




Cvec Cmat::operator [](int col)
{
	return column(col+1);
}


bool Cmat::operator ==(const Cmat &a) const {
  return (m[0][0]==a.m[0][0] && m[0][1]==a.m[0][1] && m[0][2]==a.m[0][2] && 
					m[1][0]==a.m[1][0] && m[1][1]==a.m[1][1] && m[1][2]==a.m[1][2] && 
					m[2][0]==a.m[2][0] && m[2][1]==a.m[2][1] && m[2][2]==a.m[2][2]);
}


Cvec operator * (const Cmat &m,const Cvec &v)
 {
  Cvec r;
  r.v[0]=m.m[0][0]*v.v[0]+m.m[1][0]*v.v[1]+m.m[2][0]*v.v[2];
  r.v[1]=m.m[0][1]*v.v[0]+m.m[1][1]*v.v[1]+m.m[2][1]*v.v[2];
  r.v[2]=m.m[0][2]*v.v[0]+m.m[1][2]*v.v[1]+m.m[2][2]*v.v[2];
  return r;

}

Cmat& Cmat::rotmat(float a,int ax) {
  switch (ax) {
  case 1 : 
    m[0][0]=1.;m[1][0]=0.;m[2][0]=0.;
    m[0][1]=0.;m[1][1]=cos(a);m[2][1]=sin(a);
    m[0][2]=0.;m[1][2]=-sin(a);m[2][2]=cos(a);
    break;
  case 2 : 
    m[0][0]=cos(a);m[1][0]=0;m[2][0]=-sin(a);
    m[0][1]=0;m[1][1]=1;m[2][1]=0;
    m[0][2]=sin(a);m[1][2]=0;m[2][2]=cos(a);
    break;
  case 3 :
    m[0][0]=cos(a);m[1][0]=sin(a);m[2][0]=0;
    m[0][1]=-sin(a);m[1][1]=cos(a);m[2][1]=0;
    m[0][2]=0;m[1][2]=0;m[2][2]=1;
    break;
  }
  return *this;
}


Cmat::Cmat(const float a,const int ax) {
  rotmat(a,ax);
}

float Cmat::determinant()
{
	return m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2])
			-m[1][0]*(m[0][1]*m[2][2]-m[2][1]*m[0][2])
			+m[2][0]*(m[0][1]*m[1][2]-m[1][1]*m[0][2]);	
}

float Cmat::subdet(int i,int j)
{
	int i1=(i+1) % 3;
	int i2=(i+2) % 3;
	int j1=(j+1) % 3;
	int j2=(j+2) % 3;
	
	return m[i1][j1]*m[i2][j2]-m[i2][j1]*m[i1][j2];
}

Cmat Cmat::tranpose()
{
	Cmat t;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++) 
			t.m[i][j]=m[j][i];
	
	return t;		
}

Cmat const Cmat::inverse()
{
    // use A^-1=adj(A)/det(A)
	Cmat inv;
	float det=determinant();
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			inv.m[i][j]=subdet(j,i)/det;
	return inv;	
}



/*
* $Log: Cmat.cpp,v $
* Revision 1.4  2009/03/06 21:22:46  thernis
* Overload == operator
*
* Revision 1.3  2009/02/09 20:51:05  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.2  2007/05/14 17:19:40  thernis
* Add CVS id and log in all files
*
*/
