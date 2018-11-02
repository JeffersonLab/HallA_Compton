#define VJones_cxx

#include "VJones.h"


VJones::VJones(TComplex x,TComplex y)
{  
  PI= acos(-1);
  deg= 180./PI;
  rad =PI/180.;

  fv1 = TComplex(x.Re(),x.Im());
  fv2 = TComplex(y.Re(),y.Im());  

}

VJones::VJones()
{  
  PI= acos(-1);
  deg= 180./PI;
  rad =PI/180.;

  fv1 = TComplex(0,0);
  fv2 = TComplex(0,0);  
}

VJones::~VJones()
{  
}

//---------------------------------------------------
// Initialization of the two complex elements
//---------------------------------------------------
void VJones::Init(TComplex x,TComplex y)
{
  fv1 = TComplex(x.Re(),x.Im());
  fv2 = TComplex(y.Re(),y.Im());
  return;
}

//---------------------------------------------------
// Initialization of the two complex elements
// with constraint of normalization to 1
//---------------------------------------------------
void VJones::InitNorm(TComplex x,TComplex y)
{
  Double_t Norm=1./sqrt(x.Rho2()+y.Rho2());   
  fv1 = x*Norm;
  fv2 = y*Norm;
  return;
}

//---------------------------------------------------
// Return complex number, either 0 or 1  (x or y)
//---------------------------------------------------
TComplex VJones::GetCmplx(Int_t iv)
{
  TComplex v1 = TComplex(0,0);
  if (iv==0) v1 = fv1;
  if (iv==1) v1 = fv2;
  return v1;
}

//---------------------------------------------------
// Return modulus of Jones vector
//---------------------------------------------------
Double_t VJones::Norm()
{
  return fv1.Rho2()+fv2.Rho2();
  //  return sqrt(fv1.Rho2()+fv2.Rho2());
}

//---------------------------------------------------
// Return phase of Jones vector (in rad)
// (relative phase Phi_y-Phi_x)
//---------------------------------------------------
Double_t VJones::Phase()
{
  return fv2.Theta()-fv1.Theta();
}

//---------------------------------------------------
// Apply to following phase shift (in rad): 
// Phi_x --> Phi_x + Delta/2
// Phi_y --> Phi_y - Delta/2
//---------------------------------------------------
void VJones::PhaseShift(Double_t Delta)
{
  fv1 *= TComplex(cos(Delta/2.),sin(Delta/2.));
  fv2 *= TComplex(cos(-Delta/2.),sin(-Delta/2.));
  return;
}

//---------------------------------------------------
// Apply a global phase shift (Delta in rad): 
// No physical effect except when dealing with inteference issues
//---------------------------------------------------
void VJones::GlobalPhase(Double_t Delta)
{
  fv1 *= TComplex(cos(Delta),sin(Delta));
  fv2 *= TComplex(cos(Delta),sin(Delta));
  return;
}

//---------------------------------------------------
// Expression of a Jones vector in a coordinate frame 
// rotated by an angle Alpha (in rad)
//---------------------------------------------------
void VJones::Rotate(Double_t Alpha)
{
  TComplex tv1 = fv1*cos(Alpha) + fv2*sin(Alpha);
  TComplex tv2 = fv1*sin(-Alpha) + fv2*cos(Alpha);
  fv1=tv1;
  fv2=tv2;
  return;
}

//---------------------------------------------------
// Return angle Alpha (in rad) of the polarization ellipse 
// (modulo 180 degrees)
//---------------------------------------------------
Double_t VJones::Alpha()
{
  Double_t Alpha,tan2Alpha,Denom; 

  Denom = fv1.Rho2()-fv2.Rho2();
  if (Denom) 
  {
    tan2Alpha = 2.*fv1.Rho()*fv2.Rho()*cos(fv2.Theta()-fv1.Theta())/Denom;
    Alpha=atan(tan2Alpha)/2.;
  }else
    Alpha=PI/4.;

  VJones E = VJones(*this);
  //Rotate the ellipse by Alpha so that one axis of the ellipse is along X
  E.Rotate(Alpha);
  //If the ellipse axis along X is the small one, rotate by 90 degrees 
  //so that Alpha is always the angle of the large axis of the ellipse
  if (E.fv1.Rho()<E.fv2.Rho()) Alpha+=PI/2.;
  
  //Keep axis direction in [0,PI]
  if (Alpha>PI) Alpha-=PI; 
  if (Alpha<0 && Alpha>=-PI) Alpha+=PI; 
  if (Alpha<-PI) Alpha+=2*PI; 
    
  return Alpha;
}

//---------------------------------------------------
// Projection of a Jones vector on the basis of the 
// circular Left/Right states:
//     1    | 1  i | |v1   |Projection on right state
//  ------- |      | |     = | 
//  sqrt(2) | 1 -i | |v2   |Projection on left state
//---------------------------------------------------
void VJones::ProjCirc()
{  
  TComplex tv1 = fv1*(1./sqrt(2)) + fv2*TComplex(0.,1./sqrt(2));
  TComplex tv2 = fv1*(1./sqrt(2)) + fv2*TComplex(0.,-1./sqrt(2));
  fv1 = tv1;
  fv2 = tv2;
  return;
}

//---------------------------------------------------
// Return the Degree Of Circular Polarization of Jones vector
// assuming the light is fully polarized. 
// Convention: DOCP(Left)=+1, DOCP(Right)=-1 
//---------------------------------------------------
Double_t VJones::DOCP()
{
  VJones Proj = VJones(*this);
  Proj.ProjCirc();
  Double_t retval = (Proj.fv2.Rho2()-Proj.fv1.Rho2())/(Proj.fv2.Rho2()+Proj.fv1.Rho2());
  return retval;
}

Double_t VJones::DOLP()
{
  VJones Proj = VJones(*this);
  Double_t retval = sqrt(1-Proj.DOCP()*Proj.DOCP());
  return retval;
}


//---------------------------------------------------
// Print out of a Jones vector
//---------------------------------------------------
void VJones::Print()
{  
//   cout << "v1: (Re,Im)=(" << fv1.Re()  << "," << fv1.Im()<<")\t";
//   cout << "(Rho,Theta)=(" << fv1.Rho() << "," << fv1.Theta()*deg<<")"<<endl;
//   cout << "v2: (Re,Im)=(" << fv2.Re()  << "," << fv2.Im()<<")\t";
//   cout << "(Rho,Theta)=(" << fv2.Rho() << "," << fv2.Theta()*deg<<")"<<endl;
//   cout<<endl;
//   cout<<"Phase Phi: \t"<< Phase()*deg <<endl;
//   cout<<"Angle Alpha: \t"<< Alpha()*deg <<endl;
//   cout<<"DOCP: \t"<< DOCP() <<endl;
//   cout<<endl;
  printf("v1: (Re,Im)=(%+.4f,%+.4f)\t(Rho,Theta)=(%+.4f,%+.4f)\n",fv1.Re(),fv1.Im(),fv1.Rho(),fv1.Theta()*deg);
  printf("v2: (Re,Im)=(%+.4f,%+.4f)\t(Rho,Theta)=(%+.4f,%+.4f)\n",fv2.Re(),fv2.Im(),fv2.Rho(),fv2.Theta()*deg);
  printf("Phase Phi  = %+.2f       \tAngle Alpha= %+.2f\n", Phase()*deg, Alpha()*deg);
  printf("DOCP       = %+.4f       \tDOLP       = %+.4f\n", DOCP(), DOLP());
  cout<<endl;

}

//---------------------------------------------------
// Propagation of a Jones vector from CIP to Exit line 
// through 2 mirrors modelized by a phase shift 
// (Delta in rad) + a rotation of the ellipse (Theta in rad)
//---------------------------------------------------
void VJones::Propa(Double_t Theta1,Double_t Delta1,Double_t Theta2,Double_t Delta2)
{
  Rotate(Theta1);
  PhaseShift(Delta1);
  Rotate(Theta2);
  PhaseShift(Delta2);
}

//---------------------------------------------------
// Global parametrization with only 1 phase and 1 rotation
//---------------------------------------------------
void VJones::Propa(Double_t Theta,Double_t Delta)
{
  Rotate(Theta);
  PhaseShift(Delta);
}

//---------------------------------------------------
// Backward propagation of a Jones vector from Exit line 
// to CIP through 2 mirrors modelized by a phase shift 
// (Delta in rad) + a rotation of the ellipse (Theta in rad)
//---------------------------------------------------
void VJones::PropaBack(Double_t Delta2,Double_t Theta2,Double_t Delta1,Double_t Theta1)
{
  PhaseShift(-Delta2);
  Rotate(-Theta2);
  PhaseShift(-Delta1);
  Rotate(-Theta1);
}

//---------------------------------------------------
//Global parametrization with only 1 phase and 1 rotation
//---------------------------------------------------
void VJones::PropaBack(Double_t Delta,Double_t Theta)
{
  PhaseShift(-Delta);
  Rotate(-Theta);
}


