
#ifndef VJones_h
#define VJones_h


#include "Rtypes.h"
#include <TComplex.h>
#include <iostream>

const Double_t PI=acos(-1.);
const Double_t deg=180./PI; 
const Double_t rad=PI/180.;


class VJones
{
//---------------------------------------------------
// Jones representation of oscillating electric field
// 2-Vector of complex numbers
//---------------------------------------------------

 public:
  
  VJones (TComplex x,TComplex y); 
  VJones (); 
  virtual ~VJones();
  
  void Init(TComplex,TComplex);
  void InitNorm(TComplex ,TComplex );
  TComplex GetCmplx(Int_t);
  Double_t Norm();
  Double_t Phase();
  void PhaseShift(Double_t );
  void GlobalPhase(Double_t );
  void Rotate(Double_t );
  Double_t Alpha();
  void ProjCirc();
  Double_t DOCP();
  Double_t DOLP();
  void Print();
  void Propa(Double_t ,Double_t ,Double_t ,Double_t );
  void Propa(Double_t ,Double_t );
  void PropaBack(Double_t ,Double_t ,Double_t ,Double_t );
  void PropaBack(Double_t ,Double_t );


 
 protected:

  TComplex fv1;
  TComplex fv2;
  Double_t PI;
  Double_t deg;
  Double_t rad;
/*   const Double_t PI=acos(-1.); */
/*   const Double_t deg=180./PI; */
/*   const Double_t rad=PI/180.; */
  
 private:


  ClassDef(VJones,0);

};


#endif

//---------------------------------------------------


