#ifndef MJones_h
#define MJones_h

#include "VJones.h"

//---------------------------------------------------
// Jones Matrix
//---------------------------------------------------
struct MJones {
	TComplex m[2][2];
};

MJones MJonesInv(MJones M);
MJones MJonesTranspose(MJones M);
MJones MJonesRotator(Double_t Alpha);
MJones MJonesProd(MJones M1, MJones M2);
MJones MJonesPlate(Double_t gamma);
MJones MJonesPolarizer(Double_t Theta);
MJones MJonesPartPol(Double_t Pol);
MJones MJonesPartPol2(Double_t Pol, Double_t Pol2);
MJones MJonesRPartPol(Double_t Pol, Double_t Theta);
MJones MJonesRPlate(Double_t Phi, Double_t Theta);
MJones MJonesGenBirefringence(Double_t A, Double_t gamma, Double_t B);
MJones MJonesGenSystem(Double_t A1, Double_t gamma1, Double_t B1, Double_t Pol, 
					   Double_t A2, Double_t gamma2, Double_t B2, Double_t Pol2);
VJones MtimesVJones(MJones M, VJones* V);
MJones MJonesMirror(Double_t Phi, Double_t Theta, Double_t Alpha);
VJones TF(VJones* E, MJones M1, MJones M2);
void MJonesPrint(MJones M);
void MJonesPrintComplex(MJones M);

#endif


/* emacs
 * Local Variables:
 * mode:C++
 * mode:font-lock
 * c-file-style: "stroustrup"
 * tab-width: 4
 * End:
 */

