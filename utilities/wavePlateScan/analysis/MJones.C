#include <TComplex.h>
#include "VJones.h"
#include "MJones.h"

//---------------------------------------------------
// Inverse of a Jones Matrix .
//---------------------------------------------------
MJones MJonesInv(MJones M)
{
	MJones W;
  
	TComplex Det,DetInv,ad,bc;
  
	ad = M.m[0][0]*M.m[1][1]; 
	bc = M.m[0][1]*M.m[1][0];
  
	Det = ad; Det -=bc;
	DetInv = TComplex(cos(Det.Theta())/Det.Rho(),-sin(Det.Theta())/Det.Rho());
  
	W.m[0][0] = M.m[1][1]*DetInv;
	W.m[0][1] = M.m[0][1]*DetInv; W.m[0][1]*=-1;
	W.m[1][0] = M.m[1][0]*DetInv; W.m[1][0]*=-1;
	W.m[1][1] = M.m[0][0]*DetInv;
  
	return W;
}

MJones MJonesTranspose(MJones M)
{
	MJones W;
	W.m[0][0] = M.m[0][0];
	W.m[0][1] = M.m[1][0];
	W.m[1][0] = M.m[0][1];
	W.m[1][1] = M.m[1][1];
	return W;
}

//---------------------------------------------------
// Jones Matrix of a rotator. Rotates the polarization 
// ellipse by an angle Alpha (in rad) 
//---------------------------------------------------
MJones MJonesRotator(Double_t Alpha)
{
	MJones M;
  
	M.m[0][0] = TComplex(cos(Alpha),0.);
	M.m[0][1] = TComplex(-sin(Alpha),0.);
	M.m[1][0] = TComplex(sin(Alpha),0.);
	M.m[1][1] = TComplex(cos(Alpha),0.);
  
	return M;
}

//---------------------------------------------------
// Multiplication of two Jones Matrices
//---------------------------------------------------
// Note that this is  M2.M1 as expected
MJones MJonesProd(MJones M1, MJones M2)
{
	TComplex CProd;
	MJones M;
  
	for (Int_t i=0; i<2; i++)
		for (Int_t j=0; j<2; j++)
			M.m[i][j] = TComplex(0,0);
  
	for (Int_t i=0; i<2; i++)
		for (Int_t j=0; j<2; j++)
			for (Int_t k=0; k<2; k++)
			{
				CProd = M1.m[i][k]*M2.m[k][j];
				M.m[i][j] += CProd;
			}
  
	return M;
}


// This routine give a vertical fast axis ??
// routine altered MMD
MJones MJonesPlate(Double_t gamma)
{
	MJones M;
  	M.m[0][0] = TComplex(cos(gamma),sin(gamma));
	M.m[0][1] = 0;
	M.m[1][0] = 0;
	M.m[1][1] = TComplex(cos(gamma),-1*sin(gamma));
	return M;
}

MJones MJonesPolarizer(Double_t Theta)
{  
	MJones M;
  	M.m[0][0] = 1;
	M.m[0][1] = 0;
	M.m[1][0] = 0;
	M.m[1][1] = 0;
	MJones R = MJonesRotator(Theta);
	MJones Rm = MJonesRotator(-1*Theta);
 	M = MJonesProd(M,R);
 	M = MJonesProd(Rm,M);
	return M;
}

MJones MJonesPartPol(Double_t Pol)
{  
	MJones M;
  	M.m[0][0] = 1;
	M.m[0][1] = 0;
	M.m[1][0] = 0;
	M.m[1][1] = 1-Pol;
	return M;
}


MJones MJonesPartPol2(Double_t Pol, Double_t Pol2)
{  
	MJones M;
  	M.m[0][0] = 1-Pol;
	M.m[0][1] = 0;
	M.m[1][0] = 0;
	M.m[1][1] = 1-Pol2;
	return M;
}

MJones MJonesRPartPol(Double_t Pol, Double_t Theta)
{  
	MJones M;
  	M.m[0][0] = 1;
	M.m[0][1] = 0;
	M.m[1][0] = 0;
	M.m[1][1] = 1-Pol;
	MJones R = MJonesRotator(Theta);
	MJones Rm = MJonesRotator(-1*Theta);
 	M = MJonesProd(M,R);
 	M = MJonesProd(Rm,M);
	return M;
}



//---------------------------------------------------
// Jones Matrix of a retardation plate with phase-shift
// Phi (in rad) and slow axis at an angle Theta (in rad) 
// w.r.t. the X axis of the current coordinate frame.
//---------------------------------------------------
MJones MJonesRPlate(Double_t Phi, Double_t Theta)
{
	MJones M = MJonesPlate(Phi);
	MJones R = MJonesRotator(Theta);
	MJones Rm = MJonesRotator(-1*Theta);
 	M = MJonesProd(M,R);
 	M = MJonesProd(Rm,M);
	return M;
}



//---------------------------------------------------
// Jones Matrix of a general collection of birefringent 
// elements as described in section 1 of Jones II
// http://www.opticsinfobase.org/abstract.cfm?URI=josa-31-7-488
// symbols as defined in paper
//---------------------------------------------------
MJones MJonesGenBirefringence(Double_t A, Double_t gamma, Double_t B)
{
	//  RB.M.RA
	MJones M = MJonesPlate(gamma);
	MJones RA = MJonesRotator(A);
	MJones RB = MJonesRotator(B);
 	M = MJonesProd(M,RA);
 	M = MJonesProd(RB,M);
	return M;
}

//---------------------------------------------------
// Jones Matrix of a general collection of birefringent and polarizing 
// elements as described in equation 24 of Jones II
// http://www.opticsinfobase.org/abstract.cfm?URI=josa-31-7-488
// symbols as defined in paper
//---------------------------------------------------
MJones MJonesGenSystem(Double_t A1, Double_t gamma1, Double_t B1, Double_t Pol, 
					   Double_t A2, Double_t gamma2, Double_t B2, Double_t Pol2)
{
	//   B2.g2.A2.P.B1.g1.A1
	MJones Mg1 = MJonesPlate(gamma1);
	MJones Mg2 = MJonesPlate(gamma2);
	MJones MRA1 = MJonesRotator(A1);
	MJones MRB1 = MJonesRotator(B1);
	MJones MRA2 = MJonesRotator(A2);
	MJones MRB2 = MJonesRotator(B2);
	MJones MP = MJonesPartPol2(Pol, Pol2);
	MJones Mprod = MJonesProd(Mg2,MRB2);
	Mprod = MJonesProd(MRA2,Mprod);
	Mprod = MJonesProd(MP,Mprod);
	Mprod = MJonesProd(MRB1,Mprod);
	Mprod = MJonesProd(Mg1,Mprod);
	Mprod = MJonesProd(MRA1,Mprod);
	return Mprod;
}






//---------------------------------------------------
// Multiplication of a Jones vector by a Jones Matrix
//---------------------------------------------------
VJones MtimesVJones(MJones M, VJones* V)
{
	VJones MV;
	TComplex tc[2];

	for (Int_t i=0; i<2; i++)
		tc[i]=TComplex(0,0);

	for (Int_t i=0; i<2; i++)
		for (Int_t j=0; j<2; j++)
			tc[i] += M.m[i][j]*V->GetCmplx(j);

	MV = VJones(tc[0],tc[1]);
	return MV;
}

//---------------------------------------------------
// Jones Matrix of a mirror. Phi and Theta parametrize 
// the phaseshift induced by the birefringence of the 
// mirror. Alpha is the ellipse rotation induced by 
// the orientation of the mirror.
// We choose M = RP*R, which is NOT equivalent to R*RP
//---------------------------------------------------
MJones MJonesMirror(Double_t Phi, Double_t Theta, Double_t Alpha)
{
	MJones RP=MJonesRPlate(Phi,Theta);
	MJones R=MJonesRotator(Alpha);
 
	return MJonesProd(RP,R);
}

//---------------------------------------------------
// New modelization of Transfer Function
//---------------------------------------------------
VJones TF(VJones* E, MJones M1, MJones M2)
{
	return MtimesVJones(MJonesProd(M1,M2),E);
}

//---------------------------------------------------
// Print out of a Jones matrix
//---------------------------------------------------
void MJonesPrint(MJones M)
{  
	printf("| (%+.4f,%+7.2f), (%+.4f,%+7.2f)|\n",  
		   M.m[0][0].Rho(),M.m[0][0].Theta()*deg,M.m[0][1].Rho(),M.m[0][1].Theta()*deg);
	printf("| (%+.4f,%+7.2f), (%+.4f,%+7.2f)|\n\n",
		   M.m[1][0].Rho(),M.m[1][0].Theta()*deg,M.m[1][1].Rho(),M.m[1][1].Theta()*deg);
}

//---------------------------------------------------
// Print out of a Jones matrix in complex form
//---------------------------------------------------
void MJonesPrintComplex(MJones M)
{  
	printf("| %+.4f %+.4f i, %+.4f %+.4f i|\n",  
		   M.m[0][0].Re(),M.m[0][0].Im(),M.m[0][1].Re(),M.m[0][1].Im());
	printf("| %+.4f %+.4f i, %+.4f %+.4f i|\n\n",
		   M.m[1][0].Re(),M.m[1][0].Im(),M.m[1][1].Re(),M.m[1][1].Im());
}



/* emacs
 * Local Variables:
 * mode:C++
 * mode:font-lock
 * c-file-style: "stroustrup"
 * tab-width: 4
 * End:
 */

