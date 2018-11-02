#include <TComplex.h>
#include "VJones.h"
#include "MJones.C"

//-----------------------------------------------------------
// CONVENTIONS:
//
// Jones vector of Left Polarization state  : (1,i)/sqrt(2) 
// Jones vector of Right Polarization state : (1,-i)/sqrt(2)
// Linear Polarization along X (1,0), along Y (0,1)
//
// XYZ is a direct coordinate frame. We choose the following 
// orientations at the CIP AND in the exit line:
//      Z is along the propagation of the photon
// 	X is vertical upward
// 	Y is horizontal to the right when looking along Z
//
//---------------------------------------------------
// Test routine to debug operations with Jones vectors 
// and matrices
//---------------------------------------------------

Double_t HWPoffset=10;
Double_t QWPoffset=0;

// this function tests a rotating QWP
Double_t testF1(Double_t x) {
	VJones H = VJones(1,0);
	MJones M = MJonesRPlate(PI/4,x*rad);
	H = MtimesVJones(M,&H);
	return H.DOCP();
}

// this function tests a rotating mirror
Double_t test_mir1(Double_t x) {
	VJones Rc = VJones(1,0);//TComplex(0,-1));
	MJones Mcav1 = MJonesMirror(PI/2,x,0);
	Rc = MtimesVJones(Mcav1,&Rc);
	return Rc.DOCP();
}

MJones enter(Double_t hwpang, Double_t qwpang) {
	MJones Mhwp = MJonesRPlate(PI/2,(hwpang-HWPoffset)*rad);
	MJones Mqwp = MJonesRPlate(PI/4,(qwpang-QWPoffset)*rad);
	MJones MEnter = MJonesProd(Mqwp,Mhwp);
	MJones Mystery = MJonesRPlate(PI/10,PI/10);
	MEnter = MJonesProd(Mystery,MEnter);
	return MEnter;
}

// The most basic functions
VJones propagate(VJones init,Double_t hwpang,Double_t qwpang) {
	MJones Mtot;
	MJones MEnter = enter(hwpang,qwpang);
	MJones MEnterInv = MJonesTranspose(MEnter);
	Mtot = MJonesProd(MEnterInv,MEnter);
	VJones fin = MtimesVJones(Mtot,&init);
	return fin;
}

VJones propagate_half(VJones init,Double_t hwpang,Double_t qwpang) {
	MJones Mtot;
	MJones MEnter = enter(hwpang,qwpang);
	VJones fin = MtimesVJones(MEnter,&init);
	return fin;
}

// VJones propagate(VJones init,Double_t x,Double_t y) {
// 	MJones Mtot;
// 	MJones Mhwp = MJonesRPlate(PI/2,(x-HWPoffset)*rad);
// 	MJones Mqwp = MJonesRPlate(PI/4,(y-QWPoffset)*rad);
// 	MJones MEnter = MJonesProd(Mqwp,Mhwp);
// 	MJones Mystery = MJonesRPlate(0,0);
// 	MJones MEnter = MJonesProd(Mystery,MEnter);
// 	MJones MEnterInv = MJonesTranspose(MEnter);
// //	MJones MEnterInv = MJonesInv(MEnter);
// 	MJones Mcav1 = MJonesMirror(PI/2,0,0); // 2nd parameter irrelevant for PI/2 
// 	Mtot = MJonesProd(Mcav1,MEnter);
// 	Mtot = MJonesProd(MEnterInv,MEnter);
// 	VJones fin = MtimesVJones(Mtot,&init);
// 	return fin;
// }

// VJones propagate_half(VJones init,Double_t x,Double_t y) {
// 	MJones Mtot;
// 	MJones Mhwp = MJonesRPlate(PI/2,(x-HWPoffset)*rad);
// 	MJones Mqwp = MJonesRPlate(PI/4,(y-QWPoffset)*rad);
// 	MJones MEnter = MJonesProd(Mhwp,Mqwp);
// //	MJones Mystery = MJonesRPlate(PI/10,-PI/10);
// 	MJones Mystery = MJonesRPlate(0,0);
//     MEnter = MJonesProd(Mystery,MEnter);
// 	VJones fin = MtimesVJones(MEnter,&init);
// 	return fin;
// }

// This function explores an eigen generator
Double_t testDOCP(Double_t hwpang,Double_t qwpang) {
	VJones H = VJones(1,0);//TComplex(0,-1));
	H = propagate(H,hwpang,qwpang);
	return H.DOCP();
}

// This function explores an eigen generator
Double_t testDOLP(Double_t hwpang,Double_t qwpang) {
	VJones H = VJones(1,0);//TComplex(0,-1));
	H = propagate(H,hwpang,qwpang);
	return H.DOLP();
}

// This function explores an eigen generator
Double_t testCavDOCP(Double_t hwpang,Double_t qwpang) {
	VJones H = VJones(1,0);//TComplex(0,-1));
	H = propagate_half(H,hwpang,qwpang);
	return H.DOCP();
}

// This function explores an eigen generator
Double_t testAngle(Double_t hwpang,Double_t qwpang) {
	VJones H = VJones(1,0);//TComplex(0,-1));
	H = propagate(H,hwpang,qwpang);
	return H.Alpha()*deg;
}

Double_t testAngleVert(Double_t hwpang,Double_t qwpang) {
	VJones H = VJones(1,0);
	H = propagate(H,hwpang,qwpang);
	MJones Mvert = MJonesPolarizer(PI/2);
	H = MtimesVJones(Mvert,&H);
	return H.Alpha()*deg;
}

Double_t testAngleHoriz(Double_t hwpang,Double_t qwpang) {
	VJones H = VJones(1,0);
	H = propagate(H,hwpang,qwpang);
	MJones Mhoriz = MJonesPolarizer(0);
	H = MtimesVJones(Mhoriz,&H);
	return H.Alpha()*deg;
}

Double_t testPower(Double_t hwpang,Double_t qwpang) {
	VJones H = VJones(1,0);//TComplex(0,-1));
	H = propagate(H,hwpang,qwpang);
	return H.Norm();
}

Double_t testPowerVert(Double_t hwpang,Double_t qwpang) {
	VJones H = VJones(1,0);//TComplex(0,-1));
	H = propagate(H,hwpang,qwpang);
	MJones Mvert = MJonesPolarizer(PI/2);
	H = MtimesVJones(Mvert,&H);
	return H.Norm();
}

Double_t testPowerHoriz(Double_t hwpang,Double_t qwpang) {
	VJones H = VJones(1,0);//TComplex(0,-1));
	H = propagate(H,hwpang,qwpang);
	MJones Mhoriz = MJonesPolarizer(0);
	H = MtimesVJones(Mhoriz,&H);
	return H.Norm();
}

void Test()
{
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadBottomMargin(0.15);
	

// 	VJones *H = new VJones(1,0);
// 	H->Print();
// 	H->Rotate(90*rad);
// 	H->Print();
	
// 	VJones *Rc = new VJones(1,TComplex(0,-1));
// //	Rc->Print();
// 	MJones Mhwp1 = MJonesRPlate(PI/2,0);
// 	MJones Mqwp1 = MJonesRPlate(PI/2,0);
// 	MJones Mcav1 = MJonesMirror(PI/2,0,0); // 2nd parameter irrelevant for PI/2 
// 	MJones Mtot = MJonesProd(Mhwp1,Mqwp1);
// 	MJonesPrint(Mcav1);
// 	MJonesPrintComplex(Mcav1);

	MJones Mvert = MJonesPolarizer(PI/2);
	MJonesPrint(Mvert);


//	VJones G = MtimesVJones(Mcav1,Rc);
//	G.Print();

// 	TF1 *test = new TF1("test","testF1(x)",0,360);
// 	test->SetTitle("DOCP off mirror");
// 	test->Draw();

// 	TF1 *test = new TF1("test","testF1(x)",0,360);
// 	test->SetTitle("DOCP off mirror");
// 	test->Draw();

// 	Double_t HWPmin = 180;
// 	Double_t HWPmax = 290;
// 	Double_t QWPmin = 200;
// 	Double_t QWPmax = 360;

	Double_t HWPmin = 0;
	Double_t HWPmax = 360;
	Double_t QWPmin = 0;
	Double_t QWPmax = 360;

	TCanvas *canvas = new TCanvas("canvas","gen canvas",20,20,1400,600);
	canvas->Divide(5,2,0.001,0.001);
	canvas->cd(1);
	TF2 *test = new TF2("test","testDOCP(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	test->SetTitle("DOCP at polarizer;HWP;QWP");
	test->Draw("colz");
	canvas->cd(2);
	TF2 *testDOLP = new TF2("testDOLP","testDOLP(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testDOLP->SetTitle("DOLP at polarizer;HWP;QWP");
	testDOLP->Draw("colz");
	canvas->cd(3);
	TF2 *testAngle = new TF2("testAngle","testAngle(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testAngle->SetTitle("Angle at polarizer;HWP;QWP");
	testAngle->Draw("colz");
	canvas->cd(3);
// 	TF2 *testPV = new TF2("testPV","testPowerVert(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
// 	testPV->SetTitle("Power after polarizer (vert);HWP;QWP");
// 	testPV->Draw("colz");
	canvas->cd(4);
	TF2 *testPH = new TF2("testPH","testPowerHoriz(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testPH->SetTitle("Power after polarizer (horiz);HWP;QWP");
	testPH->Draw("colz");
	canvas->cd(5);
	TF2 *testCavDOCP = new TF2("testCavDOCP","testCavDOCP(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testCavDOCP->SetTitle("DOCP at cavity;HWP;QWP");
	testCavDOCP->Draw("colz");


	canvas->cd(6);
	gPad->SetTheta(0);
	gPad->SetPhi(35);
//	TF2 *test = new TF2("test","testDOCP(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	test->SetTitle("DOCP at polarizer;HWP;QWP");
	test->Draw("surf1");
	canvas->cd(7);
	gPad->SetTheta(0);
	gPad->SetPhi(35);
//	TF2 *testDOLP = new TF2("testDOLP","testDOLP(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testDOLP->SetTitle("DOLP at polarizer;HWP;QWP");
	testDOLP->Draw("surf1");
	canvas->cd(8);
	gPad->SetTheta(0);
	gPad->SetPhi(35);
//	TF2 *testAngle = new TF2("testAngle","testAngle(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testAngle->SetTitle("Angle at polarizer;HWP;QWP");
	testAngle->Draw("surf1");
// 	canvas->cd(8);
// 	gPad->SetTheta(0);
// 	gPad->SetPhi(35);
// 	testPV->SetTitle("Power after polarizer (vert);HWP;QWP");
// 	testPV->Draw("surf1");
	canvas->cd(9);
	gPad->SetTheta(0);
	gPad->SetPhi(35);
	testPH->SetTitle("Power after polarizer (horiz);HWP;QWP");
	testPH->Draw("surf1");
	canvas->cd(10);
	gPad->SetTheta(0);
	gPad->SetPhi(35);
//	TF2 *testCavDOCP = new TF2("testCavDOCP","testCavDOCP(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testCavDOCP->SetTitle("DOCP at cavity;HWP;QWP");
	testCavDOCP->Draw("surf1");


	Int_t padnum=1;
	TCanvas *powcanvas = new TCanvas("powcanvas","pow canvas",20,20,1400,600);
	powcanvas->Divide(3,2,0.001,0.001);
	powcanvas->cd(padnum); padnum++;
	TF2 *testPB = new TF2("testPB","testPower(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testPB->SetTitle("Power before polarizer;HWP;QWP");
	testPB->Draw("surf1");
	powcanvas->cd(padnum); padnum++;
	gPad->SetTheta(0);
	gPad->SetPhi(35);
	TF2 *testPV = new TF2("testPV","testPowerVert(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testPV->SetTitle("Power after polarizer (vert);HWP;QWP");
	testPV->Draw("surf1");
	powcanvas->cd(padnum); padnum++;
	gPad->SetTheta(0);
	gPad->SetPhi(35);
	TF2 *testPH = new TF2("testPH","testPowerHoriz(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testPH->SetTitle("Power after polarizer (horiz);HWP;QWP");
	testPH->Draw("surf1");
	powcanvas->cd(padnum); padnum++;
	gPad->SetTheta(0);
	gPad->SetPhi(35);
	TF2 *testPSum = new TF2("testPSum","testPowerHoriz(x,y)+testPowerVert(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testPSum->SetTitle("Power after polarizer (sum);HWP;QWP");
	testPSum->Draw("surf1");
	powcanvas->cd(padnum); padnum++;
	TF2 *testAngleV = new TF2("testAngleV","testAngleVert(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testAngleV->SetTitle("Angle after polarizer (vert);HWP;QWP");
	testAngleV->Draw("surf1");
	powcanvas->cd(padnum); padnum++;
	TF2 *testAngleH = new TF2("testAngleH","testAngleHoriz(x,y)",HWPmin,HWPmax,QWPmin,QWPmax);
	testAngleH->SetTitle("Angle after polarizer (horiz);HWP;QWP");
	testAngleH->Draw("surf1");
}



/* emacs
 * Local Variables:
 * mode:C++
 * mode:font-lock
 * c-file-style: "stroustrup"
 * tab-width: 4
 * End:
 */

 
