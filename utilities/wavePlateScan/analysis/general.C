#include <TComplex.h>
#include "VJones.h"
#include "MJones.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TF2.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"


#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>

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

Int_t numparams=0;
const Int_t maxparams=12;
const Int_t NMAX=35000;  //absolute maximum number of points

// The measured data is stored in 4 arrays
Double_t HWP[NMAX], QWP[NMAX], RRPD[NMAX], dRRPD[NMAX];
Int_t numdatapoints;
Int_t numfcncalls;

Double_t gDOCP, gDOLP, gangle;
Double_t params[maxparams], dparams[maxparams];
//Double_t polarizer_angle = -0.1483; // exit line analyzer
//Double_t polarizer_angle = 0;     // this is for the RRPD
Double_t polarizer_angle = PI/2; // this is for the RPD


// In this scenario, the model power is scaled to match the ADC values of the data
// Also the model HWP and QWP angles are made to mimic that of the data 
// The 4 parameters are those required to accomplish this change but the actual 
// change is done in the fcn routine and the Jones algebra remains pure


void drawdata(TString,	TString, TString, Double_t *, Int_t );
void printmodels();
Int_t getnumparams(Int_t modelnum);
Int_t general_setup(TMinuit *, Int_t, Int_t detector, TString &);
Double_t getmodelpower(Double_t xp, Double_t yp, Double_t *params, Int_t modelnum);

// model 1
Double_t RRPD_platerr_bir(Double_t , Double_t , Double_t, Double_t, Double_t , Double_t , Double_t , Bool_t );
Double_t getDOLP_platerr_bir(Double_t *x, Double_t *par); // used to maximize DOCP
Double_t get_modelpower_platerr_bir(Double_t xp, Double_t yp, Double_t *params); 
void fcn_platerr_bir(Int_t &npar, Double_t *gin, Double_t &returnval, Double_t *par, Int_t iflag);
void setup_platerr_bir(TMinuit *gMinuit, Int_t detector);

// model 2
Double_t RRPD_plate_bir_pol(Double_t HWPang, Double_t QWPang, Double_t A, Double_t gamma,  Double_t B, 
							Double_t HWPfac, Double_t QWPfac, Double_t Pol, Double_t Theta, Bool_t moredata=0);
Double_t getDOLP_plate_bir_pol(Double_t *x, Double_t *par); // used to maximize DOCP
void fcn_plate_bir_pol(Int_t &npar, Double_t *gin, Double_t &returnval, Double_t *par, Int_t iflag);
Double_t get_modelpower_plate_bir_pol(Double_t xp, Double_t yp, Double_t *params) ;
void setup_plate_bir_pol(TMinuit *gMinuit, Int_t detector); 

// model 3
Double_t RRPD_anal(Double_t , Double_t , Double_t, Double_t, Double_t , Double_t , Double_t , Double_t , Bool_t );
Double_t getDOLP_anal(Double_t *x, Double_t *par); // used to maximize DOCP
Double_t get_modelpower_anal(Double_t xp, Double_t yp, Double_t *params); 
void fcn_anal(Int_t &npar, Double_t *gin, Double_t &returnval, Double_t *par, Int_t iflag);
void setup_anal(TMinuit *gMinuit, Int_t detector);





//---------------------------------------------------
// Read data of reflected power measurements and fit 
// the parameters, phases and rotations
// Print results
//---------------------------------------------------
void Minimize(TString DataName="1", Int_t modelnum=-1, Int_t detector=1, TString addcomment="", 
			  TString runpointsfile="data/runpoints.txt",
			  TString TF_Directory="/home/compton/gaskelld/entrancefunction/TransferFunction/data/output", 
			  TString ActualFile="AllData.txt")
{
	Bool_t fail=0;
	if (modelnum==-1) {
		printf("need a model number\n");
 		fail=1;
	}
	if (DataName=="1") {
		printf("need arguments\n");
		fail=1;
	}
	if (fail) {
 		printf("\nuse:\n\tMinimize(DataName, modelnum, [detector], [comment], [TF_Directory], [ActualFile])\n");
		printf("defaults:\n");
		printf("\tdetector     = %i  (1=RPD, 2=RRPD, 3=analyzer)\n", detector);
		printf("\taddcomment   = \"%s\"  (optional additional comment)\n", addcomment.Data());
		printf("\tTF_Directory = %s\n", TF_Directory.Data());
		printf("\tActualFile   = %s\n\n", ActualFile.Data());
		printf("inputfile = TF_Directory/DataName/ActualFile\n\n");

		printmodels();
 		return;
 	} 

	TString detectorname;
	if (detector==1) {
		polarizer_angle=PI/2;
		detectorname="RPD";
	} else {
		if (detector==2){
			polarizer_angle=PI/2;
			detectorname="RRPD";
		} else {
			if (detector==3){
				polarizer_angle=0;
				detectorname="Analyzer";
			} else {
				printf("bad detector number\n");
			}
		}
	}
	printf("using %s,polarizer_angle=%f\n",detectorname.Data(),polarizer_angle);

	TString infilename = TF_Directory+"/"+DataName+"/"+ActualFile;
	printf("trying to open %s\n",infilename.Data());

	numdatapoints=0;
	numfcncalls=0;

	TString comment;
	Int_t counter=0;
	Double_t x1,x2,x3,x4;
	Double_t RRPDmin=0, RRPDmax=0;
	
	ifstream in(infilename.Data());
	if (!in.is_open())
	{
		printf("Parameter file not found: %s\n",infilename.Data());
		return;
	} else {
		printf("reading data...\n");
		while (counter<NMAX) {
			// Read data file
			in>>x1>>x2>>x3>>x4;
			if (!in.good()) break;
			QWP[counter]=x1; HWP[counter]=x2; RRPD[counter]=x3; dRRPD[counter]=x4; 
			if (counter==0) {
				RRPDmin=RRPD[counter];
				RRPDmax=RRPD[counter];
			} else {
				if (RRPDmin>RRPD[counter]) RRPDmin=RRPD[counter];
				if (RRPDmax<RRPD[counter]) RRPDmax=RRPD[counter];
			}
			printf("%8.2f %8.2f %10.2f %6.2f\n", HWP[counter], QWP[counter], RRPD[counter], dRRPD[counter]);
			counter++;
		}
		in.close();
	}
	numdatapoints = counter;
	printf("read %i lines.\n",numdatapoints);
	printf("min = %.2f, max = %.2f, diff=%.2f\n",RRPDmin,RRPDmax,RRPDmax-RRPDmin);

	TStopwatch timer;
	printf("Starting timer\n");
	timer.Start();

	numfcncalls=0;
	printf("numdatapoints=%i,  numfcncalls=%i\n",numdatapoints,numfcncalls);
	//--------------------------
	//MINIMIZATION USING TMinuit
	//--------------------------
  
	//Initialize TMinuit with a maximum number of params
	numparams = getnumparams(modelnum);
	printf("There are %i parameters\n",numparams);
	TMinuit *gMinuit = new TMinuit(numparams);
	// Set the minimization function.
	// use external diambuguation routine.
	if (general_setup(gMinuit, modelnum, detector, comment) <0) {
		printf("bad model number\n");
		return;
	}
	comment += addcomment;
	//Array and flag for minuit parameters
	Double_t arglist[10];
	Int_t ierflg = 0;
	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Now perform the minimization
	Int_t maxsteps=20000;
	arglist[0] = maxsteps;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	// Now write the output
	FILE *outFile, *resFile, *qwpresFile, *hwpresFile, *modelFile, *datacalFile, 
		*paramFile, *docpFile, *dolpFile, *angleFile, *DOCPcorrelFile;
	TString plotsbasestr = Form("%s_%s",DataName.Data(),comment.Data());
	char command[1255], outdir[1255];
	int result;
	sprintf(outdir,"output/%s",plotsbasestr.Data());
	sprintf(command,"mkdir %s",outdir);
	result=system(command);

//	TString filebasestring = Form("%s/%s",outdir,plotsbasestr.Data());
	TString filebasestring = Form("%s/%s",outdir,comment.Data());

	TString outfilename=Form("%s.out",filebasestring.Data());
	outFile = fopen (outfilename.Data(),"w");
	TString resfilename=Form("%s.residuals.out",filebasestring.Data());
	resFile = fopen (resfilename.Data(),"w");
	TString qwpresfilename=Form("%s.qwpresiduals.out",filebasestring.Data());
	qwpresFile = fopen (qwpresfilename.Data(),"w");
	TString hwpresfilename=Form("%s.hwpresiduals.out",filebasestring.Data());
	hwpresFile = fopen (hwpresfilename.Data(),"w");
	TString modelfilename=Form("%s.model.out",filebasestring.Data());
	modelFile = fopen (modelfilename.Data(),"w");
	TString datacalfilename=Form("%s.datacal.out",filebasestring.Data());
	datacalFile = fopen (datacalfilename.Data(),"w");
	TString paramfilename=Form("%s.params.out",filebasestring.Data());
	paramFile = fopen (paramfilename.Data(),"w");
	TString modelDOCPfilename=Form("%s.docp.out",filebasestring.Data());
	docpFile = fopen (modelDOCPfilename.Data(),"w");
	TString modelDOLPfilename=Form("%s.dolp.out",filebasestring.Data());
	dolpFile = fopen (modelDOLPfilename.Data(),"w");
	TString LPanglefilename=Form("%s.angle.out",filebasestring.Data());
	angleFile = fopen (LPanglefilename.Data(),"w");
	TString DOCPcorrelfilename=Form("%s.docpvsrrpd.out",filebasestring.Data());
	DOCPcorrelFile = fopen (DOCPcorrelfilename.Data(),"w");


// Int_t Eval(Int_t npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t flag)

//  Evaluate the minimisation function
//   Input parameters:
//     npar:    number of currently variable parameters
//     par:     array of (constant and variable) parameters
//     flag:    Indicates what is to be calculated (see example below)
//     grad:    array of gradients
//   Output parameters:
//     fval:    The calculated function value.
//     grad:    The (optional) vector of first derivatives).

	for (Int_t i=0; i<numparams; i++) {
		gMinuit->GetParameter(i,params[i],dparams[i]);
	}
	Double_t evalval(0);
	Double_t* grad=NULL;
	//Int_t evalint = 
	gMinuit->Eval(gMinuit->GetNumFreePars(), grad, evalval, params, 1);

	printf("Writing parameters.\nDouble_t vstart[] = {");
	fprintf(paramFile,"chi2  %f\n", evalval);
	for (Int_t i=0; i<numparams; i++) {
		fprintf (paramFile, "%i %10.4f +- %8.4f\n", i, params[i],dparams[i]);
		printf ("%.4f", params[i]);	
		if (i<numparams-1) printf (", ");	
	}
	printf ("};   // chi2 = %.1f \nDouble_t vstep[] = {",evalval);	
	for (Int_t i=0; i<numparams; i++) {
		printf ("%7.4f", dparams[i]);
		if (i<numparams-1) printf (", ");		
	}
	printf ("};\n");
	
	printf ("Errors:  eplus,  eminus,  eparab,  gcc\n");
	for (Int_t i=0; i<numparams; i++) {
		Double_t eplus,  eminus,  eparab,  gcc;
		gMinuit->mnerrs(i,  eplus,  eminus,  eparab,  gcc);
		fprintf (paramFile,"%i  %.4f  %.4f  %.4f  %.4f \n", i,  eplus,  eminus,  eparab,  gcc);
		printf ("%2i       %.4f  %.4f  %.4f  %.4f \n", i,  eplus,  eminus,  eparab,  gcc);	
	}


	printf("Writing %i data points to output files.\n", numdatapoints);
	for (Int_t i=0;i<numdatapoints; i++) 
	{

 		Double_t modelpower = getmodelpower(HWP[i],QWP[i],params,modelnum);

		fprintf (outFile, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.6f  %12.6f  %12.4f  \n", 
				 HWP[i], QWP[i],  RRPD[i], dRRPD[i], modelpower, RRPD[i] - modelpower, gDOCP, gDOLP, gangle);
		fprintf (resFile, "%12.3f  %12.3f  %12.6f\n", 
				 HWP[i], QWP[i], RRPD[i] - modelpower);
		fprintf (qwpresFile, "%12.3f  %12.6f  %12.6f\n", 
				 QWP[i], RRPD[i] - modelpower,dRRPD[i]);
		fprintf (hwpresFile, "%12.3f  %12.6f  %12.6f\n", 
				 HWP[i], RRPD[i] - modelpower,dRRPD[i]);
		fprintf (modelFile, "%12.3f  %12.3f  %12.6f\n", 
				 HWP[i], QWP[i], modelpower);
		fprintf (datacalFile, "%12.3f  %12.3f  %12.6f\n", 
				 HWP[i], QWP[i], RRPD[i]);
		fprintf (docpFile, "%12.3f  %12.3f  %12.6f\n", 
				 HWP[i], QWP[i], gDOCP);
		fprintf (dolpFile, "%12.3f  %12.3f  %12.6f\n", 
				 HWP[i], QWP[i], gDOLP);
		fprintf (angleFile, "%12.3f  %12.3f  %12.6f\n", 
				 HWP[i], QWP[i], gangle);
		fprintf (DOCPcorrelFile, "%12.4f  %12.4f  %12.4f  %12.4f\n", 
//				 RRPD[i], gDOCP, dRRPD[i], 0.0);
				 gDOCP, RRPD[i], 0.0, dRRPD[i]);
// /u/group/hallc/documents/beamline/bcm
	}
	printf("Writing to %s\n",outfilename.Data());
	fclose (outFile);  
//	system(Form("head %s", outfilename.Data()));  
	printf("Writing to %s\n",resfilename.Data());
	fclose (resFile);  
//	system(Form("head %s", resfilename.Data()));  
	printf("Writing to %s\n",qwpresfilename.Data());
	fclose (qwpresFile);  
//	system(Form("head %s", qwpresfilename.Data()));  
	printf("Writing to %s\n",hwpresfilename.Data());
	fclose (hwpresFile);  
	printf("Writing to %s\n",modelfilename.Data());
	fclose (modelFile);  
//	system(Form("head %s", modelfilename.Data()));  
	printf("Writing to %s\n",datacalfilename.Data());
	fclose (datacalFile);  
//	system(Form("head %s", datacalfilename.Data()));  
	printf("Writing to %s\n",modelDOCPfilename.Data());
	fclose (docpFile);  
//	system(Form("head %s", modelDOCPfilename.Data()));  
	printf("Writing to %s\n",modelDOLPfilename.Data());
	fclose (dolpFile);  
//	system(Form("head %s", modelDOLPfilename.Data()));  
	printf("Writing to %s\n",LPanglefilename.Data());
	fclose (angleFile);  
//	system(Form("head %s", LPanglefilename.Data()));  
	printf("Writing to %s\n",paramfilename.Data());
	fclose (paramFile);  
	system(Form("cat %s", paramfilename.Data()));  
	printf("Writing to %s\n", DOCPcorrelfilename.Data());
	fclose (DOCPcorrelFile);  
	//   return;
	printf("Time at the end of job = %f seconds\n",timer.CpuTime());

	FILE *runpointFile, *runpointDOCPFile,  *runpointDOLPFile, *runpointangleFile, *runpointDOLPangleFile, *runpointminmaxangleFile;
	TString runpointfilename=Form("%s.runpoint.out",filebasestring.Data());
	runpointFile = fopen (runpointfilename.Data(),"w");
	runpointfilename=Form("%s.runpointDOCP.out",filebasestring.Data());
	runpointDOCPFile = fopen (runpointfilename.Data(),"w");
	runpointfilename=Form("%s.runpointDOLP.out",filebasestring.Data());
	runpointDOLPFile = fopen (runpointfilename.Data(),"w");
	runpointfilename=Form("%s.runpointangle.out",filebasestring.Data());
	runpointangleFile = fopen (runpointfilename.Data(),"w");
	runpointfilename=Form("%s.runpointDOLPangle.out",filebasestring.Data());
	runpointDOLPangleFile = fopen (runpointfilename.Data(),"w");
	runpointfilename=Form("%s.runpointminmaxangle.out",filebasestring.Data());
	runpointminmaxangleFile = fopen (runpointfilename.Data(),"w");

	printf("opening %s\n",runpointsfile.Data());
	TGraph *runpoints = new TGraph(runpointsfile.Data());
	Int_t markerstyle = 5;
	runpoints->SetMarkerStyle(markerstyle);
	Int_t numrunpoints = runpoints->GetN();
	printf("there are %i runpoints\n",numrunpoints);
	Int_t DOCPsignint;
	for (Int_t i=0; i<numrunpoints; i++) {
		Double_t xp,yp;
		runpoints->GetPoint(i,xp,yp);
		Double_t modelpower = getmodelpower(xp,yp,params,modelnum);

		printf( "runpoint %3i: HWP=%6.2f  QWP=%6.2f  from model:   RRPD=%4.1f  DOCP=%+.4f  DOLP=%.4f  angle=%5.1f degrees\n", 
			   i,xp,yp,modelpower,gDOCP,gDOLP,gangle);
		fprintf(runpointFile,
				"runpoint %3i: HWP=%6.2f  QWP=%6.2f  from model:   RRPD=%4.1f  DOCP=%+.4f  DOLP=%.4f  angle=%5.1f degrees\n",
				i,xp,yp,modelpower,gDOCP,gDOLP,gangle);
		fprintf(runpointDOCPFile,"%.2f  %+.4f\n",yp,gDOCP);
		fprintf(runpointDOLPFile,"%.2f  %+.4f\n",yp,gDOLP);
		fprintf(runpointangleFile, "%.2f  %5.1f\n",yp,gangle);
		fprintf(runpointDOLPangleFile, "%6.3f  %7.2f\n",gDOLP,gangle);
		if (gDOCP<0) DOCPsignint=1; else DOCPsignint=0;
		fprintf(runpointminmaxangleFile, "%i  %7.4f  %7.4f  %7.2f\n",DOCPsignint,1-gDOLP,1+gDOLP,gangle);
	
// 		fprintf(runpointDOCPFile,"%.2f  %.2f  %+.4f\n",xp,yp,gDOCP);
// 		fprintf(runpointDOLPFile,"%.2f  %.2f  %+.4f\n",xp,yp,gDOLP);
// 		fprintf(runpointangleFile, "%.2f  %.2f  %5.1f\n",xp,yp,gangle);
	}
	fclose (runpointFile);  
	fclose (runpointDOCPFile);  
	fclose (runpointDOLPFile);  
	fclose (runpointangleFile);  
	fclose (runpointDOLPangleFile);  
	fclose (runpointminmaxangleFile);  

	// Here maximize the DOCP in a certain range: y is QWP, x is HWP
	Double_t ymin=270, ymax=360, xmin=0, xmax=90;
	printf("\nMinimization of the DOCP\n");
	TF2 *f2 = new TF2("f2",getDOLP_platerr_bir,0.,360.,0.,360.,0);
//	TF2 *f2 = new TF2("f2",getDOLP_plate_bir_pol,0.,360.,0.,360.,0);
	Double_t xminpos=0, yminpos=0;
	f2->SetRange(xmin, ymin, xmax, ymax);
	f2->GetMinimumXY(xminpos, yminpos);
//	Double_t out = f2->Eval(xminpos, yminpos);  // unnecessary, just get from global variable
	printf ("minimum in range (%.1f,%.1f) and (%.1f,%.1f) is:\nHWP = %.3f,  QWP = %.3f, DOLP = %f %%, DOCP = %f %%\n", 
			xmin, xmax, ymin, ymax, xminpos, yminpos, gDOLP*100, gDOCP*100);
	ymin=40; ymax=50; xmin=0; xmax=45;
//	ymin=220; ymax=250; xmin=230; xmax=250;
	f2->SetRange(xmin, ymin, xmax, ymax);
	f2->GetMinimumXY(xminpos, yminpos);
//	out = f2->Eval(xminpos, yminpos);
	printf ("minimum in range (%.1f,%.1f) and (%.1f,%.1f) is:\nHWP = %.3f, QWP = %.3f, DOLP = %f %%, DOCP = %f %%\n", 
			xmin, xmax, ymin, ymax, xminpos, yminpos, gDOLP*100, gDOCP*100);
// 	if (polarizer_angle != 0) {
// 		printf("\nWARNING, polarizer is for vertical i.e. RPD at value %.4f\n",polarizer_angle);
// 	}

	printf("\nusing %s, polarizer_angle=%f\n",detectorname.Data(),polarizer_angle);

// 	TString line = Form("\nnow run:\n.x macros/graph2D_hist.C(\"%s\",\"%s\",\"%s\",\"%s\",0)",
// 						filebasestring.Data(),plotsbasestr.Data(),runpointsfile.Data(),detectorname.Data());
// 	printf("%s\n",line.Data());

// 	line = Form("\nnow run:\ndrawdata(\"%s\",\"%s\",\"%s\")",
// 						plotsbasestr.Data(),runpointsfile.Data(),detectorname.Data());
// 	printf("%s\n",line.Data());

	drawdata(plotsbasestr.Data(),runpointsfile.Data(),detectorname.Data(),params,modelnum);


}

//--------------------------------------------------------------
//--------------------------------------------------------------
// **** The following are drawing routines
//--------------------------------------------------------------
//--------------------------------------------------------------


// void makehistograms() {
// 	printf("making histograms for plots\n");
// 	TH1F *hHWP = new TH1F();

// 	for (Int_t i=0;i<numdatapoints; i++) 
		
// 	}
// }



void drawdata(TString plotnamestr, TString runpointsfile,  TString powname, Double_t *params, Int_t modelnum) {

// 	Int_t xbin = 71;
// 	Double_t xmin = 0;
// 	Double_t xmax = 355;
// 	Int_t ybin = 71;
// 	Double_t ymin = 0;
// 	Double_t ymax = 355;

// 	Int_t xbin = 72;
// 	Double_t xmin = 0;
// 	Double_t xmax = 360;
// 	Int_t ybin = 72;
// 	Double_t ymin = 0;
// 	Double_t ymax = 360;

//	Int_t xbin = 18;
//	Double_t xmin = 0;
//	Double_t xmax = 90;
//	Int_t ybin = 10;
//	Double_t ymin = 40;
//	Double_t ymax = 50;

//	Int_t xbin = 36;
//	Double_t xmin = 0;
//	Double_t xmax = 360;
//	Int_t ybin = 36;
//	Double_t ymin = 0;
//	Double_t ymax = 360;

//	Int_t xbin = 18;
//	Double_t xmin = 0;
//	Double_t xmax = 180;
//	Int_t ybin = 18;
//	Double_t ymin = 0;
//	Double_t ymax = 180;

//	Int_t xbin = 9;
//	Double_t xmin = 0;
//	Double_t xmax = 90;
//	Int_t ybin = 9;
//	Double_t ymin = 0;
//	Double_t ymax = 90;

	Int_t xbin = 24;
	Double_t xmin = 7.5;
	Double_t xmax = 367.5;
//	Int_t ybin = 24;
//	Double_t ymin = 7.5;
//	Double_t ymax = 367.5;

	Int_t ybin = 12;
	Double_t ymin = 187.5;
	Double_t ymax = 367.5;


// 	Int_t xbin = 179;
// 	Double_t xmin = 180.5;
// 	Double_t xmax = 359.5;
// 	Int_t ybin = 180;
// 	Double_t ymin = 179.5;
// 	Double_t ymax = 359.5;

	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetPadRightMargin(0.2);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetTitleSize(0.05,"xyzh");
	gStyle->SetLabelSize(0.05,"xyzh");
	Float_t xoffset=1.5;
	Float_t yoffset=1.2;
	Float_t zoffset=1.4;
	gStyle->SetTitleOffset(xoffset);
	gStyle->SetTitleOffset(yoffset);
	gStyle->SetTitleOffset(zoffset);


	char dirname[255], plotname[255];
	sprintf(dirname,"plots");
	if (plotnamestr=="1") {
		sprintf(plotname,"modeltest_def");
	} else {
		sprintf(plotname,"%s",plotnamestr.Data());
	}

	gStyle->SetPalette(1);

	char filename[255];
	sprintf(filename,"%s",runpointsfile.Data());
	printf("opening %s\n",filename);
	TGraph *runpoints = new TGraph(filename);
	Int_t markerstyle = 5;
	runpoints->SetMarkerStyle(markerstyle);
	printf("There are %i runpoints\n", runpoints->GetN());

	Int_t canvaswidth=800, canvasheight=650;
	char drawstyle2d[255];
	sprintf(drawstyle2d,"colz");


	TH2F *datah = new TH2F("datah",Form("measured data;HWP angle  (degrees);QWP angle  (degrees);%s power  (ADC units)", powname.Data()),
						   xbin,xmin,xmax,ybin,ymin,ymax);
	TH2F *modelh = new TH2F("modelh",Form("fit to data;HWP angle  (degrees);QWP angle  (degrees);%s power  (ADC units)",
										  powname.Data()),xbin,xmin,xmax,ybin,ymin,ymax);
	TH2F *resh = new TH2F("resh",Form("residuals;HWP angle  (degrees);QWP angle  (degrees);%s power  (ADC units)",
									  powname.Data()),xbin,xmin,xmax,ybin,ymin,ymax);
	TH2F *docp = new TH2F("docp","cavity DOCP from model;HWP angle  (degrees);QWP angle  (degrees);DOCP at cavity mirror",
						   xbin,xmin,xmax,ybin,ymin,ymax);
	TH2F *dolp = new TH2F("dolp","cavity DOLP from model;HWP angle  (degrees);QWP angle  (degrees);DOLP at cavity mirror",
						   xbin,xmin,xmax,ybin,ymin,ymax);
	TH2F *angle = new TH2F("angle","cavity LP angle from model;HWP angle  (degrees);QWP angle  (degrees);LP angle  (degrees)",
						   xbin,xmin,xmax,ybin,ymin,ymax);



	Double_t resmin = -95;
	Double_t resmax = 95;
	for (Int_t i=0;i<numdatapoints; i++) {
 		Double_t modelpower = getmodelpower(HWP[i],QWP[i],params,modelnum);
		Double_t res = RRPD[i] - modelpower;
		if (i==1) {
			resmin = resmax = res;
		} else {
			resmin = min(resmin,res);
			resmax = max(resmax,res);
		}
	}
// 	resmin -= 5;
// 	resmax += 5;
	printf("residuals: min=%f  max=%f\n",resmin,resmax);
	TH1F *reshist = new TH1F("reshist",Form("fit residuals;%s power (arb. units)",powname.Data()),100,resmin-0.2,resmax+0.2);
	Double_t docparr[numdatapoints];
	Double_t anglearr[numdatapoints];
	Double_t dolparr[numdatapoints];

	for (Int_t i=0;i<numdatapoints; i++) {
 		Double_t modelpower = getmodelpower(HWP[i],QWP[i],params,modelnum);
		Double_t res = RRPD[i] - modelpower;

		Int_t gbin = datah->FindBin(HWP[i],QWP[i]);
		
		datah->SetBinContent(gbin,RRPD[i]);
		modelh->SetBinContent(gbin,modelpower);
		resh->SetBinContent(gbin,res);
		reshist->Fill(res);

		docp->SetBinContent(gbin,gDOCP);
		dolp->SetBinContent(gbin,gDOLP);
		angle->SetBinContent(gbin,gangle);


		docparr[i]=gDOCP;
		dolparr[i]=gDOLP;
		anglearr[i]=gangle;
		if (i<10) {
			printf("bin %3i  %5.1f %5.1f pow: %5.1f %5.1f %5.1f pol: %f %f %f\n",
				   gbin,HWP[i],QWP[i],RRPD[i],modelpower,res,gDOCP,gDOLP,gangle);
		}
	}

	TGraphErrors *powervsdocp = new TGraphErrors(numdatapoints,docparr,RRPD);
	powervsdocp->SetNameTitle("powervsdocp",Form("%s vs DOCP;DOCP at cavity mirror from model;Measured %s power  (ADC units)",
											   powname.Data(),powname.Data()));
	TGraphErrors *anglevsdocp = new TGraphErrors(numdatapoints,docparr,anglearr);
	anglevsdocp->SetNameTitle("anglevsdocp",Form("LP angle vs DOCP;DOCP at cavity mirror from model;LP angle   (degrees)"));
	TGraphErrors *anglevsdolp = new TGraphErrors(numdatapoints,dolparr,anglearr);
	anglevsdolp->SetNameTitle("anglevsdolp",Form("LP angle vs DOLP;DOLP at cavity mirror from model;LP angle   (degrees)"));
	TGraphErrors *samps = new TGraphErrors(numdatapoints,HWP,QWP);
	samps->SetNameTitle("powervsdocp",Form("Samples;HWP   (degrees);QWP   (degrees)"));

	Double_t modelmin =   modelh->GetMinimum();
	Double_t modelmax =   modelh->GetMaximum();
	Double_t datamin =   datah->GetMinimum();
	Double_t datamax =   datah->GetMaximum();
	Double_t minmin = min(modelmin,datamin);
	Double_t maxmax = max(modelmax,datamax);
	datah->SetMinimum(minmin);
	datah->SetMaximum(maxmax);
	modelh->SetMinimum(minmin);
	modelh->SetMaximum(maxmax);

	TCanvas *canvas2d = new TCanvas("canvas2d","canvas2d",20,0,canvaswidth,canvasheight);
	canvas2d->Divide(2,2,0.001,0.001);
	canvas2d->cd(1); 
	//gStyle->SetPalette(1);
	datah->DrawCopy(drawstyle2d);
	runpoints->Draw("same,p");
	canvas2d->cd(2); 
	modelh->Draw(drawstyle2d);
	runpoints->Draw("same,p");
	canvas2d->cd(3); 
	resh->Draw(drawstyle2d);
	runpoints->Draw("same,p");
	canvas2d->cd(4); 
	reshist->SetStats(1);
	reshist->Draw();
	canvas2d->Print(Form("%s/%s_1powers.png",dirname,plotname));

	TCanvas *canvaspol = new TCanvas("canvaspol","canvaspol",70,0,canvaswidth,canvasheight);
	canvaspol->Divide(2,2,0.001,0.001);
	canvaspol->cd(1); 
	docp->Draw(drawstyle2d);
	runpoints->Draw("same,p");
	canvaspol->cd(2); 
	dolp->Draw(drawstyle2d);
	runpoints->Draw("same,p");
	canvaspol->cd(3); 
	angle->Draw(drawstyle2d);
	runpoints->Draw("same,p");
	canvaspol->cd(4); 
	powervsdocp->Draw("ap");
	canvaspol->Print(Form("%s/%s_2pols.png",dirname,plotname));

	TCanvas *canvascorrel = new TCanvas("canvascorrel","canvascorrel",120,0,canvaswidth,canvasheight);
 	canvascorrel->Divide(2,2,0.001,0.001);
 	canvascorrel->cd(1); 
	samps->Draw("ap");
 	canvascorrel->cd(2); 
 	anglevsdolp->Draw("ap");
	//	canvascorrel->cd(3); 
 	//runpoints->Draw("same,p");
	canvascorrel->cd(3); 
	anglevsdocp->Draw("ap");
	canvascorrel->Print(Form("%s/%s_3samps.png",dirname,plotname));

}


//--------------------------------------------------------------
//--------------------------------------------------------------
// **** The foloowing are general functions
//--------------------------------------------------------------
//--------------------------------------------------------------

void printmodels() {
	printf("model 1: \"platerr_bir\",    9 parameters\n");
	printf("model 2: \"plate_bir_pol\", 11 parameters\n");
	printf("model 3: \"anal\",          10 parameters\n");

}

Int_t getnumparams(Int_t modelnum) {
	Int_t numparams = 0;
	if (modelnum==1) { numparams=9; }
	if (modelnum==2) { numparams=11; }
	if (modelnum==3) { numparams=10; }
	return numparams;
}

Int_t general_setup(TMinuit *gMinuit, Int_t modelnum, Int_t detector, TString &comment) {
	Int_t returnval = 0;
	if (modelnum==1) {
		setup_platerr_bir(gMinuit,detector);
		comment="platerr_bir";
	} else {
		if (modelnum==2) {
			setup_plate_bir_pol(gMinuit,detector);
			comment="plate_bir_pol";
		} else {
			if (modelnum==3) {
				setup_anal(gMinuit,detector);
				comment="anal";
			} else returnval = -1;
		}
	}
	return returnval;
}

Double_t getmodelpower(Double_t xp, Double_t yp, Double_t *params, Int_t modelnum) {
	Double_t modelpower=0;
	if (modelnum==1) modelpower = get_modelpower_platerr_bir(xp,yp,params);
	if (modelnum==2) modelpower = get_modelpower_plate_bir_pol(xp,yp,params);
	if (modelnum==3) modelpower = get_modelpower_anal(xp,yp,params);
	return modelpower;
}


//--------------------------------------------------------------
//--------------------------------------------------------------
// **** The following are the models
//--------------------------------------------------------------
//--------------------------------------------------------------


//---------------------------------------------------------------------
// Model including plate error on the QWP-HWP and a birefringent element
//---------------------------------------------------------------------
Double_t RRPD_platerr_bir(Double_t HWPang, Double_t QWPang, Double_t A,  
						   Double_t gamma,  Double_t B, Double_t HWPfac, Double_t QWPfac, Bool_t moredata) {	
	VJones Vinit = VJones(0,1);
	MJones Mhwp = MJonesRPlate(PI/2*HWPfac,HWPang*rad);
	MJones Mqwp = MJonesRPlate(PI/4*QWPfac,QWPang*rad);
	MJones MEigState = MJonesProd(Mhwp,Mqwp);
	MJones MGeneral = MJonesGenBirefringence(A,gamma,B);
	MJones MEnter = MJonesProd(MGeneral,MEigState);
	MJones MEnterInv = MJonesTranspose(MEnter);
	MJones Mtot = MJonesProd(MEnterInv,MEnter);
	VJones Vfin = MtimesVJones(Mtot,&Vinit);
	VJones Vmid = MtimesVJones(MEnter,&Vinit);
	MJones Mpizer = MJonesPolarizer(polarizer_angle);
	VJones H = MtimesVJones(Mpizer,&Vfin); 
	// let's return all the information we might need from the same function so that the
	// possiblity for errors is reduced.
	if (moredata) { // if required, put some other data in global variables
		gDOCP = Vmid.DOCP();
		gDOLP = Vmid.DOLP();
		gangle = Vmid.Alpha()*deg;
	}
	return H.Norm();
}

Double_t getDOLP_platerr_bir(Double_t *x, Double_t *par) // used to maximize DOCP
{
	RRPD_platerr_bir(x[0]-params[2],x[1]-params[3],params[4],params[5],params[6],params[7],params[8],1);
	Double_t localDOLP = gDOLP;
	return localDOLP;
}

Double_t get_modelpower_platerr_bir(Double_t xp, Double_t yp, Double_t *params) 
{
	Double_t modelpower = params[0] + (RRPD_platerr_bir(xp-params[2],yp-params[3],params[4],params[5],
														params[6],params[7],params[8],1)*params[1]);
	return modelpower;
}

void fcn_platerr_bir(Int_t &npar, Double_t *gin, Double_t &returnval, Double_t *par, Int_t iflag)
{
	//calculate chisquare
	Double_t chisq = 0;
	Double_t delta, delta2;
	for (Int_t i=0;i<numdatapoints; i++) 
	{
		Double_t normpow = RRPD_platerr_bir(HWP[i]-par[2],QWP[i]-par[3],par[4],par[5],par[6],par[7],par[8],0);
		Double_t pow = par[0] + normpow*par[1];
 		delta  = (RRPD[i] - pow)/dRRPD[i]/100;
		delta2 = delta*delta;
		chisq += delta2;
// 		if (numfcncalls==0) {
// 			printf("%4i %4i  %.2f+-%5.2f    %12.2f   %8.2f %8.2f   %12.2f  \n",
// 				   numdatapoints, i, RRPD[i],dRRPD[i], pow, delta,delta2,chisq);
// 		}
	}
	printf("%3i  ",numfcncalls);
	for (Int_t j=0; j<numparams; j++) {
		printf("%9.3f, ",par[j]);
	}
	printf("  %10.4f\n",chisq);
	returnval = chisq;
	numfcncalls++;
	return;
}

void setup_platerr_bir(TMinuit *gMinuit, Int_t detector) {
	gMinuit->SetFCN(fcn_platerr_bir);
	Int_t ierflg = 0;
	char names[maxparams][7] = {"Pedest","Normal","HWPang","QWPang","A",     "gamma", "B", "HWPfac","QWPfac"};

	if (detector==1) // RPD, vertical polarization
	{
		//  Jul 13 scan  
		//	Double_t vstart[maxparams] = {0.0000, 352.5867, 0.0000, -24.2959, 0.3467, -0.1014, 0.0000, 0.9687, 1.0074};
		//	Double_t vstep[maxparams] = {0.0000, 0.7744, 0.0000, 0.0497, 0.0060, 0.0012, 0.0000, 0.0027, 0.0040};
		//  Jul 13 scan: allowing HWP to vary and using analyzer start vals
		// 	Double_t vstart[maxparams] = {0.0000, 353.1638, -12.5504, -49.3996, 0.7847, -0.1013, 0.0000, 0.9642, 1.0073};
		// 	Double_t vstep[maxparams] =  {0.0000, 0.7928, 2.0144, 4.0295, 0.0705, 0.0012, 0.0000, 0.0027, 0.0039};
		// Jul 25: SuperHiRes
		Double_t vstart[] = {0.0000, 0.8, 90, 90.9299, -0.9217, -0.1050, 0.0000, 0.9659, 1.0120};
		Double_t vstep[] = {0.0000, 0.3288, 0.3498, 0.7001, 0.0123, 0.0002, 0.0000, 0.0004, 0.0006};
		for (Int_t i = 0; i<numparams; i++) {
			gMinuit->mnparm(i, names[i], vstart[i], vstep[i],0,0,ierflg);
		}

	} else {
		if (detector==2) // RRPD, horizontal polarization
		{
			// Feb 22 scan
			//	Double_t vstart[maxparams] = {0.0000, 916.5060, 0.0000, -24.2168, 0.3796, -0.1024, 0.0000, 0.9650, 1.0148};
			//	Double_t vstep[maxparams] =  {0.0000, 1.5301, 0.0000, 0.0403, 0.0048, 0.0010, 0.0000, 0.0020, 0.0037};
			//	Double_t vstart[maxparams] = {0.0000, 967.6997, 45.0000, -24.1234, 0.3405, -0.1020, 0.0000, 0.9689, 1.0107};
			//	Double_t vstep[maxparams] =  {0.0000, 0.7949, 0.0000, 0.0196, 0.0024, 0.0005, 0.0000, 0.0010, 0.0016};
			// Jul 13 scan
			// 	Double_t vstart[maxparams] = {0.0000, 967.6997, 45.0000, -24.1234, 0.3405, -0.1020, 0.0000, 0.9689, 1.0107};
			// 	Double_t vstep[maxparams] =  {0.0000, 0.7949, 0.1000, 0.0196, 0.0024, 0.0005, 0.0000, 0.0010, 0.0016};
			// Jul 13 scan: allowing HWP to vary and using analyzer start vals
			// 	Double_t vstart[maxparams] = {0.000,937.554,  0.000,-24.502,  0.564, -0.098,  0.000,  0.963, 1.008};
			// 	Double_t vstep[maxparams] = {0.0000, 4.0507, 0.0000, 0.2385, 0.0305, 0.0058, 0.0000, 0.0152, 0.0254};
			// Jul 25: SuperHiRes
// 			Double_t vstart[] = {0.0000, 2513.4834, 78.7207, 43.3911, -0.8414, -0.1059, 0.0000, 0.9649, 1.0150};  // chi2 21887.3
// 			Double_t vstep[] = {0.0000,  0.3840,  0.1206,  0.2412,  0.0042,  0.0001,  0.0000,  0.0001,  0.0002};
			Double_t vstart[] = {0.03, 0.8, 0.0, 0.0,5.0, 0.0, 0.0000, 1.0, 1.0};   // chi2 18445.4
			Double_t vstep[] = {0.01, 0.1, 0.0004, 0.0001, 0.0004, 0.0001, 0.0000, 0.0001, 0.0002};
			for (Int_t i = 0; i<numparams; i++) {
				gMinuit->mnparm(i, names[i], vstart[i], vstep[i],0,0,ierflg);
			}
		} else {
			if (detector==3) // Analyzer
			{
				
			} 
		}
	}
// 	gMinuit->mnparm(0, "Pedest", vstart[0], vstep[0],0,0,ierflg);
// 	gMinuit->mnparm(1, "Normal", vstart[1], vstep[1],0,0,ierflg);
//  	gMinuit->mnparm(2, "HWPang", vstart[2], vstep[2],0,0,ierflg);
//  	gMinuit->mnparm(3, "QWPang", vstart[3], vstep[3],0,0,ierflg);
//  	gMinuit->mnparm(4, "A",      vstart[4], vstep[4],0,0,ierflg);
//  	gMinuit->mnparm(5, "gamma",  vstart[5], vstep[5],0,0,ierflg);
//   	gMinuit->mnparm(6, "B",      vstart[6], vstep[6],0,0,ierflg);
//  	gMinuit->mnparm(7, "HWPfac", vstart[7], vstep[7],0,0,ierflg);
//    	gMinuit->mnparm(8, "QWPfac", vstart[8], vstep[8],0,0,ierflg);
}





//---------------------------------------------------------------------------
// Include plate error on the QWP-HWP, a birefringent element and a polarizer
//---------------------------------------------------------------------------
Double_t RRPD_plate_bir_pol(Double_t HWPang, Double_t QWPang, Double_t A, Double_t gamma,  Double_t B, 
							Double_t HWPfac, Double_t QWPfac, Double_t Pol, Double_t Theta, Bool_t moredata) {	
	VJones Vinit = VJones(0,1);
	MJones Mhwp = MJonesRPlate(PI/2*HWPfac,HWPang*rad);
	MJones Mqwp = MJonesRPlate(PI/4*QWPfac,QWPang*rad);
	MJones MEigState = MJonesProd(Mhwp,Mqwp);
	MJones MGeneral = MJonesGenBirefringence(A,gamma,B);
	MJones MPartPol = MJonesRPartPol(Pol, Theta);
	MGeneral = MJonesProd(MPartPol,MGeneral);
	MJones MEnter = MJonesProd(MGeneral,MEigState);
	MJones MEnterInv = MJonesTranspose(MEnter);
	MJones Mtot = MJonesProd(MEnterInv,MEnter);
	VJones Vfin = MtimesVJones(Mtot,&Vinit);
	VJones Vmid = MtimesVJones(MEnter,&Vinit);
	MJones Mpizer = MJonesPolarizer(polarizer_angle);
//	VJones H = Vmid;                        // this is for analyzed light
	VJones H = MtimesVJones(Mpizer,&Vfin); // this is for reflected light
	// return all the information we might need from the same function so that the
	// possiblity for errors is reduced.
	if (moredata) { // if required, put some other data in global variables
		gDOCP = Vmid.DOCP();
		gDOLP = Vmid.DOLP();
		gangle = Vmid.Alpha()*deg;
	}
	return H.Norm();
}

Double_t getDOLP_plate_bir_pol(Double_t *x, Double_t *par) // used to maximize DOCP
{
	RRPD_plate_bir_pol(x[0]-params[2],x[1]-params[3],params[4],params[5],params[6],
					   params[7],params[8],params[9],params[10],1);
	Double_t localDOLP = gDOLP;
	return localDOLP;
}

void fcn_plate_bir_pol(Int_t &npar, Double_t *gin, Double_t &returnval, Double_t *par, Int_t iflag)
{
	//calculate chisquare
	Double_t chisq = 0;
	Double_t delta, delta2;
	for (Int_t i=0;i<numdatapoints; i++) 
	{
		Double_t normpow = RRPD_plate_bir_pol(HWP[i]-par[2],QWP[i]-par[3],par[4],par[5],par[6],
											  par[7],par[8],par[9],par[10],0);
		Double_t pow = par[0] + normpow*par[1];
 		delta  = (RRPD[i] - pow)/dRRPD[i]/100;
		delta2 = delta*delta;
		chisq += delta2;
// 		if (numfcncalls==0) {
// 			printf("%4i %4i  %.2f+-%5.2f    %12.2f   %8.2f %8.2f   %12.2f  \n",
// 				   numdatapoints, i, RRPD[i],dRRPD[i], pow, delta,delta2,chisq);
// 		}
	}
	printf("%3i  ",numfcncalls);
	for (Int_t j=0; j<numparams; j++) {
		printf("%8.3f, ",par[j]);
	}
	printf("  %.3f\n",chisq);
	returnval = chisq;
	numfcncalls++;
	return;
}

Double_t get_modelpower_plate_bir_pol(Double_t xp, Double_t yp, Double_t *params) 
{
	Double_t modelpower = params[0] + (RRPD_plate_bir_pol(xp-params[2],yp-params[3],params[4],params[5],
														  params[6],params[7],params[8],params[9],params[10],1)*params[1]);
	return modelpower;
}

void setup_plate_bir_pol(TMinuit *gMinuit, Int_t detector) {
	Int_t ierflg = 0;

	gMinuit->SetFCN(fcn_plate_bir_pol);

	Double_t vstart[maxparams] = {0.0000, 17739.6885, 0.0000, 0, 2.3819, -0.1718, 0.0000, 1.0259, 1, 1.0000,0};
	Double_t vstep[maxparams] =  {0.0000, 15.0352, 0.0000, 0.0288, 0.0004, 0.0008, 0.0000, 0.0004, 0.0009, 0.0005, 0.0004};

	gMinuit->mnparm(0, "Pedest", vstart[0], vstep[0],0,0,ierflg);
	gMinuit->mnparm(1, "Normal", vstart[1], vstep[1],0,0,ierflg);
 	gMinuit->mnparm(2, "HWPang", vstart[2], vstep[2],0,0,ierflg);
 	gMinuit->mnparm(3, "QWPang", vstart[3], vstep[3],0,0,ierflg);
 	gMinuit->mnparm(4, "A",      vstart[4], vstep[4],0,0,ierflg);
 	gMinuit->mnparm(5, "gamma",  vstart[5], vstep[5],0,0,ierflg);
  	gMinuit->mnparm(6, "B",      vstart[6], vstep[6],0,0,ierflg);
 	gMinuit->mnparm(7, "HWPfac", vstart[7], vstep[7],0,0,ierflg);
   	gMinuit->mnparm(8, "QWPfac", vstart[8], vstep[8],0,0,ierflg);
 	gMinuit->mnparm(9, "Pol",    vstart[9], vstep[9],0,0,ierflg);
   	gMinuit->mnparm(10, "Theta", vstart[10], vstep[10],0,0,ierflg);
}



//-----------------------------------------------------------------------------
// Analyzer model: including plate error on the QWP-HWP, a birefringent element 
//-----------------------------------------------------------------------------
Double_t RRPD_anal(Double_t HWPang, Double_t QWPang, Double_t A, Double_t gamma,  Double_t B, 
				   Double_t HWPfac, Double_t QWPfac, Double_t analyzer_angle, Bool_t moredata) {	
	VJones Vinit = VJones(1,0);
	MJones Mhwp = MJonesRPlate(PI/2*HWPfac,HWPang*rad);
	MJones Mqwp = MJonesRPlate(PI/4*QWPfac,QWPang*rad);
	MJones MEigState = MJonesProd(Mqwp,Mhwp);
	MJones MGeneral = MJonesGenBirefringence(A,gamma,B);
	MJones MEnter = MJonesProd(MGeneral,MEigState);
	VJones Vmid = MtimesVJones(MEnter,&Vinit);
	MJones Mpizer = MJonesPolarizer(analyzer_angle);
//	MJones Manal = MJonesProd(Mpizer,MEnter);
//	VJones H = Vmid;                        // this is for analyzed light
	VJones H = MtimesVJones(Mpizer,&Vmid); // this is for reflected light
	// return all the information we might need from the same function so that the
	// possiblity for errors is reduced.
	if (moredata) { // if required, put some other data in global variables
		gDOCP = Vmid.DOCP();
		gDOLP = Vmid.DOLP();
		gangle = Vmid.Alpha()*deg;
	}
	return H.Norm();
}

Double_t getDOLP_anal(Double_t *x, Double_t *par) // used to maximize DOCP
{
	RRPD_anal(x[0]-params[2],x[1]-params[3],params[4],params[5],params[6],
			  params[7],params[8],params[9],1);
	Double_t localDOLP = gDOLP;
	return localDOLP;
}

void fcn_anal(Int_t &npar, Double_t *gin, Double_t &returnval, Double_t *par, Int_t iflag)
{
	//calculate chisquare
	Double_t chisq = 0;
	Double_t delta, delta2;
	for (Int_t i=0;i<numdatapoints; i++) 
	{
		Double_t normpow = RRPD_anal(HWP[i]-par[2],QWP[i]-par[3],par[4],par[5],par[6],
											  par[7],par[8],par[9],0);
		Double_t pow = par[0] + normpow*par[1];
 		delta  = (RRPD[i] - pow)/dRRPD[i]/100;
		delta2 = delta*delta;
		chisq += delta2;
// 		if (numfcncalls==0) {
// 			printf("%4i %4i  %.2f+-%5.2f    %12.2f   %8.2f %8.2f   %12.2f  \n",
// 				   numdatapoints, i, RRPD[i],dRRPD[i], pow, delta,delta2,chisq);
// 		}
	}
	printf("%3i  ",numfcncalls);
	for (Int_t j=0; j<numparams; j++) {
		printf("%8.3f, ",par[j]);
	}
	printf("  %.3f\n",chisq);
	returnval = chisq;
	numfcncalls++;
	return;
}

Double_t get_modelpower_anal(Double_t xp, Double_t yp, Double_t *params) 
{
	Double_t modelpower = params[0] + (RRPD_anal(xp-params[2],yp-params[3],params[4],params[5],
														  params[6],params[7],params[8],params[9],1)*params[1]);
	return modelpower;
}

void setup_anal(TMinuit *gMinuit, Int_t detector) {
	Int_t ierflg = 0;

	gMinuit->SetFCN(fcn_anal);

//	Double_t vstart[maxparams] = {0.0000, 17739.6, 0.0000, -24.0383, 10.1010, -43.6760, 0.0000, 1.0259, 2.9918, -0.1574};
//	Double_t vstep[maxparams] =  {0.0000, 15.2016, 1.0000, 0.0295, 0.0072, 0.0127, 0.0000, 0.0004, 0.0009, 0.0059};

	Double_t vstart[maxparams] = {0.0000, 17768.9541, -12.5504, -49.3996, 0.7847, -0.1013, 0.0000, 0.9642, 1.0073, -0.089};
	Double_t vstep[maxparams] =  {0.1000, 15.3219, 0.3527, 0.7059, 0.0086, 0.0208, 0.1000, 0.0004, 0.0009, 0.0004};

	gMinuit->mnparm(0, "Pedest", vstart[0], vstep[0],0,0,ierflg);
	gMinuit->mnparm(1, "Normal", vstart[1], vstep[1],0,0,ierflg);
 	gMinuit->mnparm(2, "HWPang", vstart[2], vstep[2],0,0,ierflg);
 	gMinuit->mnparm(3, "QWPang", vstart[3], vstep[3],0,0,ierflg);
 	gMinuit->mnparm(4, "A",      vstart[4], vstep[4],0,0,ierflg);
 	gMinuit->mnparm(5, "gamma",  vstart[5], vstep[5],0,0,ierflg);
  	gMinuit->mnparm(6, "B",      vstart[6], vstep[6],0,0,ierflg);
 	gMinuit->mnparm(7, "HWPfac", vstart[7], vstep[7],0,0,ierflg);
   	gMinuit->mnparm(8, "QWPfac", vstart[8], vstep[8],0,0,ierflg);
   	gMinuit->mnparm(9, "polang", vstart[9], vstep[9],0,0,ierflg);
}







/* emacs
 * Local Variables:
 * mode:C++
 * mode:font-lock
 * c-file-style: "stroustrup"
 * tab-width: 4
 * End:
 */

 
