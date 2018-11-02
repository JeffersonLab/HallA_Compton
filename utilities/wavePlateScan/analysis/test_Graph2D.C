

#include <algorithm>


void test_Graph2D(TString datanamestr="1", TString plotnamestr="1", 
		  TString runpointsfile = "data/runpoints.txt", 
		  TString powname = "RRPD", Bool_t plotbig=0)
{
  if (datanamestr=="1") {
    printf("\n\ttest_Graph2D(datanamestr,[plotnamestr]\n\n");
    return;
  }

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.05,"xyzh");
  gStyle->SetLabelSize(0.05,"xyzh");

  char dirname[255], plotname[255];
  sprintf(dirname,"plots");
  if (plotnamestr=="1") {
    sprintf(plotname,"modeltest_def");
  } else {
    sprintf(plotname,"%s",plotnamestr.Data());
  }


//   const UInt_t NumStops = 3;
//   const UInt_t PalDepth = 20;
//   Int_t CyclePalette[PalDepth];
//   Double_t Red[NumStops]   = { 0.0, 1.0, 0.0 };
//   Double_t Green[NumStops] = { 0.0, 0.0, 0.0 };
//   Double_t Blue[NumStops]  = { 1.0, 0.0, 1.0 };
//   Double_t Stops[NumStops] = { 0.0, 0.5, 1.0 };

//   const UInt_t NumStops = 5;
//   const UInt_t PalDepth = 40;
//   Int_t CyclePalette[PalDepth];
//   Double_t Red[NumStops]   = { 0.0, 0.5, 1.0, 0.5, 0.0 };
//   Double_t Green[NumStops] = { 0.5, 0.0, 0.5, 1.0, 0.5 };
//   Double_t Blue[NumStops]  = { 1.0, 0.5, 0.0, 0.5, 1.0 };
//   Double_t Stops[NumStops] = { 0.0, 0.25, 0.5, 0.75, 1.0 };
//   Int_t FI = TColor::CreateGradientColorTable(NumStops, Stops, Red, Green, Blue, PalDepth);
//   for (int i=0;i<PalDepth;i++) CyclePalette[i] = FI+i;

//   gStyle->SetPalette(PalDepth, CyclePalette);

  char filename[255];
  sprintf(filename,"%s",runpointsfile.Data());
  printf("opening %s\n",filename);
  TGraph *runpoints = new TGraph(filename);
  runpoints->SetMarkerStyle(3);
  printf("There are %i runpoints\n", runpoints->GetN());

  sprintf(filename,"output/%s.datacal.out",datanamestr.Data());
  printf("opening %s\n",filename);
  TGraph2D *data = new TGraph2D(filename);
  data->SetNameTitle("data",Form("measured data;HWP angle  (degrees);QWP angle  (degrees);%s power  (ADC units)",powname.Data()));
  Int_t nbins = data->GetN();
  Int_t rootnbins = sqrt(nbins);
  printf("we have %i bins, %i\n",nbins,rootnbins);
  data->SetNpx(rootnbins);
  data->SetNpy(rootnbins);

  Float_t xoffset=1.5;
  Float_t yoffset=1.2;
  Float_t zoffset=1.4;

  data->GetXaxis()->SetTitleOffset(xoffset);
  data->GetYaxis()->SetTitleOffset(yoffset);
  data->GetZaxis()->SetTitleOffset(zoffset);
  Double_t datamin =   data->GetZmin();
  Double_t datamax =   data->GetZmax();

  sprintf(filename,"output/%s.model.out",datanamestr.Data());
  printf("opening %s\n",filename);
  TGraph2D *model = new TGraph2D(filename);
  model->SetNameTitle("model",Form("fit to data;HWP angle  (degrees);QWP angle  (degrees);%s power  (ADC units)",powname.Data()));
  model->GetXaxis()->SetTitleOffset(xoffset);
  model->GetYaxis()->SetTitleOffset(yoffset);
  model->GetZaxis()->SetTitleOffset(zoffset);
  Double_t modelmin =   model->GetZmin();
  Double_t modelmax =   model->GetZmax();
  Double_t minmin = min(modelmin,datamin);
  Double_t maxmax = max(modelmax,datamax);

  model->GetZaxis()->SetRangeUser(minmin,maxmax);
  data->GetZaxis()->SetRangeUser(minmin,maxmax);

  printf("min and max:  data: %f  %f   model: %f  %f   both: %f  %f\n", 
	 datamin, datamax, modelmin, modelmax, minmin, maxmax);

  sprintf(filename,"output/%s.docp.out",datanamestr.Data());
  printf("opening %s\n",filename);
  TGraph2D *docp = new TGraph2D(filename);
  docp->SetNameTitle("docp","cavity DOCP from model;HWP angle  (degrees);QWP angle  (degrees);DOCP at cavity mirror");
  docp->GetXaxis()->SetTitleOffset(xoffset);
  docp->GetYaxis()->SetTitleOffset(yoffset);
  docp->GetZaxis()->SetTitleOffset(zoffset);

  sprintf(filename,"output/%s.dolp.out",datanamestr.Data());
  printf("opening %s\n",filename);
  TGraph2D *dolp = new TGraph2D(filename);
  dolp->SetNameTitle("dolp","cavity DOLP from model;HWP angle  (degrees);QWP angle  (degrees);DOLP at cavity mirror");
  dolp->GetXaxis()->SetTitleOffset(xoffset);
  dolp->GetYaxis()->SetTitleOffset(yoffset);
  dolp->GetZaxis()->SetTitleOffset(zoffset);


  sprintf(filename,"output/%s.angle.out",datanamestr.Data());
  printf("opening %s\n",filename);
  TGraph2D *angle = new TGraph2D(filename);
  angle->SetNameTitle("angle","cavity LP angle from model;HWP angle  (degrees);QWP angle  (degrees);LP angle  (degrees)");
  angle->GetXaxis()->SetTitleOffset(xoffset);
  angle->GetYaxis()->SetTitleOffset(yoffset);
  angle->GetZaxis()->SetTitleOffset(zoffset);

  sprintf(filename,"output/%s.residuals.out",datanamestr.Data());
  printf("opening %s\n",filename);
  TGraph2D *res = new TGraph2D(filename);
  res->SetNameTitle("res",Form("fit residuals;HWP angle  (degrees);QWP angle  (degrees);%s power  (ADC units)",powname.Data()));
  res->GetXaxis()->SetTitleOffset(xoffset);
  res->GetYaxis()->SetTitleOffset(yoffset);
  res->GetZaxis()->SetTitleOffset(zoffset);

  sprintf(filename,"output/%s.qwpresiduals.out",datanamestr.Data());
  printf("opening %s\n",filename);
  TGraph *qwpres = new TGraph(filename);
  Double_t resmin = -95;
  Double_t resmax = 95;
  printf("all data loaded\n");

  // Find the real min and max (this takes too long) 
  printf("finding residuals min and max\n");
  for (Int_t i=1; i<=qwpres->GetN(); i++) {
    //printf("%i\n",i);
    Double_t xp,yp;
    qwpres->GetPoint(i,xp,yp);
    if (i==1) {
      resmin = resmax = yp;
    } else {
      resmin = min(resmin,yp);
      resmax = max(resmax,yp);
    }
  }
  resmin -= 5;
  resmax += 5;

  printf("residuals: min=%f  max=%f\n",resmin,resmax);
  TH1F *reshist = new TH1F("reshist",Form("fit residuals;%s power (ADC units)",powname.Data()),200,resmin-5,resmax+5);
  for (Int_t i=1; i<=qwpres->GetN(); i++) {
    Double_t xp,yp;
    qwpres->GetPoint(i,xp,yp);
    reshist->Fill(yp);
  }

  Int_t canvaswidth=800, canvasheight=650;
  if (plotbig) {
    canvaswidth=1000, canvasheight=800;
  }

  char drawstyle[255];
  sprintf(drawstyle,"colz");

  TCanvas *canvas2d = new TCanvas("canvas2d","canvas2d",20,0,canvaswidth,canvasheight);
  canvas2d->Divide(2,2,0.001,0.001);
  canvas2d->cd(1); 
  data->Draw(drawstyle);
  runpoints->Draw("same,p");
  canvas2d->cd(2); 
  model->Draw(drawstyle);
  runpoints->Draw("same,p");
  canvas2d->cd(3); 
  res->Draw(drawstyle);
  runpoints->Draw("same,p");
  canvas2d->cd(4); 
  reshist->Draw();
  canvas2d->Print(Form("%s/%s.png",dirname,plotname));




  TCanvas *canvaspol = new TCanvas("canvaspol","canvaspol",40,0,canvaswidth,canvasheight);
  canvaspol->Divide(2,2,0.001,0.001);
  canvaspol->cd(1); 
  docp->Draw(drawstyle);
  runpoints->Draw("same,p");
  canvaspol->cd(2); 
  dolp->Draw(drawstyle);
  runpoints->Draw("same,p");
  canvaspol->cd(3); 
  angle->Draw(drawstyle);
  runpoints->Draw("same,p");

  sprintf(filename,"output/%s.docpvsrrpd.out",datanamestr.Data());
  printf("opening %s\n",filename);
  TGraphErrors *docpvsrrpd = new TGraphErrors(filename,"%lg %lg %lg %lg");
  //docpvsrrpd->SetMarkerStyle(3);
  docpvsrrpd->SetNameTitle("docpvsrrpd",Form("DOCP vs %s;DOCP at cavity mirror from model;Measured %s power  (ADC units)",
					     powname.Data(),powname.Data()));
  //docpvsrrpd->GetYaxis()->SetTitleOffset(yoffset);
  //docpvsrrpd->GetZaxis()->SetTitleOffset(zoffset);
  canvaspol->cd(4); 
  docpvsrrpd->Draw("ap");

  canvaspol->Print(Form("%s/%s_pol.png",dirname,plotname));

  // * It turns out that this method is too inaccurate and can't be used.
  // * the actual values must be determined from the fit itself
//   TH2D *docpHist =  docp->GetHistogram();
//   TH2D *dolpHist =  dolp->GetHistogram();
//   TH2D *angleHist =  angle->GetHistogram();
//   printf("DOCP: min: %.4f  max: %.4f\n",docpHist->GetMinimum(),docpHist->GetMaximum());
//   for (Int_t i=0; i<runpoints->GetN(); i++) {
//     Double_t xp,yp;
//     runpoints->GetPoint(i,xp,yp);
//     Double_t docpout = docpHist->GetBinContent(docpHist->FindBin(xp,yp));
//     Double_t dolpout = dolpHist->GetBinContent(dolpHist->FindBin(xp,yp));
//     Double_t angleout = angleHist->GetBinContent(angleHist->FindBin(xp,yp));
//     printf("runpoint %i: HWP=%.1f  QWP=%.1f  DOCP=%+.4f  DOCP=%.4f  angle=%.1f\n", 
// 	   i,xp,yp,docpout,dolpout,angleout);
//   }

  TCanvas *testcanvas = new TCanvas("testcanvas","testcanvas",60,0,canvaswidth,canvasheight);
  TH2D *docpHist =  docp->GetHistogram();
  docpHist->Draw();
//   TH2D *dolpHist =  dolp->GetHistogram();
//   TH2D *angleHist =  angle->GetHistogram();



  // Set up for the surf plots
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.2);
//   gStyle->SetPadBottomMargin(0.15);
//   Double_t theta = 30;
//   Double_t phi = 25;
  Double_t theta = 0;
  //  Double_t phi = 39; // this is for the smaller 
  Double_t phi = 26;
//   TCanvas *canvas2dsurf = new TCanvas("canvas2dsurf","canvas2dsurf",40,0,800,680);
//   canvas2dsurf->Divide(2,2,0.001,0.001);
//   canvas2dsurf->cd(1); 
//   gPad->SetTheta(theta);
//   gPad->SetPhi(phi);
//   data->Draw("surf1");
//   canvas2dsurf->cd(2); 
//   gPad->SetTheta(theta);
//   gPad->SetPhi(phi);
//   model->Draw("surf1");
//   canvas2dsurf->cd(3); 
//   gPad->SetTheta(theta);
//   gPad->SetPhi(phi);
//   res->Draw("surf1");
//   canvas2dsurf->Print(Form("%s/%s_surf.png",dirname,plotname));

}

