


void compare_method() {
	gROOT->Reset();

	TCanvas *canvasall;
	if (gROOT->IsBatch()) {
		canvasall = new TCanvas("canvasall","Summary",40,80,1200,1400);
	} else {
		canvasall = new TCanvas("canvasall","Summary",40,80,1000,1000);
	}
	gROOT->SetStyle("Plain");  
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetLabelSize(0.05,"hxyz");
	gStyle->SetTitleSize(0.05,"hxyz");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadRightMargin(0.12);
	gStyle->SetOptFit(1);

	TTree *etree = new TTree("etree","entrance tree");
	char etreefile[1000];
	sprintf(etreefile,"/home/cdaq/users/dalton/compton/entrancefunction/output/TFscanopen_checkmethod_constsub_platerr_bir/platerr_bir.out");
	printf("opening %s\n",etreefile);
	etree->ReadFile(etreefile,"eHWP:eQWP:eRRPD:edRRPD:emodelpower:eres:eDOCP:eDOLP:eangle");
	etree->SetMarkerStyle(21);
   

	TTree *atree = new TTree("atree","actual tree");
	char atreefile[1000];
	sprintf(atreefile,"/home/cdaq/users/dalton/compton/TransferFunction/data/output/TFscanopen_checkmethod/TFscanopen_checkmethod_outputforfit.txt");
	printf("opening %s\n",atreefile);

	atree->ReadFile(atreefile,"HWP:QWP:DOCP:dDOCP:DOLP:dDOLP:angle:dangle:min:max:mean:resmean:resrms");
	atree->SetMarkerStyle(20);
	atree->SetMarkerSize(0.5);
	atree->AddFriend(etree);

	atree->SetAlias("goodangle","180+25.8982-angle-180*(eangle<40&&angle<40)");

	Double_t marginval = 0.14;
	canvasall->Divide(3,4,0.001,0.001);
	Int_t padnum=0;

//	gStyle->SetOptStat(0);
	gStyle->SetStatW(0.5);
//	gStyle->SetStatW(0.4);


	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("eDOCP:DOCP>>hist%i",padnum),"");
	hist1->SetTitle("DOCP;Direct measurement;Entrance fit");
	hist1->Draw();
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("eDOLP:DOLP>>hist%i",padnum),"");
	hist2->SetTitle("DOLP;Direct measurement;Entrance fit");
	hist2->Draw();
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("eangle:goodangle>>hist%i",padnum),"");
	hist3->SetTitle("angle;Direct measurement;Entrance fit");
	hist3->Draw();


//	gStyle->SetOptStat(2200);
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("abs(eDOCP)-DOCP>>hist%i",padnum),"","e");
	hist4->SetTitle("DOCP difference;Absolute DOCP difference");
	hist4->SetStats();
// 	gaus->SetLineWidth(1);
// 	gaus->SetLineColor(kRed);
// 	hist4->Fit("gaus");
// 	gaus->SetLineWidth(1);
// 	gaus->SetLineColor(kRed);
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("eDOLP-DOLP>>hist%i",padnum),"","e");
	hist5->SetTitle("DOLP difference;Absolute DOLP difference");
	hist5->SetStats();
// 	hist5->Fit("gaus");
// 	gaus->SetLineWidth(1);
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("eangle-goodangle>>hist%i",padnum),"","e");
	hist6->SetTitle("angle difference;Absolute angle difference  (degrees)");
	hist6->SetStats();
// 	hist6->Fit("gaus");
// 	gaus->SetLineWidth(1);


//	gStyle->SetOptStat(0);
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("abs(eDOCP)-DOCP:DOCP>>hist%i",padnum),"");
	hist7->SetTitle("DOCP difference;DOCP  (direct);Absolute DOCP difference");
	hist7->Draw();
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("eDOLP-DOLP:DOLP>>hist%i",padnum),"");
	hist8->SetTitle("DOLP difference;DOLP  (direct);Absolute DOLP difference");
	hist8->Draw();
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("eangle-goodangle:DOLP>>hist%i",padnum),"");
	hist9->SetTitle("angle difference;DOLP  (direct);Absolute angle difference  (degrees)");
	hist9->Draw();
	canvasall->Update();


	atree->SetMarkerStyle(21);
	atree->SetMarkerSize(1);
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("QWP:HWP:abs(eDOCP)-DOCP>>hist%i",padnum),"","colz");
	hist10->SetTitle("DOCP difference;HWP;QWP");
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("QWP:HWP:eDOLP-DOLP>>hist%i",padnum),"","colz");
	hist11->SetTitle("DOLP difference;HWP;QWP");
	padnum++; canvasall->cd(padnum);
	gPad->SetRightMargin(marginval);
	atree->Draw(Form("QWP:HWP:eangle-goodangle>>hist%i",padnum),"","colz");
	hist12->SetTitle("angle difference;HWP;QWP");

	canvasall->Update();
	if (gROOT->IsBatch()) {
		canvasall->Print("plots/compare_method_big.png");
	} else {
		canvasall->Print("plots/compare_method.png");
	}



}

/* emacs
 * Local Variables:
 * mode:C++
 * mode:font-lock
 * c-file-style: "stroustrup"
 * tab-width: 4
 * End:
 */
