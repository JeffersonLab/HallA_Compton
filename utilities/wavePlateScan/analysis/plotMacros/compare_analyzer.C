

void compare_analyzer() {

	gROOT->SetStyle("Plain");
	gStyle->SetPadRightMargin(0.04);
	gStyle->SetPadLeftMargin(0.12);
	//  gStyle->SetTitleOffset(0.9,"x");
	gStyle->SetTitleOffset(1.5,"y");
//   gStyle->SetLabelSize(0.055,"hxyz");
//   gStyle->SetTitleSize(0.055,"hxyz");

	char plotname[255];
	sprintf(plotname,"plots/compare_analyzer.png");

	TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
	TGraph *docp_QWP_fit = new TGraph("output/TFscanexit_Analyzer_200mW_constsub_analallparam.runpointDOCP.out","%lg %lg");
	TGraphErrors *docp_QWP_meas1 = new TGraphErrors("/home/cdaq/users/dalton/compton/TransferFunction/data/output/TFscanpol_THWP18000_200mW/TFscanpol_THWP18000_200mW_DOCP.txt", "%lg %lg %lg");
	TGraphErrors *docp_QWP_meas2 = new TGraphErrors("/home/cdaq/users/dalton/compton/TransferFunction/data/output/TFscanpol_THWP19992_200mW/TFscanpol_THWP19992_200mW_DOCP.txt", "%lg %lg %lg");
	TGraphErrors *docp_QWP_meas3 = new TGraphErrors("/home/cdaq/users/dalton/compton/TransferFunction/data/output/TFscanpol_THWP22002_200mW/TFscanpol_THWP22002_200mW_DOCP.txt", "%lg %lg %lg");
	canvas->cd();
	docp_QWP_fit->SetMarkerStyle(20);
	docp_QWP_fit->SetMarkerColor(kBlue);
	docp_QWP_meas1->SetMarkerStyle(21);
	docp_QWP_meas1->SetMarkerColor(kGreen);
	docp_QWP_meas2->SetMarkerStyle(22);
	docp_QWP_meas2->SetMarkerColor(kGreen);
	docp_QWP_meas3->SetMarkerStyle(23);
	docp_QWP_meas3->SetMarkerColor(kGreen);
	TMultiGraph *mg_docp = new TMultiGraph();
	mg_docp->Add(docp_QWP_fit,"p");
	mg_docp->Add(docp_QWP_meas1,"p");
 	mg_docp->Add(docp_QWP_meas2,"p");
 	mg_docp->Add(docp_QWP_meas3,"p");
	mg_docp->SetTitle("State 1: DOCP in exit line;QWP angle;DOCP");
	mg_docp->Draw("a");

	leg3 = new TLegend(0.7,0.85,0.99,0.99);
	//   leg->SetHeader("The Legend Title");
	leg3->AddEntry(docp_QWP_fit,"fit results","p");
	leg3->AddEntry(docp_QWP_meas1,"direct measurement","p");
	//leg3->SetBorderSize(0);
	leg3->Draw();
	
	canvas->Print("plots/analyzer_check.png");


	TCanvas *canvas_angle = new TCanvas("canvas_angle","canvas_angle",800,600);
	TGraph *angle_QWP_fit = new TGraph("output/TFscanexit_Analyzer_200mW_constsub_analallparam.runpointangle.out","%lg %lg");
	TGraphErrors *angle_QWP_meas1 = new TGraphErrors("/home/cdaq/users/dalton/compton/TransferFunction/data/output/TFscanpol_THWP18000_200mW/TFscanpol_THWP18000_200mW_angle.txt", "%lg %lg %lg");
	TGraphErrors *angle_QWP_meas2 = new TGraphErrors("/home/cdaq/users/dalton/compton/TransferFunction/data/output/TFscanpol_THWP19992_200mW/TFscanpol_THWP19992_200mW_angle.txt", "%lg %lg %lg");
	TGraphErrors *angle_QWP_meas3 = new TGraphErrors("/home/cdaq/users/dalton/compton/TransferFunction/data/output/TFscanpol_THWP22002_200mW/TFscanpol_THWP22002_200mW_angle.txt", "%lg %lg %lg");
	canvas_angle->cd();
	angle_QWP_fit->SetMarkerStyle(20);
	angle_QWP_fit->SetMarkerColor(kBlue);
	angle_QWP_meas1->SetMarkerStyle(21);
	angle_QWP_meas1->SetMarkerColor(kGreen);
	angle_QWP_meas2->SetMarkerStyle(22);
	angle_QWP_meas2->SetMarkerColor(kGreen);
	angle_QWP_meas3->SetMarkerStyle(23);
	angle_QWP_meas3->SetMarkerColor(kGreen);
	TMultiGraph *mg_angle = new TMultiGraph();
	mg_angle->Add(angle_QWP_fit,"p");
	mg_angle->Add(angle_QWP_meas1,"p");
 	mg_angle->Add(angle_QWP_meas2,"p");
 	mg_angle->Add(angle_QWP_meas3,"p");
	mg_angle->SetTitle("State 1: angle in exit line;QWP angle;angle");
	mg_angle->Draw("a");

	leg3 = new TLegend(0.7,0.85,0.99,0.99);
	//   leg->SetHeader("The Legend Title");
	leg3->AddEntry(docp_QWP_fit,"fit results","p");
	leg3->AddEntry(docp_QWP_meas1,"direct measurement","p");
	//leg3->SetBorderSize(0);
	leg3->Draw();
	canvas_angle->Print("plots/analyzer_check.png");

}




/* emacs
 * Local Variables:
 * mode:C++
 * mode:font-lock
 * c-file-style: "stroustrup"
 * tab-width: 4
 * End:
 */

 
