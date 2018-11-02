

plot_runpoints ()
{
	TGraph *anglesgr = new TGraph("output/TFscanHiResEntrance_Jul13_2012_RRPD_constsub_platerr_birhwp10/platerr_birhwp10.runpointDOLPangle.out","%lg %lg");
	anglesgr->SetMarkerStyle(20);
	anglesgr->SetMarkerColor(kBlue);
	//  anglesgr->Draw("ap");

	TGraph *anglesgr2 = new TGraph("output/TFscanHiResEntrance_Jul13_2012_RRPD_constsub_platerr_bir/platerr_bir.runpointDOLPangle.out","%lg %lg");
	anglesgr2->SetMarkerStyle(20);
	anglesgr2->SetMarkerColor(kRed);

	TGraph *anglesgr3 = new TGraph("output/TFscanHiResEntrance_Jul13_2012_RPD_constsub_platerr_bir/platerr_bir.runpointDOLPangle.out","%lg %lg");
	anglesgr3->SetMarkerStyle(20);
	anglesgr3->SetMarkerColor(kGreen);

	TGraph *anglesgr4 = new TGraph("output/TFscanHiResEntrance_Jul13_2012_RRPD_constsub_platerr_birhwp45var/platerr_birhwp45var.runpointDOLPangle.out",
								   "%lg %lg");
	anglesgr4->SetMarkerStyle(21);
//	anglesgr4->SetMarkerColor(kGreen);

	TMultiGraph *mg = new TMultiGraph();
	mg->Add(anglesgr,"p");
	mg->Add(anglesgr2,"p");
	mg->Add(anglesgr3,"p");
	mg->Add(anglesgr4,"p");
	mg->Draw("a");
}

/* emacs
 * Local Variables:
 * mode:C++
 * mode:font-lock
 * c-file-style: "stroustrup"
 * tab-width: 4
 * End:
 */

 
 
