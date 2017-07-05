#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#define LINEWIDTH 2
#define MARKERSIZE 1.5

void plot_cc()
{
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetMarkerSize(1.5);
	TCanvas *c1=new TCanvas("c1","",0,0,1600,900);
	c1->SetLogy();
	c1->SetLeftMargin(0.10);
	c1->SetRightMargin(0.05);
	c1->SetTopMargin(0.05);
	c1->SetBottomMargin(0.15);
	//define TLegend
	TLegend *leg = new TLegend(0.45, 0.77, 0.95, 0.95);
	leg->SetNColumns(3);
	//	leg->SetFillColor(10);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03453388);
	leg->SetFillStyle(0);
	leg->SetTextFont(132);
	// below is c12c12 part

	// exp data for 12C+12C
	TGraphErrors* gr_becker=new TGraphErrors("../becker.dat","%lg %lg %lg");
	TGraphErrors* gr_c12c12ko=new TGraphErrors("../c12c12_kovar.smod","%lg %lg %lg");
	TGraphErrors* gr_c12c12pa=new TGraphErrors("../c12c12_aguilera.smod","%lg %lg %lg");
	// parital data from Spillane with rel_err<40%
	TGraphErrors* gr_c12c12sp=new TGraphErrors("../c12c12_spillane_err_lessthan_0.4.smod","%lg %lg %lg");

	convert_smod_2_s(gr_becker,12.0/2);
	convert_smod_2_s(gr_c12c12ko,12.0/2);
	convert_smod_2_s(gr_c12c12pa,12.0/2);
	convert_smod_2_s(gr_c12c12sp,12.0/2);	
	TF1* func_hind_c12c12 = new TF1("func_hind_c12c12",hindrance_s,0.5,5.0,5);
	func_hind_c12c12->SetParameters(-1.32,52.93,3.68,2.55e16,6.0);
	func_hind_c12c12->SetLineStyle(kDashed);
	func_hind_c12c12->SetLineColor(kRed);

	TMultiGraph* mg_c12c12_smod=new TMultiGraph;
	mg_c12c12_smod->Add(gr_becker, "PZ");	
	mg_c12c12_smod->Add(gr_c12c12ko, "PZ");
	mg_c12c12_smod->Add(gr_c12c12pa, "PZ");
	mg_c12c12_smod->Add(gr_c12c12sp, "PZ");
	//below is the c12c13 part
	//exp data for 12C+13C
	TGraphErrors* gr_c12c13dayras=new TGraphErrors("../c12c13_stokstad_mod.smod","%lg %lg %lg");
	TGraphErrors* gr_c12c13nd=new TGraphErrors("../c12c13_nd_experr.smod","%lg %lg %lg");
	TGraphErrors* gr_c12c13ko=new TGraphErrors("../c12c13_kovar.smod","%lg %lg %lg");
	TGraphErrors* gr_c12c13das=new TGraphErrors("../c12c13_dasmahapatra.smod","%lg %lg %lg");
	TGraphErrors* gr_c12c13anl=new TGraphErrors("../c12c13_rehm_0.05err.smod","%lg %lg %lg %lg");
	TGraphErrors* gr_c12c13imp=new TGraphErrors("../c12c13_imp_experr.smod","%lg %lg %lg");	

	convert_smod_2_s(gr_c12c13dayras,12*13/25.0);
	convert_smod_2_s(gr_c12c13nd,12*13/25.0);
	convert_smod_2_s(gr_c12c13das,12*13/25.);
	convert_smod_2_s(gr_c12c13imp,12*13/25.);
	TMultiGraph* mg_c12c13_smod=new TMultiGraph;
	mg_c12c13_smod->Add(gr_c12c13dayras, "PZ");
	mg_c12c13_smod->Add(gr_c12c13nd, "PZ");
	mg_c12c13_smod->Add(gr_c12c13imp, "PZ");
	mg_c12c13_smod->Add(gr_c12c13das, "PZ");

	TF1* func_hind_c12c13 = new TF1("func_hind_c12c13",hindrance_s,0.5,5.0,5);
	func_hind_c12c13->SetParameters(-2.32,59.37,3.45,6.23e16,12*13/25.);
	func_hind_c12c13->SetLineStyle(kDashed);
	func_hind_c12c13->SetLineColor(kRed);	


	{	// set colors for experimental data
		// c12+c12
		gr_becker->SetLineWidth(2);
		gr_becker->SetLineColor(kMagenta);
		gr_becker->SetMarkerColor(kMagenta);
		gr_becker->SetMarkerStyle(24);

		gr_c12c12ko->SetLineWidth(2);
		gr_c12c12ko->SetLineColor(1);
		gr_c12c12ko->SetMarkerColor(1);
		gr_c12c12ko->SetMarkerStyle(24);

		gr_c12c12pa->SetLineWidth(2);
		gr_c12c12pa->SetLineColor(kGreen+2);
		gr_c12c12pa->SetMarkerColor(kGreen+2);
		gr_c12c12pa->SetMarkerStyle(24);	

		gr_c12c12sp->SetLineWidth(2);
		gr_c12c12sp->SetLineColor(kBlue);
		gr_c12c12sp->SetMarkerColor(kBlue);
		gr_c12c12sp->SetMarkerStyle(24);

		//c12+c13
		gr_c12c13nd->SetLineWidth(2);
		gr_c12c13nd->SetLineColor(kRed);
		gr_c12c13nd->SetMarkerColor(kRed);
		gr_c12c13nd->SetMarkerStyle(20);

		gr_c12c13dayras->SetLineWidth(2);
		gr_c12c13dayras->SetLineColor(kBlue);
		gr_c12c13dayras->SetMarkerColor(kBlue);
		gr_c12c13dayras->SetMarkerStyle(20);

		gr_c12c13ko->SetLineWidth(2);
		gr_c12c13ko->SetLineColor(49);
		gr_c12c13ko->SetMarkerColor(49);
		gr_c12c13ko->SetMarkerStyle(20);

		gr_c12c13das->SetLineWidth(2);
		gr_c12c13das->SetLineColor(kGreen+2);
		gr_c12c13das->SetMarkerColor(kGreen+2);
		gr_c12c13das->SetMarkerStyle(20);

		gr_c12c13anl->SetLineWidth(2);
		gr_c12c13anl->SetLineColor(6);
		gr_c12c13anl->SetMarkerColor(6);
		gr_c12c13anl->SetMarkerStyle(33);

		gr_c12c13imp->SetLineWidth(2);
		gr_c12c13imp->SetLineColor(kBlack);
		gr_c12c13imp->SetMarkerColor(kBlack);
		gr_c12c13imp->SetMarkerStyle(20);

		// use plot_c12c13_exp.C to determine the normalization ratio	
		// DASMAHAPATA/Dayras=1.11 +/- 0.07
		scale_gr_err(gr_c12c13das,1.0/1.11);
		cout<<"DASMAHAPATA scaling factor : "<<1.0/1.11<<endl;
		// IMP/nd=1.119 +/- 0.23, nd/Dayras=0.93 +/-0.10
		scale_gr_err(gr_c12c13imp,1/1.119*1/0.93);
		cout<<"IMP scaling factor : "<<1.0/(0.93*1.119)<<endl;	
		// use plot_c12c13_exp.C to determine the normalization ratio
		// nd/Dayras=0.93 +/- 0.10	
		scale_gr_err(gr_c12c13nd,1.0/0.93);
		cout<<"ND scaling factor : "<<1.0/0.93<<endl;	
	}	

	// below is c13c13 part	
	TGraphErrors* gr_c13c13tren_smod=new TGraphErrors("c13c13_trentalange.smod","%lg %lg %lg");
	scale_gr_err(gr_c13c13tren_smod,1.0/0.84);
	gr_c13c13tren_smod->SetMarkerStyle(22);
	convert_smod_2_s(gr_c13c13tren_smod,13.0/2);

	Double_t xmin_global=1, xmax_global=6.5;
	TH2F *h2_sfactor=new TH2F("h2_sfactor","",100,xmin_global, xmax_global,1000,0.3e15,20e17);
	h2_sfactor->GetXaxis()->SetTitle("E_{c.m.} (MeV)");
	h2_sfactor->GetXaxis()->CenterTitle(true);
	h2_sfactor->GetXaxis()->SetLabelFont(132);
	h2_sfactor->GetXaxis()->SetLabelSize(0.07);
	h2_sfactor->GetXaxis()->SetTitleSize(0.08);
	h2_sfactor->GetXaxis()->SetTitleOffset(0.78);
	h2_sfactor->GetXaxis()->SetTitleFont(132);
	h2_sfactor->GetYaxis()->SetTitle("S factor (MeV b)");
	h2_sfactor->GetYaxis()->CenterTitle(true);
	h2_sfactor->GetYaxis()->SetLabelFont(132);
	h2_sfactor->GetYaxis()->SetLabelSize(0.07);
	h2_sfactor->GetYaxis()->SetTitleSize(0.08);
	h2_sfactor->GetYaxis()->SetTitleOffset(0.6);
	h2_sfactor->GetYaxis()->SetTitleFont(132);


	TMultiGraph* mg_c13c13_smod=new TMultiGraph;
	mg_c13c13_smod->Add(gr_c13c13tren_smod, "PZ");

	h2_sfactor->Draw();

	mg_c12c12_smod->Draw("p");
	func_hind_c12c12->Draw("same");

	Double_t xmin_fit, xmax_fit;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	cout<<"Fitting 12C+13C ........."<<endl<<endl;
	xmin_fit=2.0, xmax_fit=4.5;
	cout<<"Fitting c12+c13 data *****************"<<endl;
	//	Fit_low_E(mg_c12c13_smod,func_hind_c12c13, 12*13/25.0, xmin_fit,xmax_fit);
	//	Fit_low_E_c12c13(mg_c12c13_smod,func_hind_c12c13, 12*13/25.0, xmin_fit,xmax_fit);
	func_hind_c12c13->SetRange(0.5,xmax_fit);
	//	func_hind_c12c13->Draw("same");

	mg_c12c13_smod->Draw("P");
	//	func_hind_c12c13->Draw("same");

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	cout<<"Fitting 13C+13C ........."<<endl<<endl;
	xmin_fit=2.0, xmax_fit=5.0;
	TF1* func_hind_c13c13 = new TF1("func_hind_c13c13",hindrance_s,xmin_fit,xmax_fit,5);
//	Fit_low_E(mg_c13c13_smod,func_hind_c13c13, 13.0/2, xmin_fit,xmax_fit);
//	func_hind_c13c13->SetRange(0.5,xmax_fit);	

//	mg_c13c13_smod->Draw("P");
	//	func_hind_c13c13->Draw("same");	

	leg->AddEntry(gr_c12c12sp, "^{12}C+^{12}C(SP07)", "P");		
	leg->AddEntry(gr_c12c12pa, "^{12}C+^{12}C(AG06)", "P");		
	leg->AddEntry(gr_becker, "^{12}C+^{12}C(BK81)", "P");		
//	leg->AddEntry(gr_c12c12ko, "^{12}C+^{12}C(Kovar)", "P");		
	leg->AddEntry(gr_c12c13dayras, "^{12}C+^{13}C(DA76)", "P");
	leg->AddEntry(gr_c12c13nd, "^{12}C+^{13}C(NO11)", "P");
	leg->AddEntry(gr_c12c13das, "^{12}C+^{13}C(Das)", "P");
	leg->AddEntry(gr_c12c13imp, "^{12}C+^{13}C(IMP17)", "P");
	leg->AddEntry(gr_c13c13tren_smod, "^{13}C+^{13}C(TR88)", "P");
	leg->AddEntry(func_hind_c12c12, "Jiang", "L");
	plot_m3y(leg);
	plot_tdhf(leg);
	plot_spp(leg);
	plot_esw(leg);

// spp from fig1
	TGraph *ufc12c12spp = new TGraph("gasques/spp.sfactor", "%lg %lg");
	ufc12c12spp->Draw("C");
	ufc12c12spp->SetLineWidth(2);
	ufc12c12spp->SetLineStyle(kDashed);
	leg->AddEntry(ufc12c12spp, "gasques", "l");
	plot_hind();
	//	cal_hind();

	//	plot_hind_mc();
	leg->Draw();
}

void plot_hind_mc() {
	Double_t xmin_global=1, xmax_global=6.5;
	TH2F* h2_mc = new TH2F("h2_mc","",100,xmin_global, xmax_global,10000,5,20);

	Double_t a0_mean=-2.32, a0_err=0.24;
	Double_t b0_mean=59.37, b0_err=2.25;
	Double_t es_mean=3.45,  es_err=0.37;
	Double_t smod_mean=6.23e16;
	//	Double_t xfit_max=3.6;
	Double_t xfit_max=4.6;
	TF1* func_hind = new TF1("func_hind_b",hindrance_s,0.5,xfit_max,5);
	func_hind->SetParameters(a0_mean,b0_mean,es_mean,smod_mean,12*13/25.);
	func_hind->FixParameter(4,12*13/25.);

	//    TGraphErrors* gr_c12c13dayras=new TGraphErrors("../c12c13_stokstad_mod_nd_imp_exp_err.smod","%lg %lg %lg");
	//    TGraphErrors* gr_c12c13dayras_orig=new TGraphErrors("../c12c13_stokstad_mod_nd_imp_exp_err.smod","%lg %lg %lg");
	TGraphErrors* gr_c12c13dayras=new TGraphErrors("../c12c13_stokstad_mod.smod","%lg %lg %lg");
	TGraphErrors* gr_c12c13dayras_orig=new TGraphErrors("../c12c13_stokstad_mod.smod","%lg %lg %lg");	
	convert_smod_2_s(gr_c12c13dayras,12*13/25.0);	
	convert_smod_2_s(gr_c12c13dayras_orig,12*13/25.0);

	// create faked experimental data
	Double_t* xx_orig = gr_c12c13dayras_orig->GetX();
	Double_t* yy_orig = gr_c12c13dayras_orig->GetY();
	Double_t* ey_orig = gr_c12c13dayras_orig->GetEY();	

	Double_t* yy=gr_c12c13dayras->GetY();

	for(int i=0;i<1000;i++) {
		for(int k=0;k<gr_c12c13dayras->GetN();k++) 
			yy[k]=gRandom->Gaus(yy_orig[k],sqrt(pow(ey_orig[k],2)+pow(yy_orig[k]*0.0,2)));

		gr_c12c13dayras->Fit(func_hind,"r0");

		for(int k=0;k<100;k++) {
			Double_t x=h2_mc->GetXaxis()->GetBinCenter(k);
			Double_t tempy=func_hind(x);
			Double_t ybin=h2_mc->GetYaxis()->FindBin(log(tempy)/log(10.));	
			Int_t nn=h2_mc->GetBinContent(k,ybin)+1;
			h2_mc->SetBinContent(k,ybin,nn);
		}			
	}

	TGraphErrors* gr_err=new TGraphErrors(100);
	Double_t* x_pt=gr_err->GetX();
	Double_t* y_pt=gr_err->GetY();
	Double_t* ey_pt=gr_err->GetEY();
	for(int i=0;i<100;i++) {
		TH1F* h1=(TH1F*)h2_mc->ProjectionY("h1",i,i);
		x_pt[i]=h2_mc->GetXaxis()->GetBinCenter(i);
		if(x_pt[i]>xfit_max) break;
		y_pt[i]=pow(10,h1->GetMean());
		ey_pt[i]=pow(10,h1->GetMean()+h1->GetRMS())-y_pt[i];
	}

	gr_err->SetMarkerStyle(29);
	gr_err->Draw("P");
}

void plot_hind() {
	// hindrance fit to Dayras data
	TF1* func_hind_b = new TF1("func_hind_b",hindrance_s,0.5,5.0,5);
	func_hind_b->SetParameters(-2.32,59.37,3.45,6.23e16,12*13/25.);
	func_hind_b->SetLineColor(kRed);
	func_hind_b->SetLineStyle(kDashed);
	func_hind_b->Draw("same");	

	// hindrance fit to nd+dayras	
	TF1* func_hind_c = new TF1("func_hind_c",hindrance_s,0.5,5.0,5);
	func_hind_c->SetParameters(-1.722,54.13,2.397,2.43e16,12*13/25.);
	func_hind_c->SetLineColor(kBlue);	
	//	func_hind_c->Draw("same");

}

void cal_hind() {
	TF1* func_hind_b = new TF1("func_hind_b",hindrance_xsec,0.5,5.0,5);
	func_hind_b->SetParameters(-2.32,59.37,3.45,6.23e16,12*13/25.);
	cout<<"Fusion xsec (hind) at 2 MeV : "<<func_hind_b->Eval(2.0)<<" mb"<<endl;
	cout<<"Fusion xsec (hind) at 3 MeV : "<<func_hind_b->Eval(3.0)<<" mb"<<endl;
	TGraph* hind_new=new TGraph("../hindrance_new.smod","%lg %*s %lg");
	TMultiGraph* mg_temp=new TMultiGraph;
	mg_temp->Add(hind_new);
	Fit_low_E(mg_temp,func_hind_b, 12*13/25.0, 1,5);
}

void plot_m3y(TLegend *leg) {
	TGraph* gr_c12c12 = new TGraph("../henning/henning_c12c12_cc_modified.dat","%lg %lg");
	convert_xsec_2_s(gr_c12c12,6.0,1e-3);
	TGraph* gr_c12c13 = new TGraph("../henning/henning_c12c13_cc_modified.dat","%lg %lg");
	convert_xsec_2_s(gr_c12c13,12*13/25.,1e-3);
	TGraph* gr_c13c13 = new TGraph("../henning/henning_c13c13_cc_modified.dat","%lg %lg");
	convert_xsec_2_s(gr_c13c13,6.5,1e-3);	

	gr_c12c12->SetLineColor(kBlue);
	gr_c12c13->SetLineColor(kBlue);
	gr_c13c13->SetLineColor(kBlue);

	gr_c12c12->SetLineWidth(2);
	gr_c12c13->SetLineWidth(2);
	gr_c13c13->SetLineWidth(2);

	gr_c12c12->Draw("C");
	gr_c12c13->Draw("C");
//	gr_c13c13->Draw("C");
	//	gr_c12c12->Print("all");
	leg->AddEntry(gr_c12c12, "M3Y+Rep", "L");
}

void plot_spp(TLegend *leg) {
	TGraph* gr_c12c12 = new TGraph("../spp/c12c12/cross.dat","%lg %lg");
	convert_xsec_2_s(gr_c12c12,6.0,1.);
	TGraph* gr_c12c13 = new TGraph("../spp/c12c13/cross.dat","%lg %lg");
	convert_xsec_2_s(gr_c12c13,12*13/25.,1.);
	TGraph* gr_c13c13 = new TGraph("../spp/c13c13/cross.dat","%lg %lg");
	convert_xsec_2_s(gr_c13c13,6.5,1.);	

	gr_c12c12->SetLineColor(kRed);
	gr_c12c13->SetLineColor(kRed);
	gr_c13c13->SetLineColor(kRed);

	gr_c12c12->SetLineWidth(2);
	gr_c12c13->SetLineWidth(2);
	gr_c13c13->SetLineWidth(2);


	gr_c12c12->Draw("C");
	gr_c12c13->Draw("C");
//	gr_c13c13->Draw("C");
	leg->AddEntry(gr_c12c12, "SPP", "L");
}

void plot_tdhf(TLegend *leg) {
	TGraph* gr_c12c12 = new TGraph("../tdhf/x_12C_12C.out","%lg %lg");
	convert_xsec_2_s(gr_c12c12,6.0,1e-3);
	TGraph* gr_c12c13 = new TGraph("../tdhf/x_12C_13C.out","%lg %lg");
	convert_xsec_2_s(gr_c12c13,12*13/25.,1e-3);

	gr_c12c12->SetLineColor(kGreen+2);
	gr_c12c13->SetLineColor(kGreen+2);

	gr_c12c12->SetLineWidth(2);
	gr_c12c13->SetLineWidth(2);
	// gr_c13c13->SetLineWidth(LINEWIDTH);

	gr_c12c12->Draw("C");
	gr_c12c13->Draw("C");
	leg->AddEntry(gr_c12c12, "TDHF", "L");
}

void plot_esw(TLegend *leg) {
	TGraph* gr_c12c12 = new TGraph("../esw/c12c12_esw_fowler.xsec","%lg %lg");
	convert_xsec_2_s(gr_c12c12,6.0,1.);
	TGraph* gr_c12c13 = new TGraph("../esw/c12c13_esw_dayras.xsec","%lg %lg");
	convert_xsec_2_s(gr_c12c13,12*13/25.,1.);
	TGraph* gr_c13c13 = new TGraph("../esw/c13c13_1.2times.xsec","%lg %lg");
	convert_xsec_2_s(gr_c13c13,6.5,1e-3);	

	gr_c12c12->SetLineColor(kBlack);
	gr_c12c13->SetLineColor(kBlack);
	gr_c13c13->SetLineColor(kBlack);

	gr_c12c12->SetLineWidth(2);
	gr_c12c13->SetLineWidth(2);
	gr_c13c13->SetLineWidth(2);


//	gr_c12c12->Draw("C");
	gr_c12c13->Draw("C");
//	gr_c13c13->Draw("C");
	leg->AddEntry(gr_c12c12, "ESW", "L");
}

void Fit_low_E(TMultiGraph* mg_temp, TF1* func, Double_t mu, Double_t xmin, Double_t xmax) {
	// Creates a Root function based on function
	func->SetLineColor(kRed);
	func->SetLineStyle(kDashed);
	// Sets initial values and parameter names
	//   func->SetParameters(-2.32,59.37,3.45,3e16,mu);
	func->SetParameters(-1.722,54.13,2.397,2.43e16,mu);
	func->SetParNames("A0","B0","Es","Smod_s","mu");
	func->SetRange(xmin,xmax);
	func->FixParameter(4, mu);  

	// set limits for A0   
	Double_t xfactor=2;
	Double_t xmean=func->GetParameter("A0"); 
	Double_t xsigma=xfactor*0.24;
	func->SetParLimits(0,xmean-xsigma,xmean+xsigma);
	// set limits for B0   
	xfactor=2;
	xmean=func->GetParameter("B0");
	xsigma=xfactor*2.25;
	func->SetParLimits(1,xmean-xsigma,xmean+xsigma);
	// Fit histogram in range defined by function   
	mg_temp->Fit(func->GetName(),"r0");	
	cout<<"First fit finished."<<endl;

	// set limits for A0   
	xfactor=25;
	xmean=func->GetParameter("A0");
	xsigma=xfactor*0.24;
	func->SetParLimits(0,xmean-xsigma,xmean+xsigma);
	// set limits for B0   
	xfactor=20;
	xmean=func->GetParameter("B0");
	xsigma=xfactor*2.25;
	func->SetParLimits(1,xmean-xsigma,xmean+xsigma);
	// set limits for energies
	func->SetParLimits(2,1.5,4.0);
	// fix parameter 5: reduced mass (mu)
	func->FixParameter(4, mu);   
	// Fit histogram in range defined by function   
	mg_temp->Fit(func->GetName(),"r0");

	func->Draw("same");
}

void Fit_low_E_c12c13(TMultiGraph* mg_temp, TF1* func, Double_t mu, Double_t xmin, Double_t xmax) {
	// Creates a Root function based on function
	func->SetLineColor(kRed);
	func->SetLineStyle(kDashed);
	// Sets initial values and parameter names
	func->SetParameters(-2.32,59.37,3.45,3e16,mu);
	func->SetParNames("A0","B0","Es","Smod_s","mu");
	func->SetRange(xmin,xmax);
	func->FixParameter(4, mu); 

	func->FixParameter(0, -2.796); 
	func->FixParameter(1, 63.521); 
	// Fit histogram in range defined by function   
	mg_temp->Fit(func->GetName(),"r0E");	
	cout<<"+++++++++++++++ Fit finished."<<endl;

	/*   
	// set limits for A0   
	xfactor=25;
	xmean=-2.32, xsigma=xfactor*0.24;
	func->SetParLimits(0,xmean-xsigma,xmean+xsigma);
	// set limits for B0   
	xfactor=20;
	xmean=49.37, xsigma=xfactor*2.25;
	func->SetParLimits(1,xmean-xsigma,xmean+xsigma);
	// set limits for energies
	func->SetParLimits(2,1.5,4.0);
	// fix parameter 5: reduced mass (mu)
	func->FixParameter(4, mu);   
	// Fit histogram in range defined by function   
	mg_temp->Fit(func->GetName(),"r0");
	*/   
	func->Draw("same");
}
void convert_smod_2_s(TGraphErrors* gr1, Double_t mu) {
	// Get the pointer for gr_1
	Double_t *x, *y, *ex, *ey;
	x=gr1->GetX();
	y=gr1->GetY();
	ex=gr1->GetEX();
	ey=gr1->GetEY();
	// define a constant relative error
	Double_t rel_err=0.02;

	Int_t nn=gr1->GetN();
	Double_t factor=31.29*36/sqrt(1000.0/mu);

	for(int i=0;i<nn;i++) {
		Double_t tempy=(87.21-factor)/sqrt(x[i])+0.46*x[i];
		//		cout<<x[i]<<" "<<tempy<<" "<<myexp(-tempy)<<endl;
		y[i] = y[i]*myexp(-tempy);
		ey[i]=ey[i]*myexp(-tempy);
		ey[i]=sqrt(ey[i]*ey[i]+pow(y[i]*rel_err,2));
	}		
}

void convert_xsec_2_s(TGraph* gr1, Double_t mu, Double_t scaling_factor) {
	// Get the pointer for gr_1
	Double_t *x, *y;
	x=gr1->GetX();
	y=gr1->GetY();

	Int_t nn=gr1->GetN();
	Double_t factor=31.29*36/sqrt(1000.0/mu);

	for(int i=0;i<nn;i++) {
		Double_t tempy=factor/sqrt(x[i]);
		y[i] = scaling_factor*y[i]*x[i]*myexp(tempy);
	}			
}

Double_t hindrance_s(Double_t *x, Double_t *par)
{
	Double_t A0=par[0];
	Double_t B0=par[1];
	Double_t Es=par[2];
	Double_t smod_0=par[3];
	Double_t mu=par[4];
	Double_t factor=31.29*36/sqrt(1000.0/mu);

	Double_t zeta=36/sqrt(mu);
	//   Es=pow((0.495*zeta-B0)/A0,2.0/3);

	Double_t E=x[0];
	Double_t nn=1.5;
	// use mb as unit   
	Double_t sigma_s=smod_0/Es/myexp(87.21/sqrt(Es)+0.46*Es)*1000;   

	Double_t f1=A0*(E-Es)-B0*(pow(Es/E,nn-1)-1)/pow(Es,nn-1)/(nn-1);
	if(f1>700) f1=700;
	Double_t f2=sigma_s*Es/E*myexp(f1);
	// convert back to standard s-factor   
	Double_t f3=f2/1000*E*myexp(factor/sqrt(E)+0.0*E);
	return f3;
}

Double_t hindrance_xsec(Double_t *x, Double_t *par)
{
	Double_t A0=par[0];
	Double_t B0=par[1];
	Double_t Es=par[2];
	Double_t smod_0=par[3];
	Double_t mu=par[4];
	Double_t factor=31.29*36/sqrt(1000.0/mu);

	Double_t zeta=36/sqrt(mu);
	//   Es=pow((0.495*zeta-B0)/A0,2.0/3);

	Double_t E=x[0];
	Double_t nn=1.5;
	// use mb as unit   
	Double_t sigma_s=smod_0/Es/myexp(87.21/sqrt(Es)+0.46*Es)*1000;   

	Double_t f1=A0*(E-Es)-B0*(pow(Es/E,nn-1)-1)/pow(Es,nn-1)/(nn-1);
	if(f1>700) f1=700;
	Double_t f2=sigma_s*Es/E*myexp(f1);
	// convert back to mb
	Double_t f3=f2;
	return f3;
}

void scale_gr_err(TGraphErrors* gr1, Double_t scale) {
	// Get the pointer for gr_1
	Double_t *x, *y, *ex, *ey;
	x=gr1->GetX();
	y=gr1->GetY();
	ex=gr1->GetEX();
	ey=gr1->GetEY();

	Int_t nn=gr1->GetN();	
	for(int i=0;i<nn;i++) {
		y[i] = y[i]*scale;
		ey[i]=ey[i]*scale;
	}		
}

Double_t myexp(Double_t x) {
	if(x>700) x=700;
	return exp(x);
}
