#define MARKERSIZE 1.5
#define LINEWIDTH 3
#define HBARC 197.326
#define MAXPOINT 100
#include <iostream>
#include "TStyle.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TColor.h"
#include "TAxis.h"
#include "Rtypes.h"

using namespace std;
void cal_uf(TGraph *, double, double, double, double, char *, double, char *, double, double, double, double, char *);
void shadow(const TGraph*, const TGraph*, const TGraph*, TMultiGraph *, short);
void shadow(const TGraph*, const TGraph*, const TGraph*, TMultiGraph *, short, short, short, double, double, double);
void cal_uf_grerr(Double_t, Double_t, Double_t, Double_t, TGraphErrors *, Double_t , char *, Double_t, Double_t, Double_t, Double_t, char *); 
Double_t ff_2_smod_ref(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t); 
Double_t ff(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t); 
void plot_exp_c12c12_data(TLegend *);
void plot_exp_c12c13_data(TLegend *);
void plot_exp_c13c13_data(TLegend *);
void modsmfeld(TGraph *, double);
void convert_xsec_2_smod(TGraph *);
void convert_smod_2_xsec(TGraphErrors *);
void grmod(TGraph *, double);
double max(double, double, double);
double min(double, double, double);
void plot_uf_normalize() 
{
	const float small = 1e-5;
	gStyle->SetOptStat(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetTitleFont(132);
	gStyle->SetLabelFont(132);
	gStyle->SetLegendFont(132);
	//define the canvas
	TCanvas *c1 = new TCanvas("c1", "UFF", 0, 0, 1000, 600);
//	c1->SetGrid();
//	c1->SetLeftMargin(0.01);
	c1->SetRightMargin(0.01);
	c1->SetTopMargin(0.05);
	c1->SetBottomMargin(0.15);

	//define the multigraph
	TMultiGraph *mg = new TMultiGraph();
	mg->Draw("AC");

	//define the legend	
	TLegend *leg = new TLegend(0.15, 0.68, 0.98, 0.95, NULL, "brNDC");
	leg->SetNColumns(3);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.047);
	leg->SetFillStyle(0);
	TLegend *leg2 = new TLegend(0.68, 0.55, 0.95, 0.98, NULL, "brNDC");
	leg2->SetNColumns(1);
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.055);
	leg2->SetFillStyle(0);

	leg->Draw();
	
	TGraph *ufc12c12_henning = new TGraph(100);
	TGraph *ufc12c13_henning = new TGraph(100);
	TGraph *ufc13c13_henning = new TGraph(100);

	// define the reference system
	const Double_t z1ref = 6, z2ref = 6, a1ref = 12, a2ref = 12;
	Double_t z1, z2, a1, a2;

	//henning's m3y+rep

	z1 = 6, z2 = 6, a1 = 12, a2 = 12;
	cal_uf(ufc12c12_henning, z1, a1, z2, a2, "henning/henning_c12c12_cc.dat", 0.73e-3, "c12c12/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");

	z1 = 6, z2 = 6, a1 = 12, a2 = 13;  
	cal_uf(ufc12c13_henning, z1, a1, z2, a2, "henning/henning_c12c13_cc.dat", 0.73e-3, "c12c13/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");

//	z1 = 6, z2 = 6, a1 = 13, a2 = 13;
//	cal_uf(ufc13c13_henning ,z1,a1,z2,a2,"henning/henning_c13c13_cc_modified.dat",0.9*1e-3,"c13c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");


	ufc12c12_henning->SetLineColor(kRed);
	ufc12c12_henning->SetLineWidth(LINEWIDTH);	
	ufc12c13_henning->SetLineColor(kGreen+2);
	ufc12c13_henning->SetLineWidth(LINEWIDTH);
//	ufc13c13_henning->SetLineColor(kBlue);
//	ufc13c13_henning->SetLineWidth(LINEWIDTH);

//	mg->Add(ufc12c12_henning);
//	mg->Add(ufc12c13_henning);
//	mg->Add(ufc13c13_henning);

// tdhf
	TGraph *ufc12c12_tdhf = new TGraph(100);
	TGraph *ufc12c13_tdhf = new TGraph(100);
//	TGraph *ufc13c13_tdhf = new TGraph(100);

	z1=6, z2=6, a1=12, a2=12;
	cal_uf(ufc12c12_tdhf,z1,a1,z2,a2,"tdhf/x_12C_12C_amc2_931.out",0.67e-3,"c12c12/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

	z1=6, z2=6, a1=12, a2=13;
	cal_uf(ufc12c13_tdhf,z1,a1,z2,a2,"tdhf/x_12C_13C_amc2_931.out",0.67e-3,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

	ufc12c12_tdhf->SetLineColor(kBlue);
	ufc12c12_tdhf->SetLineWidth(LINEWIDTH);	
	ufc12c13_tdhf->SetLineColor(kMagenta);
	ufc12c13_tdhf->SetLineWidth(LINEWIDTH);
//	ufc13c13_tdhf->SetLineColor(kBlue-6);
//	ufc13c13_tdhf->SetLineWidth(LINEWIDTH);

//	mg->Add(ufc12c12_tdhf);
//	mg->Add(ufc12c13_tdhf);
//	mg->Add(ufc13c13_tdhf);

	// esw models
	TGraph *ufc12c13_esw = new TGraph(100);
	TGraph *ufc13c13_esw = new TGraph(100);
	
	z1 = 6, z2 = 6, a1 = 12, a2 = 13;  
	cal_uf(ufc12c13_esw,z1,a1,z2,a2,"esw/c13c12_all_imp.xsec",1e-3,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

//	cal_uf(ufc12c13_esw,z1,a1,z2,a2,"esw/c12c13_esw_imp_stock.xsec",1e-3,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

	z1 = 6, z2 = 6, a1 = 13, a2 = 13;
//	cal_uf(ufc13c13_esw, z1, a1, z2, a2, "esw/c13c13_1.2times_emin_0.1.xsec_cut", 1e-3, "c13c13/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");

	// plot esw 
	ufc12c13_esw->SetLineColor(kBlack);
	ufc12c13_esw->SetLineWidth(5);
	ufc12c13_esw->SetLineStyle(kDashed);
//	ufc13c13_esw->SetLineColor(kMagenta);
//	ufc13c13_esw->SetLineWidth(5);
//	ufc13c13_esw->SetLineStyle(2);
	mg->Add(ufc12c13_esw);
//	mg->Add(ufc13c13_esw);

	//kns
	TGraph *ufc12c12kns = new TGraph("kns.smod", "%lg %lg");
	modsmfeld(ufc12c12kns, 87.21);
	ufc12c12kns->SetLineColor(kBlue);
	ufc12c12kns->SetLineWidth(LINEWIDTH);
	ufc12c12kns->SetLineStyle(kDashed);
	mg->Add(ufc12c12kns, "C");

	//spp
	TGraph *ufc12c12spp = new TGraph("c12c12/cross_back.dat", "%lg %lg");
	convert_xsec_2_smod(ufc12c12spp);
	grmod(ufc12c12spp, 4.21 /6.24);
//	TGraph *ufc12c12spp = new TGraph("gasques/spp.sfactor", "%lg %lg");
//	sfactor2smod(ufc12c12spp);
//	TGraph *ufc12c12spp = new TGraph(100);
//	z1 = 6, z2 = 6, a1 = 12, a2 = 12;  
//	cal_uf(ufc12c12spp,z1,a1,z2,a2,"c12c12/cross.dat",1, "c12c12/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	ufc12c12spp->SetLineColor(kRed);
	ufc12c12spp->SetLineWidth(LINEWIDTH);
	ufc12c12spp->SetLineStyle(kDashed);
	mg->Add(ufc12c12spp, "C");


	// get the variation of the model prediction
	for(int i = 0; i < ufc12c12_henning->GetN(); i++)
	{          
		Double_t tempx[3], tempy[3];
		tempx[0] = 0, tempx[1] = 1, tempx[2] = 2;
		tempy[0] = * (ufc12c12_henning->GetY()+i);
		tempy[1] = * (ufc12c13_henning->GetY()+i);
		tempy[2] = * (ufc13c13_henning->GetY()+i);
		TGraph *gr_temp = new TGraph(3, tempx, tempy);
		Double_t diff = gr_temp->GetMaximum() - gr_temp->GetMinimum();
		cout<<*(ufc12c12_henning->GetX() + i)<<" "
			<<tempy[0]<<" "
			<<tempy[1]<<" "
			<<tempy[2]<<" "
			<<diff<<endl;
	}

	mg->GetYaxis()->SetRangeUser(0 ,1.25e17);
	mg->GetXaxis()->SetRangeUser(0,7);
	mg->GetXaxis()->SetTitle("E_{c.m.}(MeV)");
	mg->GetXaxis()->SetTitleSize(0.09);
	mg->GetXaxis()->CenterTitle();
	mg->GetXaxis()->SetLabelSize(0.047);
	mg->GetXaxis()->SetTitleOffset(0.62);
	mg->GetYaxis()->SetTitle("S* (MeVb)");
	mg->GetYaxis()->SetTitleSize(0.08);
	mg->GetYaxis()->SetTitleFont(132);
	mg->GetYaxis()->SetLabelSize(0.047);
	mg->GetYaxis()->SetLabelFont(132);
	mg->GetYaxis()->CenterTitle();
	mg->GetYaxis()->SetTitleOffset(0.4);

	//the shadow area is  up and down limits
//	shadow(ufc12c12_henning, ufc12c13_henning, ufc12c13_henning, mg, kGreen+2);
	shadow(ufc12c12_henning, ufc12c13_henning, ufc12c13_henning, mg, kGreen+2, kGreen + 2, kRed, 0., 3., 1);
//	shadow(ufc12c12_tdhf, ufc12c13_tdhf, ufc12c12_tdhf, mg, kBlue);
	shadow(ufc12c12_tdhf, ufc12c13_tdhf, ufc12c12_tdhf, mg, kBlue, kMagenta, kBlue, 0., 3., 1);
	
	//tex
	TLatex *tex = new TLatex(1.,30e+015,"M3Y+REP");
	tex->Draw();
	tex->SetTextFont(132);
	tex->SetTextSizePixels(34);
	TLatex *tex2 = new TLatex(1.,58e+015,"TDHF");
	tex2->Draw();
	tex2->SetTextFont(132);
	tex2->SetTextSizePixels(34);

	leg->AddEntry(ufc12c12_henning, "^{12}C+^{12}C(M3Y+Rep)", "l");
	leg->AddEntry(ufc12c13_henning, "^{12}C+^{13}C(M3Y+Rep)", "l");
	leg->AddEntry(ufc12c12_tdhf, "^{12}C+^{12}C(TDHF)", "l");
	leg->AddEntry(ufc12c13_tdhf, "^{12}C+^{13}C(TDHF)", "l");
	leg->AddEntry(ufc12c13_esw, "^{12}C+^{13}C(ESW)", "l");
	leg->AddEntry(ufc12c12spp, "^{12}C+^{12}C(SPP)", "l");
	leg->AddEntry(ufc12c12kns, "^{12}C+^{12}C(KNS)", "l");
//	leg->AddEntry(ufc13c13, "^{13}C+^{13}C(M3Y+Rep)", "l");
//	leg->AddEntry(ufc13c13_tdhf, "^{13}C+^{13}C(TDHF)", "l");
//	leg->AddEntry(ufc13c13_esw, "^{13}C+^{13}C(ESW)", "l");
	leg->Draw();	

	plot_exp_c12c12_data(leg);
	plot_exp_c12c13_data(leg);
//	plot_exp_c13c13_data(leg);

	c1->Update();
} 

void test() 
{
	double fftemp = ff(6, 12*938.918/931.494, 6, 12*938.918/931.494, 0.55, 6.68e-32 / 1000) ;
	cout << "ff = " << fftemp << endl;
	cout << "ff = " << ff(6, 12, 6, 12, 0.55, 8.87e-32 / 1000) << endl;
	cout << "Smod" << ff_2_smod_ref(6, 12, 6, 12, 0.55 / 36, fftemp) << endl;
	cout << "Smod ref" << 8.87e-32 / 1000 * 0.55 *exp(87.26/sqrt(0.55) + 0.46 * 0.55) << endl;
}

void plot_exp_c12c12_data(TLegend *leg)
{
	//define the reference system
	Double_t z1ref = 6, z2ref = 6,a1ref = 12,a2ref = 12;

	TGraphErrors* gr_c12c12_bk = new TGraphErrors("exp_data/c12c12/becker.smod", "%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c12_bk);
	TGraphErrors* gr_c12c12_sp = new TGraphErrors("exp_data/c12c12/c12c12_spillane_err_lessthan_0.4.smod", "%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c12_sp);
	TGraphErrors* gr_c12c12_ag = new TGraphErrors("exp_data/c12c12/aguilera.smod", "%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c12_ag);

	gr_c12c12_sp->SetLineWidth(2);
	gr_c12c12_sp->SetLineColor(kBlue+2);
	gr_c12c12_sp->SetMarkerColor(kBlue+2);
	gr_c12c12_sp->SetMarkerStyle(24);
	gr_c12c12_sp->SetMarkerSize(MARKERSIZE);	

	gr_c12c12_ag->SetLineWidth(2);
	gr_c12c12_ag->SetLineColor(kGreen+2);
	gr_c12c12_ag->SetMarkerColor(kGreen+2);
	gr_c12c12_ag->SetMarkerStyle(24);
	gr_c12c12_ag->SetMarkerSize(MARKERSIZE);

	gr_c12c12_bk->SetLineWidth(2);
	gr_c12c12_bk->SetLineColor(kMagenta);
	gr_c12c12_bk->SetMarkerColor(kMagenta);
	gr_c12c12_bk->SetMarkerStyle(24);
	gr_c12c12_bk->SetMarkerSize(MARKERSIZE);

	cal_uf_grerr(6, 12, 6, 12, gr_c12c12_sp, 1, "c12c12/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");
	gr_c12c12_sp->Draw("PZ");
	cal_uf_grerr(6, 12, 6, 12, gr_c12c12_ag, 1, "c12c12/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");
	gr_c12c12_ag->Draw("LPZ");
	cal_uf_grerr(6, 12, 6, 12, gr_c12c12_bk, 1, "c12c12/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");
	gr_c12c12_bk->Draw("LPZ");
	leg->AddEntry(gr_c12c12_sp, "^{12}C+^{12}C(SP07)", "p");
	leg->AddEntry(gr_c12c12_ag, "^{12}C+^{12}C(AG06)", "p");
	leg->AddEntry(gr_c12c12_bk, "^{12}C+^{12}C(BK81)", "p");
}
void plot_exp_c12c13_data(TLegend *leg)
{
	//define the reference system
	Double_t z1ref, z2ref, a1ref, a2ref;
	z1ref = 6, z2ref = 6,a1ref = 12,a2ref = 12;
	//the smod caculation use the factor 87.21
	TGraphErrors* gr_c12c13_st = new TGraphErrors("exp_data/c12c13/stockstad.smod", "%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c13_st);
	TGraphErrors* gr_c12c13_no = new TGraphErrors("exp_data/c12c13/c12c13_nd_experr.smod", "%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c13_no);	
	TGraphErrors* gr_c12c13_imp = new TGraphErrors("exp_data/c12c13/imp_experr.smod", "%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c13_imp);

	gr_c12c13_st->SetLineWidth(2);
	gr_c12c13_st->SetLineColor(kBlue);
	gr_c12c13_st->SetMarkerColor(kBlue);
	gr_c12c13_st->SetMarkerSize(1);
	gr_c12c13_st->SetMarkerStyle(20);
	gr_c12c13_st->SetMarkerSize(MARKERSIZE);	

	gr_c12c13_no->SetLineWidth(2);
	gr_c12c13_no->SetLineColor(kRed);
	gr_c12c13_no->SetMarkerColor(kRed);
	gr_c12c13_no->SetMarkerStyle(20);
	gr_c12c13_no->SetMarkerSize(MARKERSIZE);

	gr_c12c13_imp->SetLineWidth(2);
	gr_c12c13_imp->SetLineColor(kBlack);
	gr_c12c13_imp->SetMarkerColor(kBlack);
	gr_c12c13_imp->SetMarkerStyle(20);
	gr_c12c13_imp->SetMarkerSize(MARKERSIZE);

	cal_uf_grerr(6, 12, 6, 13, gr_c12c13_st, 1, "c12c13/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");
	gr_c12c13_st->Draw("PZ");
	cal_uf_grerr(6, 12, 6, 13, gr_c12c13_no, 1.075, "c12c13/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");
	gr_c12c13_no->Draw("PZ");
	cal_uf_grerr(6, 12, 6, 13, gr_c12c13_imp, 1.41/1.50, "c12c13/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");
	gr_c12c13_imp->Draw("PZ");
	leg->AddEntry(gr_c12c13_st, "^{12}C+^{13}C(DA76)", "p");
	leg->AddEntry(gr_c12c13_no, "^{12}C+^{13}C(NO11)", "p");
	leg->AddEntry(gr_c12c13_imp, "^{12}C+^{13}C(this work)", "p");
}
void plot_exp_c13c13_data(TLegend *leg)
{
	//	define the reference system
	Double_t z1ref, z2ref, a1ref, a2ref;
	z1ref = 6, z2ref = 6,a1ref = 12,a2ref = 12;

	TGraphErrors* gr_c13c13_trent = new TGraphErrors("exp_data/c13c13/trentalange.smod", "%lg %lg %lg");
	convert_smod_2_xsec(gr_c13c13_trent);

	gr_c13c13_trent->SetLineWidth(2);
	gr_c13c13_trent->SetLineColor(kBlack);
	gr_c13c13_trent->SetMarkerColor(kBlack);
	gr_c13c13_trent->SetMarkerStyle(22);
	gr_c13c13_trent->SetMarkerSize(MARKERSIZE);

	cal_uf_grerr(6, 13, 6, 13, gr_c13c13_trent, 1.2, "c13c13/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");
	gr_c13c13_trent->Draw("P");
	leg->AddEntry(gr_c13c13_trent, "^{13}C+^{13}C(TR88)", "p");
}

// convert Smod into xsec in unit of barn
void convert_smod_2_xsec(TGraphErrors* gr1)
{
	// Get the pointer for gr_1
	Double_t *x, *y, *ex, *ey;
	x = gr1->GetX();
	y = gr1->GetY();
	ex = gr1->GetEX();
	ey = gr1->GetEY();

	Int_t nn = gr1->GetN();

	for(int i = 0; i < nn; i++)
	{
		Double_t tempy = 87.26 / sqrt(x[i]) + 0.46 * x[i];
		y[i] = y[i] * exp(-tempy) / x[i];
		ey[i] = ey[i] * exp(-tempy) / x[i];
	}		
}

void convert_xsec_2_smod(TGraphErrors* gr1)
{
	// Get the pointer for gr_1
	Double_t *x, *y, *ex, *ey;
	x = gr1->GetX();
	y = gr1->GetY();
	ex = gr1->GetEX();
	ey = gr1->GetEY();

	Int_t nn = gr1->GetN();

	for(int i = 0; i < nn; i++)
	{
		Double_t tempy = 87.26 / sqrt(x[i]) + 0.46 * x[i];
		y[i] = y[i] * exp(tempy) * x[i];
		ey[i] = ey[i] * exp(tempy) * x[i];
	}		
}

void convert_xsec_2_smod(TGraph* gr1)
{
	// Get the pointer for gr_1
	Double_t *x, *y;
	x = gr1->GetX();
	y = gr1->GetY();

	Int_t nn = gr1->GetN();

	for(int i = 0; i < nn; i++)
	{
		Double_t tempy = 87.26 / sqrt(x[i]) + 0.46 * x[i];
		y[i] = y[i] * exp(tempy) * x[i];
	}		
}

void cal_uf_grerr(Double_t z1, Double_t a1, Double_t z2, Double_t a2, TGraphErrors* uf, Double_t f2_factor, char* f1, Double_t  z1ref, Double_t a1ref, Double_t z2ref, Double_t a2ref, char* f0) 
{
	TGraphErrors* FF2 = uf;
	TGraph* FF1 = new TGraph(f1);
	TGraph* FF0 = new TGraph(f0);

	cout<<FF2->GetN()<<endl;
	cout<<f1<<" "<<FF1->GetN()<<endl;
	cout<<f0<<" "<<FF0->GetN()<<endl;

	Double_t *x, *y,*ex,*ey;
	x=FF1->GetX();
	y=FF1->GetY();
	for(int i=0;i<FF1->GetN();i++) {
		y[i]=ff(z1,a1,z2,a2,x[i],y[i]);
		x[i]=fx(z1,a1,z2,a2,x[i]);
	}

	x=FF0->GetX();
	y=FF0->GetY();
	for(int i=0;i<FF0->GetN();i++) {
		y[i]=ff(z1ref,a1ref,z2ref,a2ref,x[i],y[i]);
		x[i]=fx(z1ref,a1ref,z2ref,a2ref,x[i]);
	}

	TSpline *spline=0;
	Int_t num = FF2->GetN();
	x=FF2->GetX();
	y=FF2->GetY();
	ey=FF2->GetEY();
	for(int i=0;i<FF2->GetN();i++) 
	{
		y[i]=y[i]*f2_factor;
		ey[i]=ey[i]*f2_factor;

		if(y[i]<1e-15)
		{
			y[i]=1e-15; 
			ey[i]=1e-15; 
			continue;
		}

		Double_t oldy=y[i];
		y[i]=ff(z1,a1,z2,a2,x[i],y[i]);
		ey[i]=ff(z1,a1,z2,a2,x[i],oldy+ey[i]);
		ey[i]=fabs(ey[i]-y[i]);
		x[i]=fx(z1,a1,z2,a2,x[i]);

		y[i]=y[i]-FF1->Eval(x[i],spline,"S");
		y[i]=y[i]+FF0->Eval(x[i],spline,"S");

		// convert UF back to reference system		
		Double_t uf_x=x[i],uf_y=y[i], uf_yerr=ey[i];
		x[i]=uf_x*z1ref*z2ref;
		y[i]=ff_2_smod_ref(z1ref, a1ref, z2ref, a2ref, uf_x, y[i]);		
		ey[i]=ff_2_smod_ref(z1ref, a1ref, z2ref, a2ref,	uf_x,uf_y+uf_yerr);
		ey[i]=fabs(ey[i]-y[i]);

//		cout<<"before conversion : "<<uf_x[i]<<" "<<uf_y<<" "<<uf_yerr<<endl;
//		cout<<"after  conversion : "<<x[i]<<" "<<y[i]<<" "<<ey[i]<<endl;
	}
}

void cal_uf(TGraph* uf,Double_t z1, Double_t a1, Double_t z2, Double_t a2, char* f2, Double_t f2_factor, char* f1, Double_t  z1ref, Double_t a1ref, Double_t z2ref, Double_t a2ref, char* f0) 
{
	TGraph* FF2=new TGraph(f2);
	TGraph* FF1=new TGraph(f1);
	TGraph* FF0=new TGraph(f0);

	cout<<f2<<" "<<FF2->GetN()<<endl;
	cout<<f1<<" "<<FF1->GetN()<<endl;
	cout<<f0<<" "<<FF0->GetN()<<endl;

	Double_t *x, *y;
	x=FF2->GetX();
	y=FF2->GetY();
	Double_t emin=x[0];
	for(int i=0;i<FF2->GetN();i++) 
	{
		Double_t tempx=x[i];
		Double_t tempy=y[i];
		y[i]=ff(z1,a1,z2,a2,tempx,y[i]*f2_factor);
		x[i]=fx(z1,a1,z2,a2,tempx);

		if(tempx<0.56) cout<<"Ecm = "<<tempx<<" xsec(mb)= "<<tempy
			<<" fx = "<<x[i]<<" ff = "<<y[i]<<endl;
	}

	x=FF1->GetX();
	y=FF1->GetY();
	for(int i=0;i<FF1->GetN();i++) 
	{
		Double_t tempx=x[i];
		y[i]=ff(z1,a1,z2,a2,tempx,y[i]);
		x[i]=fx(z1,a1,z2,a2,tempx);
	}

	x=FF0->GetX();
	y=FF0->GetY();
	for(int i=0;i<FF0->GetN();i++)
	{
		Double_t tempx=x[i];
		y[i]=ff(z1ref,a1ref,z2ref,a2ref,tempx,y[i]);
		x[i]=fx(z1ref,a1ref,z2ref,a2ref,tempx);
	}

	Double_t emax = 10;
	Int_t num = uf->GetN();
	Double_t de = (emax - emin) / num;

	x = uf->GetX();
	y = uf->GetY();

	TSpline *spline=0;
	for(int i=0;i<num;i++) 
	{
		double fe=emin+de*i;
		Double_t fx_temp=fx(z1,a1,z2,a2,fe);
		x[i]=fx_temp;
		y[i]=FF2->Eval(fx_temp,spline,"S");
		y[i]=y[i]-FF1->Eval(fx_temp,spline,"S");
		y[i]=y[i]+FF0->Eval(fx_temp,spline,"S");
		// convert obtained UF to Smod for the reference system	
		x[i]=fx_temp*z1ref*z2ref;
		if(x[i]<0.56) cout<<"Ecm= "<<x[i]<<" uf = "<<y[i]<<endl;
		y[i]=ff_2_smod_ref(z1ref, a1ref, z2ref, a2ref, fx_temp, y[i]);
		if(x[i]<0.56) cout<<" smod = "<<y[i]<<endl;
	}
}

Double_t fx(Double_t z1, Double_t a1, Double_t z2, Double_t a2,Double_t ecm) 
{
	Double_t temp=pow(ecm/z1/z2,1);
	return temp;
}

Double_t ff(Double_t z1, Double_t a1, Double_t z2, Double_t a2,Double_t ecm, Double_t sigma) 
{
	Double_t mu=a1*a2/(a1+a2);
	Double_t fk=sqrt(2*mu*931.5*ecm)/HBARC;
	// sigma in barn
	Double_t pene=sigma*1e2*fk**2/3.1415926;
	Double_t temp=log(pene);		
	Double_t temp=temp/sqrt(mu*z1*z2);
	return temp;
}

Double_t ff_2_smod_ref(Double_t z1ref, Double_t a1ref, Double_t z2ref, Double_t a2ref, Double_t fe_val, Double_t uf_val) 
{	
	Double_t ecm_ref=fe_val*z1ref*z2ref;
	Double_t mu_ref=a1ref*a2ref/(a1ref+a2ref);

	Double_t temp=uf_val;
	temp=temp*sqrt(mu_ref*z1ref*z2ref);
	// get log(sigma*Ecm)
	temp=temp-log(2*mu_ref*931.5/3.1415926/(HBARC*HBARC)*100.0);
	// convert to Smod
	temp=exp(temp+31.29*z1ref*z2ref*sqrt(mu_ref/(ecm_ref*1000))+0.46*ecm_ref);

	return temp;
}

void sfactor2smod(TGraph *gr)
{
	int nn;
	double *x, *y;
	x = gr->GetX();
	y = gr->GetY();
	nn = gr->GetN();
	for(int i = 0; i < nn; ++i)
	{
		y[i] = y[i] * exp(0.46*x[i]);
	}

}

void modsmfeld(TGraph *gr, double r)
{
	int nn;
	double *x, *y;
	x = gr->GetX();
	y = gr->GetY();
	nn = gr->GetN();
	for(int i = 0; i < nn; ++i)
	{
		y[i] = y[i] * exp(87.26-r);
	}

}

void grmod(TGraph *gr, double r)
{
	int n = gr->GetN();
	double *x = gr->GetX();
	double *y = gr->GetY();
	for(int i = 0;  i < n; ++i)
		y[i] *= r;
}

void shadow(const TGraph* tge1,const TGraph* tge2,const TGraph* tge3, TMultiGraph *mg, short color)
{
	Double_t ymax[100], ymin[100],ydown[100];
	Double_t xmax, xmin, dx;
	Double_t x0[100] = {0}, y0[100] = {0}, y1[100], y2[100], y3[100];
	Double_t *x1, *x2, *x3;
	x1 = tge1->GetX();
	x2 = tge2->GetX();
	x3 = tge3->GetX();
	xmax = min(x1[100-1], x2[100-1], x3[100-1]);
	xmin = max(x1[0], x2[0], x3[0]);
	dx = (xmax - xmin) / 100;
	TSpline *spline=0;
	for(int i = 0; i < 100; i++)
	{
		x0[i] = xmin + dx * i;
		y1[i] = tge1->Eval(x0[i], spline, "S");
		y2[i] = tge2->Eval(x0[i], spline, "S");
		y3[i] = tge3->Eval(x0[i], spline, "S");

		ymax[i] = max(y1[i], y2[i], y3[i]);
		ymin[i] = min(y1[i], y2[i], y3[i]);
		ydown[i] = ymax[i] - ymin[i];
	}
	TGraphAsymmErrors *gmax = new TGraphAsymmErrors(100, x0, ymax, y0, y0, ydown, y0);
	TGraph *gmin = new TGraph(100, x0, ymin);
	gmax->SetLineWidth(1);
	gmax->SetLineStyle(2);
//	gmax->SetLineColorAlpha(color, 0.4);
	gmin->SetLineWidth(1);
	gmin->SetLineStyle(2);
//	gmin->SetLineColorAlpha(color, 0.4);
	gmax->SetFillStyle(3002);
//	gmax->SetFillColorAlpha(color, 0.4);
	gmax->SetFillColor(color);
	gmin->SetFillColor(color);
	mg->Add(gmax, "CE3");
	mg->Add(gmin);
}

void shadow(const TGraph* tge1, const TGraph* tge2,const  TGraph* tge3, TMultiGraph *mg, short color, short color1, short color2, double xl, double xh, double xcoeff)
{
	Double_t ymax[100], ymin[100],ydown[100];
	Double_t xmax, xmin, dx;
	Double_t x0[100] = {0}, y0[100] = {0}, y1[100], y2[100], y3[100];
	Double_t *x1, *x2, *x3;
	x1 = tge1->GetX();
	x2 = tge2->GetX();
	x3 = tge3->GetX();
	xmax = min(x1[100-1], x2[100-1], x3[100-1]);
	xmax = xmax<xh?xmax:xh;
	xmin = max(x1[0], x2[0], x3[0]);
	xmin = xmin>xl?xmin:xl;
	dx = (xmax - xmin) / 100;
	TSpline *spline=0;
	for(int i = 0; i < 100; i++)
	{
		x0[i] = xmin + dx * i;
		y1[i] = tge1->Eval(x0[i], spline, "S");
		y2[i] = tge2->Eval(x0[i], spline, "S");
		y3[i] = tge3->Eval(x0[i], spline, "S");

		ymax[i] = max(y1[i], y2[i], y3[i]) * xcoeff;
		ymin[i] = min(y1[i], y2[i], y3[i]) * xcoeff;
		ydown[i] = ymax[i] - ymin[i];
	}
	TGraphAsymmErrors *gmax = new TGraphAsymmErrors(100, x0, ymax, y0, y0, ydown, y0);
	TGraph *gmin = new TGraph(100, x0, ymin);
	gmax->SetLineWidth(LINEWIDTH);
	gmax->SetLineStyle(1);
//	gmax->SetLineColorAlpha(color,0.8);
	gmin->SetLineWidth(LINEWIDTH);
	gmin->SetLineStyle(1);
//	gmin->SetLineColorAlpha(color,0.8);
	gmax->SetFillStyle(3002);
//	gmax->SetFillColorAlpha(color,0.4);
	gmax->SetLineColor(color1);
	gmax->SetFillColor(color);
	gmin->SetLineColor(color2);
	mg->Add(gmax, "CE3");
	mg->Add(gmin);
}

double max(double a,double b,double c) 
{ 
	if (a > b && a > c) return a; 
	if (b > a && b > c) return b; 

	return c; 
}

double min(double a,double b,double c) 
{ 
	if (a < b && a < c) return a; 
	if (b < a && b < c) return b; 

	return c; 
}
