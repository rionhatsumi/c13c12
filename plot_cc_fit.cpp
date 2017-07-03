#define MARKERSIZE 1.5
#define LINEWIDTH 3
#define HBARC 197.326
#include <iostream>
#include "TStyle.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"

int plot_cc_fit()
{
	using namespace std;
	const float small = 1e-5;
	gStyle->SetOptStat(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetTitleFont(132);
	gStyle->SetLabelFont(132);
	gStyle->SetLegendFont(132);
	Int_t n;
	//define the canvas
	TCanvas *c1 = new TCanvas("c1", "UFF", 0, 0, 1200, 1000);
	//	c1->SetGrid();
	c1->SetLeftMargin(0.01);
	c1->SetRightMargin(0.01);
	c1->SetTopMargin(0.01);
	c1->SetBottomMargin(0.3);
	TPad *p1 = new TPad("p1", "m3y", 1e-5, 0.49999, 0.9999, 0.9999);
	TPad *p2 = new TPad("p2", "tdhf", 1e-5, 1e-5, 0.9999, 0.49999);
	p1->Draw();
	p1->SetBottomMargin(small);
	p2->Draw();
	p2->SetTopMargin(small);
	p1->cd();

	TMultiGraph *mg = new TMultiGraph();
	mg->Draw("AC");

	//define the legend	
	TLegend *leg = new TLegend(0.12, 0.7, 0.3, 0.9, NULL, "brNDC");
	leg->SetNColumns(1);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03453388);
	leg->SetFillStyle(0);
	TLegend *leg4 = new TLegend(0.7, 0.65, 0.95, 0.9,NULL,"brNDC");
	leg4->SetNColumns(1);
	leg4->SetBorderSize(0);
	leg4->SetTextSize(0.03453388);
	leg4->SetFillStyle(0);

	TGraph *ufc12c12 = new TGraph(100);
	TGraph *ufc12c13 = new TGraph(100);
	TGraph *ufc13c13 = new TGraph(100);
	TGraph *ufc12c13_esw = new TGraph(100);
	TGraph *ufc13c13_esw = new TGraph(100);

	// define the reference system
	Double_t z1ref = 6, z2ref = 6, a1ref = 12, a2ref = 12;
	Double_t z1, z2, a1, a2;

	//henning's m3y+rep

	z1 = 6, z2 = 6, a1 = 12, a2 = 12;
	cal_uf(ufc12c12, z1, a1, z2, a2, "henning/henning_c12c12_cc_modified.dat", 1e-3, "c12c12/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");

	z1 = 6, z2 = 6, a1 = 12, a2 = 13;  
	cal_uf(ufc12c13, z1, a1, z2, a2, "henning/henning_c12c13_cc_modified.dat", 1e-3, "c12c13/cross.dat", z1ref, a1ref, z2ref, a2ref, "c12c12/cross.dat");

	z1 = 6, z2 = 6, a1 = 13, a2 = 13;
	cal_uf(ufc13c13,z1,a1,z2,a2,"henning/henning_c13c13_cc_modified.dat",1e-3,"c13c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");


	float linewith = 3.0;
	ufc12c12->SetLineColor(kRed);
	ufc12c12->SetLineWidth(linewith);	
	ufc12c13->SetLineColor(kGreen+2);
	ufc12c13->SetLineWidth(linewith);
	ufc13c13->SetLineColor(kBlue);
	ufc13c13->SetLineWidth(linewith);

	mg->Add(ufc12c12);
	mg->Add(ufc12c13);
	mg->Add(ufc13c13);

	return 0;
}
void cal_uf(TGraph* uf,Double_t z1, Double_t a1, Double_t z2, Double_t a2, char* f2, Double_t f2_factor, char* f1, 
		Double_t  z1ref, Double_t a1ref, Double_t z2ref, Double_t a2ref, char* f0) {
	using namespace std;
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
	for(int i=0;i<FF2->GetN();i++) {
		y[i]=ff(z1,a1,z2,a2,x[i],y[i]*f2_factor);
		x[i]=fx(z1,a1,z2,a2,x[i]);
	}

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

	Double_t emax=10;
	Int_t num = uf->GetN();
	Double_t de=(emax-emin)/num;

	x=uf->GetX();
	y=uf->GetY();

	TSpline *spline=0;
	for(int i=0;i<num;i++) {
		double fe=emin+de*i;
		Double_t fx_temp=fx(z1,a1,z2,a2,fe);
		x[i]=fx_temp;
		y[i]=FF2->Eval(fx_temp,spline,"S");
		y[i]=y[i]-FF1->Eval(fx_temp,spline,"S");
		y[i]=y[i]+FF0->Eval(fx_temp,spline,"S");
		// convert obtained UF to Smod for the reference system	
		x[i]=fx_temp*z1ref*z2ref;
		y[i]=ff_2_smod_ref(z1ref, a1ref, z2ref, a2ref, fx_temp, y[i]);
	}
}

Double_t fx(Double_t z1, Double_t a1, Double_t z2, Double_t a2,Double_t ecm) {
	Double_t temp=pow(ecm/z1/z2,1);
	return temp;
}
Double_t ff(Double_t z1, Double_t a1, Double_t z2, Double_t a2,Double_t ecm, Double_t sigma) {
	Double_t mu=a1*a2/(a1+a2);
	Double_t fk=sqrt(2*mu*931.5*ecm)/HBARC;
	// sigma in barn
	Double_t pene=sigma*1e2*fk**2/3.14;
	Double_t temp=log(pene);		
	Double_t temp=temp/sqrt(mu*z1*z2);
	return temp;
}

Double_t ff_2_smod_ref(Double_t z1ref, Double_t a1ref, Double_t z2ref, Double_t a2ref,
		Double_t fe_val, Double_t uf_val) {	
	Double_t ecm_ref=fe_val*z1ref*z2ref;
	Double_t mu_ref=a1ref*a2ref/(a1ref+a2ref);

	Double_t temp=uf_val;
	temp=temp*sqrt(mu_ref*z1ref*z2ref);
	// get log(sigma*Ecm)
	temp=temp-log(2*mu_ref*931.5/3.14/(HBARC*HBARC)*100.0);
	// convert to Smod
	temp=exp(temp+31.29*z1ref*z2ref*sqrt(mu_ref/(ecm_ref*1000))+0.46*ecm_ref);

	return temp;
}

void shadow(const TGraph* tge1,const TGraph* tge2,const TGraph* tge3, TMultiGraph *mg)
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
	gmin->SetLineWidth(1);
	gmin->SetLineStyle(2);
	gmax->SetFillStyle(3004);
	mg->Add(gmax, "CE3");
	mg->Add(gmin);
}
void shadow(const TGraph* tge1, const TGraph* tge2,const  TGraph* tge3, TMultiGraph *mg, double xl, double xh, double xcoeff)
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
	gmax->SetLineWidth(1);
	gmax->SetLineStyle(2);
	gmin->SetLineWidth(1);
	gmin->SetLineStyle(2);
	gmax->SetFillStyle(3002);
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
void fit()
{
}
