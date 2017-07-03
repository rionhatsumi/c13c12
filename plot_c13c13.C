void plot_c13c13() {
	gROOT->Reset();
	gStyle->SetOptStat(0);
// below is c12c12 part
// exp data for 12C+12C
	TGraphErrors* gr_becker=new TGraphErrors("exp_data/c12c12/becker.smod","%lg %lg %lg");
	TGraphErrors* gr_c12c12ko=new TGraphErrors("../c12c12_kovar.smod","%lg %lg %lg");
	TGraphErrors* gr_c12c12pa=new TGraphErrors("exp_data/c12c12/aguilera.smod","%lg %lg %lg");
	// parital data from Spillane with rel_err<40%
	TGraphErrors* gr_c12c12sp=new TGraphErrors("exp_data/c12c12/spillane_err_lessthan_0.4.smod","%lg %lg %lg");
	
	convert_smod_2_s(gr_becker,12.0/2);
	convert_smod_2_s(gr_c12c12ko,12.0/2);
	convert_smod_2_s(gr_c12c12pa,12.0/2);
	convert_smod_2_s(gr_c12c12sp,12.0/2);	
	TF1* func_hind_c12c12 = new TF1("func_hind_c12c12",hindrance_s,0.5,5.0,5);
	func_hind_c12c12->SetParameters(-1.32,52.93,3.68,2.55e16,6.0);
	func_hind_c12c12->SetLineStyle(kDashed);
	func_hind_c12c12->SetLineColor(kRed);
	
	TMultiGraph* mg_c12c12_smod=new TMultiGraph;
	mg_c12c12_smod->Add(gr_becker);	
	mg_c12c12_smod->Add(gr_c12c12ko);
	mg_c12c12_smod->Add(gr_c12c12pa);
	mg_c12c12_smod->Add(gr_c12c12sp);
	
// below is the c12c13 part
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
	TMultiGraph* mg_c12c13_smod=new TMultiGraph;
	mg_c12c13_smod->Add(gr_c12c13dayras);
	mg_c12c13_smod->Add(gr_c12c13nd);
	mg_c12c13_smod->Add(gr_c12c13das);

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
	gr_c12c13das->SetMarkerStyle(33);

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
	gr_c13c13tren_smod->SetMarkerStyle(20);
	convert_smod_2_s(gr_c13c13tren_smod,13.0/2);
 
	Double_t xmin_global=1, xmax_global=6.5;
	TH2F *h2_sfactor=new TH2F("h2_sfactor","",100,xmin_global, xmax_global,1000,1e15,2e17);
	h2_sfactor->GetXaxis()->SetTitle("E_{c.m.} (MeV)");
	
	TMultiGraph* mg_c13c13_smod=new TMultiGraph;
	mg_c13c13_smod->Add(gr_c13c13tren_smod);
	
	Double_t xmin_fit=2.0, xmax_fit=5.0;
	TF1* func_hind_c13c13 = new TF1("func_hind_c13c13",hindrance_s,xmin_fit,xmax_fit,5);
	Fit_low_E(mg_c13c13_smod,func_hind_c13c13, 13.0/2, xmin_fit,xmax_fit);
	func_hind_c13c13->SetRange(0.5,xmax_fit);
	
	TCanvas *c1=new TCanvas("c1","",0,0,800,400);
	c1->SetLogy();
	h2_sfactor->Draw();
	
	mg_c12c12_smod->Draw("p");
	func_hind_c12c12->Draw("same");
	
	
	mg_c12c13_smod->Draw("p");
	xmax_fit=4.5;
	cout<<"Fitting c12+c13 data *****************"<<endl;
	Fit_low_E(mg_c12c13_smod,func_hind_c12c13, 12*13/25.0, xmin_fit,xmax_fit);
	func_hind_c12c13->SetRange(0.5,xmax_fit);
	func_hind_c12c13->Draw("same");
	
	mg_c13c13_smod->Draw("P");
	func_hind_c13c13->Draw("same");
	
}

void Fit_low_E(TMultiGraph* mg_temp, TF1* func, Double_t mu, Double_t xmin, Double_t xmax) {
// Creates a Root function based on function
   func->SetLineColor(kRed);
   func->SetLineStyle(kDashed);
// Sets initial values and parameter names
   func->SetParameters(-2.32,59.37,3.45,3e16,mu);
   func->SetParNames("A0","B0","Es","Smod_s","mu");
   func->SetRange(xmin,xmax);
   
   func->FixParameter(0,-2.32);
   func->FixParameter(1,59.37);
// Fit histogram in range defined by function   
   mg_temp->Fit(func->GetName(),"r0");	
   
// set limits for A0   
   Double_t xfactor=25;
   Double_t xmean=-2.32, xsigma=xfactor*0.24;
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
}

void convert_smod_2_s(TGraphErrors* gr1, Double_t mu) {
// Get the pointer for gr_1
	Double_t *x, *y, *ex, *ey;
	x=gr1->GetX();
	y=gr1->GetY();
	ex=gr1->GetEX();
	ey=gr1->GetEY();
	
	Int_t nn=gr1->GetN();
	Double_t factor=31.29*36/sqrt(1000.0/mu);
	
	for(int i=0;i<nn;i++) {
		Double_t tempy=(87.21-factor)/sqrt(x[i])+0.46*x[i];
//		cout<<x[i]<<" "<<tempy<<" "<<exp(-tempy)<<endl;
		y[i] = y[i]*exp(-tempy);
		ey[i]=ey[i]*exp(-tempy);
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
   Double_t sigma_s=smod_0/Es/exp(87.21/sqrt(Es)+0.46*Es)*1000;   
   
   Double_t f1=A0*(E-Es)-B0*(pow(Es/E,nn-1)-1)/pow(Es,nn-1)/(nn-1);
   Double_t f2=sigma_s*Es/E*exp(f1);
// convert back to standard s-factor   
   Double_t f3=f2/1000*E*exp(factor/sqrt(E)+0.0*E);
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
