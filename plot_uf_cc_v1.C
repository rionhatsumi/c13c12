#define MARKERSIZE 1.5
#define LINEWIDTH 3
#define hbarc      197.326

void plot_uf_cc_v1() 
{
  gStyle->SetOptStat(0);

	Int_t n;
	TCanvas *c1 = new TCanvas("c1", "UFF", 100, 100, 1800, 1000);
	c1->SetGrid();
	c1->SetLeftMargin(0.07);
	c1->SetRightMargin(0.02);
	c1->SetTopMargin(0.05);
	TMultiGraph *mg = new TMultiGraph();
	TGraph *ufc12c12=new TGraph(100);
	TGraph *ufc12c13=new TGraph(100);
	TGraph *ufc13c13=new TGraph(100);

	TGraph *ufc12c13_esw=new TGraph(100);
	TGraph *ufc13c13_esw=new TGraph(100);

	// define the reference system
	Double_t z1ref, z2ref, a1ref, a2ref;
	z1ref=6, z2ref=6,a1ref=12,a2ref=12;

	Double_t z1, z2, a1, a2;

	Int_t theory_opt;
	cout<<"Theory models 1) M3Y+Rep; 2) TDHF : ";
	cin>>theory_opt;
	/*
	// test part
	z1=6, z2=6, a1=12, a2=12;
	cal_uf(ufc12c12,z1,a1,z2,a2,"c12c12/cross.dat",1,"c12c12/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	ufc12c12->Draw("AC");

	z1=6, z2=6, a1=12, a2=13;
	cal_uf(ufc12c13,z1,a1,z2,a2,"c12c13/cross.dat",1,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

	z1=6, z2=6, a1=13, a2=13;
	cal_uf(ufc13c13,z1,a1,z2,a2,"c13c13/cross.dat",1,"c13c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	*/

	if(theory_opt ==1 ) {

	//henning's m3y+rep
	z1=6, z2=6, a1=12, a2=12;
	cal_uf(ufc12c12,z1,a1,z2,a2,"henning/henning_c12c12_cc_modified.dat",1e-3,"c12c12/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

	z1=6, z2=6, a1=12, a2=13;  
	cal_uf(ufc12c13,z1,a1,z2,a2,"henning/henning_c12c13_cc_modified.dat",1e-3,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

	z1=6, z2=6, a1=13, a2=13;
	cal_uf(ufc13c13,z1,a1,z2,a2,"henning/henning_c13c13_cc_modified.dat",1e-3,"c13c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	} else if(theory_opt==2) {	
	// tdhf
	z1=6, z2=6, a1=12, a2=12;
	cal_uf(ufc12c12,z1,a1,z2,a2,"tdhf/x_12C_12C.out",1e-3,"c12c12/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

	z1=6, z2=6, a1=12, a2=13;
	cal_uf(ufc12c13,z1,a1,z2,a2,"tdhf/x_12C_13C.out",1e-3,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
//	cal_uf(ufc12c13_up,z1,a1,z2,a2,"tdhf/x_12C_13C.out",1.2e-3,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
//	cal_uf(ufc12c13_dn,z1,a1,z2,a2,"tdhf/x_12C_13C.out",0.8e-3,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

	z1=6, z2=6, a1=12, a2=14;
	cal_uf(ufc13c13,z1,a1,z2,a2,"tdhf/x_12C_14C.out",1e-3,"c13c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	}

	// esw models
	z1=6, z2=6, a1=12, a2=13;  
//	cal_uf(ufc12c13_esw,z1,a1,z2,a2,"esw/c12c13_esw_dayras.xsec",1,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
        cal_uf(ufc12c13_esw,z1,a1,z2,a2,"esw/c12c13_esw_imp_stock.xsec",1e-3,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	z1=6, z2=6, a1=13, a2=13;
	cal_uf(ufc13c13_esw,z1,a1,z2,a2,"esw/c13c13_1.2times.xsec",1e-3,"c13c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");

	// get the variation of the model prediction
	for(int i=0;i<ufc12c12->GetN();i++) {          
          Double_t tempx[3], tempy[3];
	  tempx[0]=0, tempx[1]=1, tempx[2]=2;
	  tempy[0]=*(ufc12c12->GetY()+i);
	  tempy[1]=*(ufc12c13->GetY()+i);
	  tempy[2]=*(ufc13c13->GetY()+i);
	  TGraph *gr_temp=new TGraph(3,tempx,tempy);
	  Double_t diff=gr_temp->GetMaximum()-gr_temp->GetMinimum();
	  cout<<*(ufc12c12->GetX()+i)<<" "
	      <<tempy[0]<<" "
	      <<tempy[1]<<" "
	      <<tempy[2]<<" "
	      <<diff<<endl;
	}

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
	mg->Draw("AC");

	//the up and down limits
	Double_t ymax[100], ymin[100],ydown[100];
	Double_t y0[100] = {0};
	Double_t *publicx, *yc12c12, *yc12c13, *yc13c13;
 	publicx = ufc12c12->GetX();
	yc12c12 = ufc12c12->GetY();
	yc12c13 = ufc12c13->GetY();
	yc13c13 = ufc13c13->GetY();
	for(int i = 0; i < 100; ++i)
	{
		ymax[i] = max(yc12c12[i], yc12c13[i], yc13c13[i]);
		ymin[i] = min(yc12c12[i], yc12c13[i], yc13c13[i]);
		ydown[i] = ymax[i] - ymin[i];
	}
	TGraphAsymmErrors *gmax = new TGraphAsymmErrors(100, publicx, ymax, y0, y0, ydown, y0);
	TGraph *gmin = new TGraph(100, publicx, ymin);
	gmax->SetLineWidth(LINEWIDTH);
	gmax->SetLineStyle(2);
	gmin->SetLineWidth(LINEWIDTH);
	gmin->SetLineStyle(2);
	gmax->SetFillStyle(3004);
	mg->Add(gmax,"CE3");
	
	mg->Add(gmin);
	mg->GetYaxis()->SetRangeUser(0 ,90e15);
	mg->GetXaxis()->SetTitle("E/z_{1}z_{2}");
	mg->GetXaxis()->SetTitleSize(0.05);
	mg->GetXaxis()->CenterTitle();
	mg->GetXaxis()->SetTitleOffset(0.9);
	mg->GetYaxis()->SetTitle("UFF");
	mg->GetYaxis()->SetTitleSize(0.05);
	mg->GetYaxis()->CenterTitle();
	mg->GetYaxis()->SetTitleOffset(0.6);
	
	TLegend *leg = new TLegend(0.80, 0.65, .98, .95);
//	leg->SetFillColor(10);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03053388);
	leg->SetFillStyle(0);
	leg->AddEntry(ufc12c12, "^{12}C+^{12}C", "l");
	leg->AddEntry(ufc12c13, "^{12}C+^{13}C", "l");
	leg->AddEntry(ufc13c13, "^{13}C+^{13}C", "l");
	leg->Draw();	
	// plot esw 
	ufc12c13_esw->SetLineColor(kRed);
	ufc12c13_esw->SetLineWidth(8);
	//	ufc12c13_esw->SetLineStyle(2);

	ufc13c13_esw->SetLineColor(kBlack);
	ufc13c13_esw->SetLineWidth(8);
	//	ufc13c13_esw->SetLineStyle(2);
	
	mg->Add(ufc12c13_esw);
	mg->Add(ufc13c13_esw);

//	ufc12c13_esw->Draw("C");
//	ufc13c13_esw->Draw("C");

	plot_exp_c12c12_data(leg);
	plot_exp_c12c13_data(leg);
	plot_exp_c13c13_data(leg);
}

void plot_exp_c12c12_data(TLegend *leg){
	// define the reference system
	Double_t z1ref, z2ref, a1ref, a2ref;
	z1ref=6, z2ref=6,a1ref=12,a2ref=12;

	TGraphErrors* gr_c12c12_bk = new TGraphErrors("exp_data/c12c12/becker.smod","%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c12_bk);
	TGraphErrors* gr_c12c12_sp = new TGraphErrors("exp_data/c12c12/c12c12_spillane_err_lessthan_0.4.smod","%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c12_sp);
	TGraphErrors* gr_c12c12_ag = new TGraphErrors("exp_data/c12c12/aguilera.smod","%lg %lg %lg");
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
	gr_c12c12_bk->SetLineColor(kRed);
	gr_c12c12_bk->SetMarkerColor(kRed);
	gr_c12c12_bk->SetMarkerStyle(24);
	gr_c12c12_bk->SetMarkerSize(MARKERSIZE);

	cal_uf_grerr(6,12,6,12,gr_c12c12_sp,1,"c12c12/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	gr_c12c12_sp->Draw("P");
	cal_uf_grerr(6,12,6,12,gr_c12c12_ag,1,"c12c12/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	gr_c12c12_ag->Draw("LP");
	cal_uf_grerr(6,12,6,12,gr_c12c12_bk,1,"c12c12/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	gr_c12c12_bk->Draw("LP");
	leg->AddEntry(gr_c12c12_sp, "sp ^{12}C+^{12}C", "p");
	leg->AddEntry(gr_c12c12_ag, "ag ^{12}C+^{12}C", "p");
	leg->AddEntry(gr_c12c12_bk, "bk ^{12}C+^{12}C", "p");
}
void plot_exp_c12c13_data(TLegend *leg) {
	//define the reference system
	Double_t z1ref, z2ref, a1ref, a2ref;
	z1ref=6, z2ref=6,a1ref=12,a2ref=12;

	TGraphErrors* gr_c12c13_st = new TGraphErrors("exp_data/c12c13/stockstad.smod","%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c13_st);
	TGraphErrors* gr_c12c13_imp = new TGraphErrors("exp_data/c12c13/imp_experr.smod","%lg %lg %lg");
	convert_smod_2_xsec(gr_c12c13_imp);

	gr_c12c13_st->SetLineWidth(2);
	gr_c12c13_st->SetLineColor(kBlue);
	gr_c12c13_st->SetMarkerColor(kBlue);
	gr_c12c13_st->SetMarkerSize(1);
	gr_c12c13_st->SetMarkerStyle(20);
	gr_c12c13_st->SetMarkerSize(MARKERSIZE);	
	
	gr_c12c13_imp->SetLineWidth(2);
	gr_c12c13_imp->SetLineColor(kRed);
	gr_c12c13_imp->SetMarkerColor(kRed);
	gr_c12c13_imp->SetMarkerStyle(20);
	gr_c12c13_imp->SetMarkerSize(MARKERSIZE);

	cal_uf_grerr(6,12,6,13,gr_c12c13_st,1,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	gr_c12c13_st->Draw("P");
	cal_uf_grerr(6,12,6,13,gr_c12c13_imp,1,"c12c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	gr_c12c13_imp->Draw("P");

}
void plot_exp_c13c13_data(TLegend *leg) {
//	define the reference system
		Double_t z1ref, z2ref, a1ref, a2ref;
	z1ref=6, z2ref=6,a1ref=12,a2ref=12;

	TGraphErrors* gr_c13c13_trent = new TGraphErrors("exp_data/c13c13/trentalange.smod","%lg %lg %lg");
	convert_smod_2_xsec(gr_c13c13_trent);

	gr_c13c13_trent->SetLineWidth(2);
	gr_c13c13_trent->SetLineColor(1);
	gr_c13c13_trent->SetMarkerColor(1);
	gr_c13c13_trent->SetMarkerStyle(22);
	gr_c13c13_trent->SetMarkerSize(MARKERSIZE);

	cal_uf_grerr(6,13,6,13,gr_c13c13_trent,1.2,"c13c13/cross.dat",z1ref,a1ref,z2ref,a2ref,"c12c12/cross.dat");
	gr_c13c13_trent->Draw("P");
}


// convert Smod into xsec in unit of barn
void convert_smod_2_xsec(TGraphErrors* gr1) {
	// Get the pointer for gr_1
	Double_t *x, *y, *ex, *ey;
	x=gr1->GetX();
	y=gr1->GetY();
	ex=gr1->GetEX();
	ey=gr1->GetEY();

	Int_t nn=gr1->GetN();

	for(int i=0;i<nn;i++) {
		Double_t tempy=87.21/sqrt(x[i])+0.46*x[i];
		y[i] = y[i]*exp(-tempy)/x[i];
		ey[i]=ey[i]*exp(-tempy)/x[i];
	}		
}

void cal_uf_grerr(Double_t z1, Double_t a1, Double_t z2, Double_t a2, TGraphErrors* uf, Double_t f2_factor, char* f1, 
		Double_t  z1ref, Double_t a1ref, Double_t z2ref, Double_t a2ref, char* f0) {
	TGraphErrors* FF2=uf;
	TGraph* FF1=new TGraph(f1);
	TGraph* FF0=new TGraph(f0);

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
	for(int i=0;i<FF2->GetN();i++) {
		y[i]=y[i]*f2_factor;
		ey[i]=ey[i]*f2_factor;

		if(y[i]<1e-15) {y[i]=1e-15; ey[i]=1e-15; continue;}

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
		y[i]=ff_2_smod_ref(z1ref, a1ref, z2ref, a2ref, 
				   uf_x, y[i]);		
		ey[i]=ff_2_smod_ref(z1ref, a1ref, z2ref, a2ref,
				    uf_x,uf_y+uf_yerr);
		ey[i]=fabs(ey[i]-y[i]);

		//		cout<<"before conversion : "<<uf_x[i]<<" "<<uf_y<<" "<<uf_yerr<<endl;
		//		cout<<"after  conversion : "<<x[i]<<" "<<y[i]<<" "<<ey[i]<<endl;
	}
}

void cal_uf(TGraph* uf,Double_t z1, Double_t a1, Double_t z2, Double_t a2, char* f2, Double_t f2_factor, char* f1, 
		Double_t  z1ref, Double_t a1ref, Double_t z2ref, Double_t a2ref, char* f0) {
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
	Double_t fk=sqrt(2*mu*931.5*ecm)/hbarc;
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
	temp=temp-log(2*mu_ref*931.5/3.14/(hbarc*hbarc)*100.0);
	// convert to Smod
	temp=exp(temp+31.29*z1ref*z2ref*sqrt(mu_ref/(ecm_ref*1000))+0.46*ecm_ref);

	return temp;
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
