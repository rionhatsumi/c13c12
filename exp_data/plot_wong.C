{
  gROOT->Reset();
  // 3.14 * RB ^ 2 * omega / E * Log(1 + Exp((E - VB) / omega))
  TF1 *f1 = new TF1("f1","10.0*[0]^ 2*[1]/(2*x)*log(1+exp(2*3.1415926*(x-[2])/[1]))/1000.0*x*exp(87.21/sqrt(x)+0.46*x)",4.9,8.0);
  f1->SetParName(0,"rb");
  f1->SetParName(1,"omega");
  f1->SetParName(2,"vb");
  f1->SetLineColor(2);

  f1->SetParameters(7.0,2.1,5.54);

  TGraphErrors gr_c12c13("c12c13_dayras.dat","%lg %lg %lg");
  gr_c12c13->SetMarkerStyle(20);
  gr_c12c13->SetMarkerColor(6);
  Double_t* ey = gr_c12c13->GetEY();
  Double_t* y = gr_c12c13->GetY();
  for(int i=0;i<gr_c12c13->GetN();i++) {
    ey[i]=sqrt(pow(ey[i],2)+pow(y[i]*0.3,2));
  }
  gr_c12c13->Draw("AP");
  gr_c12c13->Fit(f1,"RE");
  cout<<"12C+13C finished."<<endl<<endl;

  //  f1->Draw();
  TGraphErrors gr_c13c13("c13c13_trentalange.dat","%lg %lg %lg");
  gr_c13c13->SetMarkerStyle(22);
  Double_t* ey = gr_c13c13->GetEY();
  Double_t* y = gr_c13c13->GetY();
  for(int i=0;i<gr_c13c13->GetN();i++) {
    ey[i]=sqrt(pow(ey[i],2)+pow(y[i]*0.3,2));
  }
  gr_c13c13->Draw("P");
  TF1 *f2=(TF1*) f1->Clone("f2");
  gr_c13c13->Fit(f2,"RE");  
  cout<<"13C+13C finished."<<endl<<endl;

//  TGraphErrors gr_c12c12("becker_kovar_peak.dat","%lg %lg %lg");
  TGraphErrors gr_c12c12("becker_smod.dat","%lg %lg %lg");
  gr_c12c12->SetMarkerStyle(21);
  gr_c12c12->SetMarkerColor(2);
  gr_c12c12->Draw("P");
  TF1 *f3=(TF1*) f1->Clone("f3");
  f3->SetParLimits(1,0.3,10);
  f3->SetRange(4.9,15);
  gr_c12c12->Fit(f3,"RE");  
  cout<<"12C+12C finished."<<endl<<endl;
}
