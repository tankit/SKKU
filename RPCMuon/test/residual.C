//root -l residual.C\(\"data_mu_2GeV_5k_new1.root\"\)
//root -l residual.C\(\"data_mu_20GeV_5k_new1.root\"\) //hresRB with res_min
//root -l residual.C\(\"data_mu_20GeV_5k_new2.root\"\) //hresRB with recHitPosX_min
//root -l residual.C\(\"data_mu_20GeV_5k_new4.root\"\) //hresRB with the closest recHitPosX_min if trajectory crossed two rolls 
//root -l residual.C\(\"data_mu_20GeV_5k_new5.root\"\) //hresRB and hresRE

//root -l residual.C\(\"data_mu_20GeV_50k_new2.root\"\) //Nevt=19243 ~ 20k (duplication exists)

//root -l residual.C\(\"data_mu_20GeV_5k_new11.root\"\) //new
//root -l residual.C\(\"data_mu_20GeV_5k_new15.root\"\) //newer from Nov 7, 2011
//root -l residual.C\(\"data_mu_20GeV_5k_new20.root\"\) //newest from Nov 8, 2011
//root -l residual.C\(\"data_mu_20GeV_5k_new22.root\"\) //latest from Nov 9, 2011
//root -l residual.C\(\"data_mu_20GeV_5k_new23.root\",1\) //Eff_RE+- added on Nov 9, 2011
//root -l residual.C\(\"data_mu_20GeV_5k_new27.root\"\) //Eff_ring fixed on Nov 10, 2011
//root -l residual.C\(\"zroot5/data_mu_20GeV_5k_new36.root\"\)

//root -l residual.C\(\"data_40k_RecHits_update.root\"\)
//root -l residual.C\(\"mc_40k_RecHits_update.root\"\)

//root -l residual.C\(\"data_40k_RecHits_update.root\",\"mc_40k_RecHits_update.root\"\) //RB<25cm and RE<20cm
//root -l residual.C\(\"data_40k_RecHits_10cm.root\",\"mc_40k_RecHits_10cm.root\"\)
//root -l residual.C\(\"data_60k_RecHits_10cm.root\",\"mc_60k_RecHits_10cm.root\"\)
//root -l residual.C\(\"mc_60k_RecHits_10cm.root\",\"data_60k_RecHits_10cm.root\"\)
//root -l residual.C\(\"data_60k_RecHits_10cm.root\",\"mc_461k_RecHits_10cm.root\"\)
//root -l residual.C\(\"mc_461k_RecHits_10cm.root\",\"data_60k_RecHits_10cm.root\"\)
//root -l residual.C\(\"data_60k_RecHits_10cm.root\",\"data_60k_GlbHits_10cm.root\"\) //data: RecHits vs. GlbHits
//root -l residual.C\(\"data_60k_GlbHits_10cm.root\",\"mc_461k_GlbHits_10cm.root\"\) //GlbHits: data vs. MC
//root -l residual.C\(\"mc_461k_RecHits_10cm.root\",\"mc_461k_GlbHits_10cm.root\"\) //MC: RecHits vs. GlbHits
////root -l residual.C\(\"mc_new_RecHits_10cm.root\",\"mc_new_GlbHits_10cm.root\"\)

//root -l residual.C\(\"data_80k_RecHits_bx_merged.root\",\"mc_461k_RecHits_10cm.root\"\) //don't use eff because of merge

// Gaussian fit function
Double_t fitf(Double_t *v, Double_t *par)
{
	Double_t arg = 0;
	if (par[2] != 0) arg = (v[0] - par[1])/par[2];
	
	Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
	return fitval;
}

// Linear or Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
	
	Double_t triangle = par[0] + par[1]*x[0];
	if(x[0]<0) triangle = par[0] - par[1]*x[0];
	//cout << x[0] << " " << par[1] << endl;
	return triangle;
}

// Exponential
Double_t background2(Double_t *x, Double_t *par) {

	Double_t expo = exp(par[0] + par[1]*x[0]);
	if(x[0]<0) expo = exp(par[0] - par[1]*x[0]);
	return expo;
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
	return background2(x,par) * fitf(x,&par[3]);
}

//void residual(TString filename, bool REside=false)
void residual(TString file1, TString file2, bool REside=false)
{

  //gROOT->SetStyle("tdrStyle");

  TFile *f[2];
  TString fileName[2];
  fileName[0] = TString(file1);
  fileName[1] = TString(file2);

  string dataset[2];
  dataset[0] = "Data";
  dataset[1] = "MC";

  const int nL = 7;
  const int n = 5; //RB stations (MB1~MB4)
  const int m = 4; //RE disks    (ME1~ME3)
  const int l = 10;
  TLegend *leg[l];
  int icol[m] = {1,3,5,2};
  //int icol[m] = {1,2,3,4};
	
  float hmax = 1.2;

  TH1F *hresidualRB[2][n], *hresidualRBlayer[nL];
  TH1F *hresidualRE[2][m], *hresaddRE[m];
  TH1F *hresidualRing[6], *hresaddREring[7];
  TH1F *hlocalXRE[m], *hlocXaddRE[m];
  

  for(int i=0; i<2; ++i) {

    f[i] = new TFile(fileName[i]);

    for(int j=0; j<n; ++j) {
      hresidualRB[i][j] = (TH1F*)f[i]->Get(Form("demo/hresidualRB%i",j));
      hresidualRB[i][j]->Sumw2(); hresidualRB[i][j]->GetXaxis()->SetTitle("Residual (cm)"); hresidualRB[i][j]->SetMarkerStyle(20); hresidualRB[i][j]->GetXaxis()->SetTitleSize(0.06); 
      hresidualRB[i][j]->SetName(Form("hresidualRB%i_%s",j,dataset[i].c_str())); hresidualRB[i][j]->SetTitle("");
      if(i==1) {
	hresidualRB[i][j]->SetFillColor(3); // hresidualRB[i][j]->SetLineWidth(2);
      }      
      
      if(j<m) {
	hresidualRE[i][j] = (TH1F*)f[i]->Get(Form("demo/hresidualRE%i",j));
	hresidualRE[i][j]->Sumw2(); hresidualRE[i][j]->GetXaxis()->SetTitle("Residual (cm)"); hresidualRE[i][j]->SetMarkerStyle(20); hresidualRE[i][j]->GetXaxis()->SetTitleSize(0.06);
	hresidualRE[i][j]->SetName(Form("hresidualRE%i_%s",j,dataset[i].c_str())); //hresidualRE[i][j]->SetTitle(""); //error on Mac
	if(i==1) hresidualRE[i][j]->SetFillColor(3);
      }
    }



  }

  //TFile *fnew = new TFile(file1);
  //gDirectory->cd("demo");

  for(int i=0; i<6; ++i) {
    
    hresidualRing[i] = (TH1F*)f[0]->Get(Form("demo/hresidualRing%i",i));
    hresidualRing[i]->GetXaxis()->SetTitle("Residual (cm)");
    
  }
  
  for(int i=0; i<nL; ++i) {
    
    hresidualRBlayer[i] = (TH1F*)f[0]->Get(Form("demo/hresidualRBlayer%i",i)); //-200,200
    
    if(i!=0) {
      hresidualRBlayer[i]->SetTitle(Form("RB layer %i",i));
    }
    
    hresidualRBlayer[i]->GetXaxis()->SetTitle("Residual (cm)");

    if(i<n) {

      //hresidualRB[i] = (TH1F*)f[0]->Get(Form("demo/hresidualRB%i",i));
      //if(i!=0) hresidualRB[i]->SetTitle(Form("RB station %i",i));
      //hresidualRB[i]->GetXaxis()->SetTitle("Residual (cm)"); hresidualRB[i]->SetMarkerStyle(20);

    }

    if(i<m) {
      //hresidualRE[i] = (TH1F*)f[0]->Get(Form("demo/hresidualRE%i",i));
      //hlocalXRE[i] = (TH1F*)f[0]>Get(Form("demo/hlocalXRE_%i",i));
      hlocalXRE[i] = (TH1F*)f[0]->Get(Form("demo/hlocalXRE1_%i",i));

      if(i!=0) {
	//hresidualRE[i]->SetTitle(Form("RE disk %i",i));
	hlocalXRE[i]->SetTitle(Form("RE disk %i",i));
      }

      //hresidualRE[i]->GetXaxis()->SetTitle("Residual (cm)"); hresidualRE[i]->SetMarkerStyle(20);
      hlocalXRE[i]->GetXaxis()->SetTitle("localX (cm)");

    }

  }

  for(int i=0; i<l; ++i) {
    //if(i>4) leg[i] = new TLegend(0.62,0.67,0.75,0.82,NULL,"brNDC"); 
    if(i>4) leg[i] = new TLegend(0.72,0.67,0.85,0.82,NULL,"brNDC");
    //else leg[i] = new TLegend(0.72,0.67,0.85,0.82,NULL,"brNDC");
	else leg[i] = new TLegend(0.15,0.67,0.4,0.82,NULL,"brNDC");
    //else leg[i] = new TLegend(0.72,0.67,0.88,0.75,NULL,"brNDC");
    leg[i]->SetBorderSize(0);
    leg[i]->SetTextFont(62);
    leg[i]->SetTextSize(0.03);
    leg[i]->SetLineColor(0);
    leg[i]->SetLineStyle(1);
    leg[i]->SetLineWidth(0.5);
    leg[i]->SetFillColor(0);
    leg[i]->SetFillStyle(1001);
  }

  TCanvas *cvv = new TCanvas("cvv","cvv: Residual in RB and RE",25,25,1000,700);
  cvv->Divide(1,2);
  //cvv->cd(1); hresidualRB[0]->Draw(); gPad->SetLogy();
  //cvv->cd(2); hresidualRE[0]->Draw(); gPad->SetLogy();
  //cvv->Divide(1,2);
  //cvv->cd(1); hresidualRB[0]->Fit("gaus"); hresidualRB[0]->Draw("e");
  //cvv->cd(2); hresidualRE[0]->Fit("gaus"); hresidualRE[0]->Draw("e");
  int normRB = hresidualRB[0][0]->Integral(), normRE = hresidualRE[0][0]->Integral();
  //cvv->cd(1); hresidualRB[0][0]->Draw("e");

  cvv->cd(1); //hresidualRB[0][0]->Draw("e");
  //hresidualRB[0][0]->DrawNormalized("e",normRB); hresidualRB[1][0]->DrawNormalized("HISTsames",normRB);
  hresidualRB[1][0]->DrawNormalized("HIST",normRB); hresidualRB[0][0]->DrawNormalized("esames",normRB); 
  gPad->SetLogy();
  leg[9]->AddEntry(hresidualRB[0][0],"Data","lp");
  leg[9]->AddEntry(hresidualRB[1][0],"MC","lf");
  //leg[9]->AddEntry(hresidualRB[0][0],"Data allHits","lp");
  //leg[9]->AddEntry(hresidualRB[1][0],"Data GlbHits","lf");
  leg[9]->Draw();
	/*
	// create a TF1 with the range from 0 to 3 and 6 parameters
	TF1 *fitFcn = new TF1("fitFcn",fitFunction,-50,50,6);
	//TF1 *fitFcn = new TF1("fitFcn",fitf,-50,50,3);
	//TF1 *fitFcn1 = new TF1("fitFcn1","pol1",-50,-10);
	//TF1 *fitFcn2 = new TF1("fitFcn2","gaus",-50,50);
	//TF1 *fitFcn3 = new TF1("fitFcn3","pol1",10,50);
	//TF1 *total = new TF1("total","pol1(0)+gaus(2)+pol1(5)",-50,50);
	TF1 *fitFcn1 = new TF1("fitFcn1","expo",-50,-10);
	TF1 *fitFcn2 = new TF1("fitFcn2","gaus",-10,10);
	TF1 *fitFcn3 = new TF1("fitFcn3","expo",10,50);
	TF1 *total = new TF1("total","expo(0)+gaus(2)",-50,50);

	fitFcn->SetNpx(500);
	fitFcn->SetLineWidth(4);
	fitFcn->SetLineColor(kMagenta);
	
	// first try without starting values for the parameters
	// This defaults to 1 for each param.
	// this results in an ok fit for the polynomial function
	// however the non-linear part (lorenzian) does not
	// respond well.
	fitFcn->SetParameters(1,1,1,1,1,1);
	//fitFcn->SetParameters(1,1,1);
	//hresidualRB[0][0]->Fit("fitFcn","");
	
	fitFcn1->SetLineColor(kMagenta); fitFcn3->SetLineColor(kMagenta);
	fitFcn2->SetLineColor(kBlue);
	hresidualRB[0][0]->Fit("fitFcn1","R");
	hresidualRB[0][0]->Fit("fitFcn2","R+");
	hresidualRB[0][0]->Fit("fitFcn3","R+");
	Double_t par[7];
	fitFcn1->GetParameters(&par[0]);
	fitFcn2->GetParameters(&par[2]);
	fitFcn3->GetParameters(&par[5]);
	total->SetParameters(par);
	hresidualRB[0][0]->Fit(total,"RWLFBM+");
	*/
	
  cvv->cd(2);
  //hresidualRE[0][0]->DrawNormalized("e",normRE); hresidualRE[1][0]->DrawNormalized("HISTsames",normRE);
  hresidualRE[1][0]->DrawNormalized("HIST",normRE); hresidualRE[0][0]->DrawNormalized("esames",normRE);
  gPad->SetLogy();

  leg[9]->Draw();

  //gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1);

  int norm[4];
  TF1 *fit[2][4];
  TCanvas *cvv2 = new TCanvas("cvv2","cvv2: Residual in stations (RB)",100,100,700,500);
  cvv2->Divide(2,2);
  for(int i=0; i<2; ++i) {
    for(int j=1; j<n; ++j) {
      cvv2->cd(j); gPad->SetLogy();
      if(i==0) {
	norm[j] = hresidualRB[i][j]->Integral(); hresidualRB[i][j]->DrawNormalized("e",norm[j]);
	fit[i][j] = new TF1(Form("fit%i_%i",i,j),"gaus",-50,50); fit[i][j]->SetLineWidth(2);
	fit[i][j]->SetParameters(norm[j],hresidualRB[i][j]->GetMean(),hresidualRB[i][j]->GetRMS());
	////fit[i][j]->SetParNames("Constant","Mean_value","Sigma");
	hresidualRB[i][j]->Fit(fit[i][j],"");
        leg[9]->Draw();
	    //if(j==1) {
		//} else hresidualRB[i][j]->Fit(fit[i][j],"");
		  
      } else {
	hresidualRB[i][j]->DrawNormalized("HISTsame",norm[j]);
	hresidualRB[i-1][j]->DrawNormalized("esame",norm[j]); //fit[i][j]->Draw("sames");
      }
      //hresidualRB[i][j]->Fit("gaus"); hresidualRB[i][j]->Draw("e"); //hresidualRB[i][j]->GetXaxis()->SetRangeUser(-10,10);
      //hresidualRB[i][j]->GetXaxis()->SetRangeUser(-15,15); //hresidualRB[i][j]->GetXaxis()->SetRangeUser(-20,20);
    }
  }

  TCanvas *cvv3 = new TCanvas("cvv3","cvv3: Residual in layers (RB)",125,125,700,500);
  cvv3->Divide(2,3);
  for(int i=1; i<nL; ++i) {
    cvv3->cd(i); gPad->SetLogy();
    hresidualRBlayer[i]->DrawCopy("");
  }
 
  TCanvas *cvv4 = new TCanvas("cvv4","cvv4: Residual in stations (RB)",150,150,700,500);
  cvv4->Divide(2,2);
  for(int i=1; i<n; ++i) {
    cvv4->cd(i); gPad->SetLogy();
    hresidualRB[0][i]->Draw("HIST"); if(i<3) hresidualRBlayer[i]->Draw(); else hresidualRBlayer[i+2]->Draw();
    if(i==1) {
      hresidualRBlayer[2]->Draw("same"); hresidualRBlayer[2]->SetFillColor(3); hresidualRBlayer[2]->SetLineColor(3); //hresidualRB[i]->Draw("same");
      hresidualRBlayer[2]->Draw(""); hresidualRBlayer[1]->Draw("same");
      //hresidualRB[0][i]->Draw("HISTsame");
	  leg[i]->AddEntry(hresidualRB[0][i],Form("RB layer %i",1),"f");
      leg[i]->AddEntry(hresidualRBlayer[2],Form("RB layer %i",2),"f");
      leg[i]->Draw();
    }
    else if(i==2) {
      hresidualRBlayer[4]->Draw("same"); hresidualRBlayer[4]->SetFillColor(3); hresidualRBlayer[4]->SetLineColor(3); //hresidualRB[i]->Draw("same");
      hresidualRBlayer[4]->Draw(""); hresidualRBlayer[3]->Draw("same");
      //hresidualRB[0][i]->Draw("HISTsame");
	  leg[i]->AddEntry(hresidualRB[0][i],Form("RB layer %i",3),"f");
      leg[i]->AddEntry(hresidualRBlayer[4],Form("RB layer %i",4),"f");
      leg[i]->Draw();
    }
    else { leg[i]->AddEntry(hresidualRB[0][i],Form("RB layer %i",i+2),"f"); leg[i]->Draw(); }

    gPad->RedrawAxis();

  }

  TCanvas *cvv5 = new TCanvas("cvv5","cvv5: RE disks",175,175,700,500);
  cvv5->Divide(1,3);
  for(int i=1; i<m; ++i) {
    cvv5->cd(i); gPad->SetLogy(); hresidualRE[0][i]->DrawCopy();
    //hresidualRE[0][i]->Fit("gaus"); hresidualRE[0][i]->DrawCopy("e");
    //hlocalXRE[i]->DrawCopy("");

    hresaddRE[i] = (TH1F*)hresidualRE[1][i]->Clone(Form("hresaddRE%i",i)); hresaddRE[i]->SetFillColor(icol[i]); //hresaddRE[i]->SetLineColor(icol[i]);
    hresaddRE[i]->SetStats(0); hresaddRE[i]->SetTitle(""); //hresaddRE[i]->SetTitle("Residuals for RE disks");
    if(i>1) hresaddRE[i]->Add(hresaddRE[i-1]);

    hlocXaddRE[i] = (TH1F*)hlocalXRE[i]->Clone(Form("hlocXaddRE%i",i)); hlocXaddRE[i]->SetFillColor(icol[i]); //hlocXaddRE[i]->SetLineColor(icol[i]);
    hlocXaddRE[i]->SetStats(0); hlocXaddRE[i]->SetTitle("recHit localX for RE disks"); 
    if(i>1) hlocXaddRE[i]->Add(hlocXaddRE[i-1]);

  }

  TCanvas *cvv6 = new TCanvas("cvv6","cvv6: RE disks (cumulative)",200,200,500,500);
  cvv6->Divide(1,3);
  for(int i=m-1; i>0; --i) {
    //cvv6->cd(i); hresaddRE[i]->Draw();
    if(i==m-1) hresaddRE[i]->DrawCopy("HIST");
    else hresaddRE[i]->DrawCopy("HISTsame");
    gPad->RedrawAxis(); gPad->SetLogy();

    leg[0]->AddEntry(hresaddRE[i],Form("RE disk %i",i),"f");
  }
  leg[0]->Draw();


  /*
  //RE ring -2
  hresaddREring[0] = (TH1F*)hresidualRing[1]->Clone(Form("hresaddREring%i",0));
  hresaddREring[0]->Add(hresidualRing[2]); hresaddREring[0]->Add(hresidualRing[3]);
  hresaddREring[0]->SetTitle("Residual for RE ring -2"); hresaddREring[0]->SetName("hresidualRing-2");
  //RE ring +2
  hresaddREring[1] = (TH1F*)hresidualRing[4]->Clone(Form("hresaddREring%i",1));
  hresaddREring[1]->Add(hresidualRing[5]); hresaddREring[1]->Add(hresidualRing[6]);
  hresaddREring[1]->SetTitle("Residual for RE ring 2"); hresaddREring[1]->SetName("hresidualRing+2");
  //RE ring -3
  hresaddREring[2] = (TH1F*)hresidualRing[7]->Clone(Form("hresaddREring%i",2));
  hresaddREring[2]->Add(hresidualRing[8]); hresaddREring[2]->Add(hresidualRing[9]);
  hresaddREring[2]->SetTitle("Residual for RE ring -3"); hresaddREring[2]->SetName("hresidualRing-3");
  //RE ring +3
  hresaddREring[3] = (TH1F*)hresidualRing[10]->Clone(Form("hresaddREring%i",3));
  hresaddREring[3]->Add(hresidualRing[11]); hresaddREring[3]->Add(hresidualRing[12]);
  hresaddREring[3]->SetTitle("Residual for RE ring 3"); hresaddREring[3]->SetName("hresidualRing+3");

  //RE ring -2,+2
  hresaddREring[4] = (TH1F*)hresaddREring[0]->Clone(Form("hresaddREring%i",4)); hresaddREring[4]->SetFillColor(3); hresaddREring[4]->SetLineColor(3);
  hresaddREring[4]->Add(hresaddREring[1]);;
  hresaddREring[4]->SetTitle("Residual for RE rings 2"); hresaddREring[4]->SetName("hresidualRing2");
  //RE ring -3,+3
  hresaddREring[5] = (TH1F*)hresaddREring[2]->Clone(Form("hresaddREring%i",5)); hresaddREring[5]->SetFillColor(3);
  hresaddREring[5]->Add(hresaddREring[3]);;
  hresaddREring[5]->SetTitle("Residual for RE rings 3"); hresaddREring[5]->SetName("hresidualRing3");

  //RE ring cumulative
  hresaddREring[6] = (TH1F*)hresaddREring[4]->Clone(Form("hresaddREring%i",6)); hresaddREring[6]->SetFillColor(0); hresaddREring[6]->SetLineColor(1);
  hresaddREring[6]->Add(hresaddREring[5]);
  hresaddREring[6]->SetTitle("Residuals for RE rings");

  for(int k=4; k<7; ++k) hresaddREring[k]->SetTitle("");

  TCanvas *cvv7 = new TCanvas("cvv7","cvv7: Residual in rings (RE)",225,225,700,500);
  cvv7->Divide(2,2);
  cvv7->cd(1); gPad->SetLogy(); hresaddREring[0]->DrawCopy(""); //RE ring -2
  cvv7->cd(2); gPad->SetLogy(); hresaddREring[1]->DrawCopy(""); //RE ring +2
  cvv7->cd(3); gPad->SetLogy(); hresaddREring[2]->DrawCopy(""); //RE ring -3
  cvv7->cd(4); gPad->SetLogy(); hresaddREring[3]->DrawCopy(""); //RE ring +3
  

  TCanvas *cvv8 = new TCanvas("cvv8","cvv8: Residual in rings (RE)",250,250,700,500);
  cvv8->Divide(1,2);
  cvv8->cd(1); gPad->SetLogy(); hresaddREring[1]->SetMarkerStyle(20);
  hresaddREring[1]->Draw("e"); hresaddREring[0]->Draw("sames"); hresaddREring[1]->SetTitle("");
  leg[5]->AddEntry(hresaddREring[1],Form("RE ring %s","+2"),"lp");
  leg[5]->AddEntry(hresaddREring[0],Form("RE ring %s","-2"),"f");
  leg[5]->Draw();
  cvv8->cd(2); gPad->SetLogy(); hresaddREring[3]->SetMarkerStyle(20);
  hresaddREring[3]->Draw("e"); hresaddREring[2]->Draw("sames"); hresaddREring[3]->SetTitle("");
  leg[6]->AddEntry(hresaddREring[3],Form("RE ring %s","+3"),"lp");
  leg[6]->AddEntry(hresaddREring[2],Form("RE ring %s","-3"),"f");
  leg[6]->Draw();

  TCanvas *cvv9 = new TCanvas("cvv9","cvv9: RE rings (cumulative)",275,275,700,500);
  hresaddREring[6]->Draw(); gPad->SetLogy(); hresaddREring[6]->SetStats(0);
  hresaddREring[4]->DrawCopy("same"); //hresaddREring[6]->DrawCopy("same");
  gPad->RedrawAxis();
  if(REside) {
    leg[7]->AddEntry(hresaddREring[6],Form("RE+"),"f");
    leg[7]->AddEntry(hresaddREring[4],Form("RE-"),"f");
  } else {
    leg[7]->AddEntry(hresaddREring[6],Form("RE ring %i",3),"f");
    leg[7]->AddEntry(hresaddREring[4],Form("RE ring %i",2),"f");
  }
  leg[7]->Draw();

  TCanvas *cvv10 = new TCanvas("cvv10","cvv10: RE rings",300,300,700,500);
  hresaddREring[4]->Draw(""); hresaddREring[5]->Draw("esames"); hresaddREring[5]->SetMarkerStyle(20); hresaddREring[4]->SetLineColor(1); gPad->SetLogy();
  if(REside) {
    leg[8]->AddEntry(hresaddREring[5],Form("RE+"),"lp"); hresaddREring[5]->SetName("hresidualRE+");
    leg[8]->AddEntry(hresaddREring[4],Form("RE-"),"f"); hresaddREring[4]->SetName("hresidualRE-");
  } else {
    leg[8]->AddEntry(hresaddREring[5],Form("RE ring %i",3),"lp");
    leg[8]->AddEntry(hresaddREring[4],Form("RE ring %i",2),"f");
  }
  leg[8]->Draw();
  */

  //TCanvas *cvv11 = new TCanvas("cvv11","cvv11: RE rings (comparison)",325,325,700,500);
  //hresidualRE[0][0]->Draw();
  //hresaddREring[6]->Draw("esames"); hresaddREring[6]->SetMarkerStyle(20);

}
