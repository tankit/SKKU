//root -l etaphi.C\(\"data_5k_RecHits_25cm_eta1_8.root\",\"data_5k_GlbHits_25cm_eta1_8.root\"\)
//root -l etaphi.C\(\"data_40k_RecHits_25cm_eta1_8.root\",\"mc_40k_RecHits_25cm_eta1_8.root\"\)
//root -l etaphi.C\(\"data_40k_RecHits_25cm_eta1_8.root\",\"mc_40k_GlbHits_test.root\"\)

//root -l etaphi.C\(\"data_20k_RecHits_test2.root\",\"mc_20k_RecHits_test2.root\"\)

//root -l etaphi.C\(\"data_40k_RecHits_update.root\",\"mc_40k_RecHits_update.root\"\) //RB<25cm and RE<20cm
//root -l etaphi.C\(\"data_40k_RecHits_10cm.root\",\"mc_40k_RecHits_10cm.root\"\)
//root -l etaphi.C\(\"data_60k_RecHits_10cm.root\",\"mc_60k_RecHits_10cm.root\"\)
//root -l etaphi.C\(\"data_60k_RecHits_10cm.root\",\"mc_60k_RecHits_10cm.root\"\,500)
//root -l etaphi.C\(\"data_60k_GlbHits_10cm.root\",\"mc_461k_GlbHits_10cm.root\",500\)
//root -l etaphi.C\(\"data_60k_RecHits_10cm.root\",\"data_60k_GlbHits_10cm.root\"\)
//root -l etaphi.C\(\"data_60k_RecHits_10cm.root\",\"mc_461k_RecHits_10cm.root\"\)
//root -l etaphi.C\(\"data_80k_RecHits_bx_merged.root\",\"mc_461k_RecHits_10cm.root\"\)
//root -l etaphi.C\(\"mc_461k_RecHits_10cm.root\",\"mc_461k_GlbHits_10cm.root\"\)
////root -l etaphi.C\(\"mc_new_RecHits_10cm.root\",\"mc_new_GlbHits_10cm.root\"\)

void etaphi(TString file1, TString file2, const int height=300, const int nhit=3, const float effmin=0.7)
{

  //gROOT->SetStyle("tdrStyle");

  //0 vs. 1: data vs. MC
  //0 vs. 1: RecHits vs. GlbHits

  TFile *f[2];
  TString fileName[2];
  fileName[0] = TString(file1);
  fileName[1] = TString(file2);
  
  string dataset[2];
  dataset[0] = "Data";
  dataset[1] = "MC";
  
  //gDirectory->cd("demo");
  
  const int l = 3;
  TLegend *leg[l];
  
  Int_t palette3[5] = {1,4,2,5,3};
  gStyle->SetPalette(5,palette3);
  
  TH2F *heffEtaPhi[2];
  TH1F *heffRPC[2][3];
  TH1F *heffEta[2], *heffPhi[2];
  TH1F *heta[2][4][7], *hphi[2][4][7], *hpt[2][4][7];
  
  int icol[7] = {1,2,3,1,4,5,6};
  int isty[7] = {24,24,22,20,23,21,25};
  
  for(int i=0; i<2; ++i) {
    
    f[i] = new TFile(fileName[i]);
    
    heffEtaPhi[i] = (TH2F*)f[i]->Get(Form("demo/heffEtaPhi%i",2));
    heffEtaPhi[i]->SetTitle(""); heffEtaPhi[i]->SetStats(0);
    heffEtaPhi[i]->GetXaxis()->SetTitle("muon track trajectory point #eta at RPC layer");
    heffEtaPhi[i]->GetYaxis()->SetTitle("muon track trajectory point #phi at RPC layer (rad)");
    
    heffEta[i] = (TH1F*)f[i]->Get(Form("demo/heffEta%i",2));
    heffEta[i]->SetTitle(""); heffEta[i]->SetStats(0);
    heffEta[i]->GetXaxis()->SetTitle("muon track trajectory point #eta at RPC layer");
    heffEta[i]->GetYaxis()->SetTitle("Acceptance of matching");
    heffEta[i]->GetYaxis()->SetTitleOffset(1.6); heffEta[i]->GetXaxis()->SetTitleOffset(1.2);
    
    heffPhi[i] = (TH1F*)f[i]->Get(Form("demo/heffPhi%i",2));
    heffPhi[i]->SetTitle(""); heffPhi[i]->SetStats(0);
    heffPhi[i]->GetXaxis()->SetTitle("muon track trajectory point #phi at RPC layer (rad)");
    heffPhi[i]->GetYaxis()->SetTitle("Acceptance of matching");
    heffPhi[i]->GetYaxis()->SetTitleOffset(1.6); heffPhi[i]->GetXaxis()->SetTitleOffset(1.2);
    
    for(int j=0; j<3; ++j) {
      heffRPC[i][j] = (TH1F*)f[i]->Get(Form("demo/heffRPC%i",j));
      heffRPC[i][j]->SetTitle(""); //heffRPC[i][j]->SetStats(0);
      if(j==0) heffRPC[i][j]->SetName(Form("heffRpc_%s",dataset[i].c_str()));
      else if(j==1) heffRPC[i][j]->SetName(Form("heffRB_%s",dataset[i].c_str()));
      else heffRPC[i][j]->SetName(Form("heffRE_%s",dataset[i].c_str()));
      heffRPC[i][j]->GetXaxis()->SetTitle("Barrel Acceptance");
      heffRPC[i][j]->GetYaxis()->SetTitle("Number of Rolls");
      
      if(i==0) heffRPC[i][j]->SetMarkerStyle(20);
      heffRPC[i][j]->GetYaxis()->SetTitleOffset(1.6); heffRPC[i][j]->GetXaxis()->SetTitleOffset(1.2);
      heffRPC[i][j]->SetMaximum(height);
    }
    
    for(int j=0; j<4; ++j) {
      for(int k=0; k<7; ++k) {
	heta[i][j][k] = (TH1F*)f[i]->Get(Form("demo/heta2_%i_%i",j,k));
	hphi[i][j][k] = (TH1F*)f[i]->Get(Form("demo/hphi2_%i_%i",j,k));
	hpt[i][j][k]  = (TH1F*)f[i]->Get(Form("demo/hpt2_%i_%i",j,k));
	
	heta[i][j][k]->SetStats(0); heta[i][j][k]->SetTitle(""); heta[i][j][k]->SetMaximum(1.1);
        heta[i][j][k]->SetMarkerColor(icol[k]); heta[i][j][k]->SetLineColor(icol[k]);
	heta[i][j][k]->SetMarkerStyle(isty[k]);
	heta[i][j][k]->GetXaxis()->SetTitle("muon #eta"); heta[i][j][k]->GetYaxis()->SetTitle("Fraction of muon tracks");
	
	float hmax = 1.1, hmin = 0;
	if(nhit==1) hmin = effmin;
	if(nhit==2) {
	  if(j<3) hmin = effmin; else if(j==3) hmin = 0.35;
	} else if(nhit==3) {
	  if(j==1) hmin = 0.7; else if(j==2) hmin = 0.35;
	} else if(nhit==4 && j==3) hmax = 0.3;

	hpt[i][j][k]->GetXaxis()->SetRangeUser(20,120);
	hpt[i][j][k]->SetStats(0); hpt[i][j][k]->SetTitle(""); hpt[i][j][k]->SetMaximum(hmax); hpt[i][j][k]->SetMinimum(hmin);
	hpt[i][j][k]->SetMarkerColor(icol[k]); hpt[i][j][k]->SetLineColor(icol[k]);
	hpt[i][j][k]->SetMarkerStyle(isty[k]); //hpt[i][j][k]->SetMarkerSize(0.8);
	hpt[i][j][k]->GetXaxis()->SetTitle("muon p_{T}"); hpt[i][j][k]->GetYaxis()->SetTitle("Fraction of muon tracks");
	//hpt[i][j][k]->GetYaxis()->SetTickLength(0.01); hpt[i][j][k]->GetYaxis()->SetTitleOffset(0.4); hpt[i][j][k]->GetXaxis()->SetTitleOffset(0.9);
	//hpt[i][j][k]->GetYaxis()->SetTitleSize(0.1); hpt[i][j][k]->GetXaxis()->SetTitleSize(0.05);

	hphi[i][j][k]->SetStats(0); hphi[i][j][k]->SetTitle(""); hphi[i][j][k]->SetMaximum(hmax); hphi[i][j][k]->SetMinimum(hmin);
	hphi[i][j][k]->SetMarkerColor(icol[k]); hphi[i][j][k]->SetLineColor(icol[k]);
	hphi[i][j][k]->SetMarkerStyle(isty[k]); //hphi[i][j][k]->SetMarkerSize(0.8);
	hphi[i][j][k]->GetXaxis()->SetTitle("muon #phi (rad)"); hphi[i][j][k]->GetYaxis()->SetTitle("Fraction of muon tracks");
	//hphi[i][j][k]->GetYaxis()->SetTickLength(0.01); hphi[i][j][k]->GetYaxis()->SetTitleOffset(0.4); hphi[i][j][k]->GetXaxis()->SetTitleOffset(0.9);
	//hphi[i][j][k]->GetYaxis()->SetTitleSize(0.1); hphi[i][j][k]->GetXaxis()->SetTitleSize(0.05);
	
	if(k==nhit) {
	  hphi[i][j][k]->SetMarkerStyle(20); hphi[i][j][k]->SetMarkerColor(1); hphi[i][j][k]->SetLineColor(1);
	  hpt[i][j][k]->SetMarkerStyle(20); hpt[i][j][k]->SetMarkerColor(1); hpt[i][j][k]->SetLineColor(1);
	}
	
      }
    }
  }
	
  for(int i=0; i<l; ++i) {
    //if(i==0) leg[i] = new TLegend(0.4,0.22,0.7,0.38,NULL,"brNDC");
    //if(i==0) leg[i] = new TLegend(0.75,0.9,1,1,NULL,"brNDC");
    if(i==0) leg[i] = new TLegend(0.7,0.8,0.85,0.9,NULL,"brNDC");
    else leg[i] = new TLegend(0.2,0.72,0.5,0.88,NULL,"brNDC");
    leg[i]->SetBorderSize(0);
    leg[i]->SetTextFont(62);
    leg[i]->SetTextSize(0.03);
    leg[i]->SetLineColor(0);
    leg[i]->SetLineStyle(1);
    leg[i]->SetLineWidth(0.5);
    leg[i]->SetFillColor(0); //if(i==0) leg[i]->SetFillColor(3);
    leg[i]->SetFillStyle(1001);
  }
  
  TCanvas *cvv = new TCanvas("cvv","cvv: Eta-Phi map",50,50,700,500);
  cvv->Divide(1,2);
  cvv->cd(1); heffEtaPhi[0]->Draw("zcol");
  cvv->cd(2); heffEtaPhi[1]->Draw("zcol");
  
  TCanvas *cvv2 = new TCanvas("cvv2","cvv2: Eta-Phi eff",75,75,1000,500);
  cvv2->Divide(2,1);
  cvv2->cd(1); heffEta[1]->Draw("HIST"); heffEta[0]->Draw("sames"); heffEta[1]->SetFillColor(3);
  heffEta[1]->SetMaximum(1.2);
  leg[0]->AddEntry(heffEta[0],"Data","lp"); leg[0]->AddEntry(heffEta[1],"MC","f");
  //leg[0]->AddEntry(heffEta[0],"Data allHits","lp"); leg[0]->AddEntry(heffEta[1],"Data GlbHits","f");  
  leg[0]->Draw();
  cvv2->cd(2); heffPhi[1]->Draw("HIST"); heffPhi[0]->Draw("sames"); heffPhi[1]->SetFillColor(3);
  heffPhi[1]->SetMaximum(1.2);
  leg[0]->Draw();
  
  //TCanvas *cvv3 = new TCanvas("cvv3","cvv3: Acceptance of matching",100,100,1000,500);
  TCanvas *cvv3 = new TCanvas("cvv3","cvv3: Acceptance of matching",100,100,500,500);
  //cvv3->Divide(2,1);
  //cvv3->cd(1); heffRPC[0]->Draw("");
  //cvv3->cd(2); heffRPC[1]->Draw("");
  //cvv3->Divide(2,2);
  //cvv3->cd(1); heffRPC[0][1]->Draw(""); //data-RB
  //cvv3->cd(2); heffRPC[0][2]->Draw(""); //data-RE
  //cvv3->cd(3); heffRPC[1][1]->Draw(""); //mc-RB
  //cvv3->cd(4); heffRPC[1][2]->Draw(""); //mc-RE
  cvv3->cd(1); heffRPC[1][1]->Draw(""); heffRPC[0][1]->Draw("esames"); heffRPC[1][1]->SetFillColor(3);
  gPad->SetTicks(); if(height==500) gPad->SetLogy();
  leg[1]->AddEntry(heffRPC[0][1],"Data","lp"); leg[1]->AddEntry(heffRPC[1][1],"MC","f");
  //leg[1]->AddEntry(heffRPC[0][1],"Data allHits","lp"); leg[1]->AddEntry(heffRPC[1][1],"Data GlbHits","f");
  leg[1]->Draw();
  //cvv3->cd(2); heffRPC[1][2]->Draw(""); heffRPC[0][2]->Draw("esames"); heffRPC[1][2]->SetFillColor(3); leg[1]->Draw();
	
  TCanvas *cvv4 = new TCanvas("cvv4","cvv4: Fraction of muon tracks vs. eta",125,125,700,500);
  cvv4->Divide(2,2);
  cvv4->cd(1); heta[0][0][1]->Draw(); heta[1][0][1]->Draw("HISTsame");
  cvv4->cd(2); heta[0][0][2]->Draw(); heta[1][0][2]->Draw("HISTsame");
  cvv4->cd(3); heta[0][0][3]->Draw(); heta[1][0][3]->Draw("HISTsame");
  cvv4->cd(4); heta[0][0][4]->Draw(); heta[1][0][4]->Draw("HISTsame");
  
  TCanvas *cvv5 = new TCanvas("cvv5","cvv5: Fraction of muon tracks vs. pt",150,150,500,700);
  cvv5->Divide(1,3);
  cvv5->cd(1); hpt[1][1][nhit]->Draw("HIST"); hpt[0][1][nhit]->Draw("esame"); hpt[1][1][nhit]->SetFillColor(3);
  cvv5->cd(2); hpt[1][2][nhit]->Draw("HIST"); hpt[0][2][nhit]->Draw("esame"); hpt[1][2][nhit]->SetFillColor(3);
  cvv5->cd(3); hpt[1][3][nhit]->Draw("HIST"); hpt[0][3][nhit]->Draw("esame"); hpt[1][3][nhit]->SetFillColor(3);
  
  TCanvas *cvv6 = new TCanvas("cvv6","cvv6: Fraction of muon tracks vs.phi",175,175,500,700);
  cvv6->Divide(1,3);
  cvv6->cd(1); hphi[1][1][nhit]->Draw("HIST"); hphi[0][1][nhit]->Draw("esame"); hphi[1][1][nhit]->SetFillColor(3);
  cvv6->cd(2); hphi[1][2][nhit]->Draw("HIST"); hphi[0][2][nhit]->Draw("esame"); hphi[1][2][nhit]->SetFillColor(3);
  cvv6->cd(3); hphi[1][3][nhit]->Draw("HIST"); hphi[0][3][nhit]->Draw("esame"); hphi[1][3][nhit]->SetFillColor(3);
  
  TCanvas *cvv7 = new TCanvas("cvv7","cvv7: Fraction of muon tracks vs. pt",200,200,700,500);
  cvv7->Divide(3,2);
  cvv7->cd(1); hpt[0][1][nhit]->Draw("");
  cvv7->cd(2); hpt[0][2][nhit]->Draw("");
  cvv7->cd(3); hpt[0][3][nhit]->Draw("");
  cvv7->cd(4); hpt[1][1][nhit]->Draw("");
  cvv7->cd(5); hpt[1][2][nhit]->Draw("");
  cvv7->cd(6); hpt[1][3][nhit]->Draw("");

  TCanvas *cvv8 = new TCanvas("cvv8","cvv8: Fraction of muon tracks vs. phi",225,225,700,500);
  cvv8->Divide(3,2);
  cvv8->cd(1); hphi[0][1][nhit]->Draw("");
  cvv8->cd(2); hphi[0][2][nhit]->Draw("");
  cvv8->cd(3); hphi[0][3][nhit]->Draw("");
  cvv8->cd(4); hphi[1][1][nhit]->Draw("");
  cvv8->cd(5); hphi[1][2][nhit]->Draw("");
  cvv8->cd(6); hphi[1][3][nhit]->Draw("");

  cvv2->Print("cvv2.pdf");

}
