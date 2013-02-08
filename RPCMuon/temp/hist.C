void hist()
{
  hist("legId != 0 && muonId_RPCMuLoose       == 1", "RPCMuLoose" );
  hist("legId != 0 && muonId_RPCMuMedium      == 1", "RPCMuMedium");
  hist("legId != 0 && muonId_RPCMuPromptTight == 1", "RPCMuPTight");

  hist("legId != 0 && muonId_TMOneStationLoose == 1", "TMOneStationLoose");
  hist("legId != 0 && muonId_TMOneStationTight == 1", "TMOneStationTight");
  hist("legId != 0 && muonId_TMTwoStationTest == 1", "TMTwoStationTest");

  //hist("legId != 0 && muonId_globalMuonMedium == 1", "globalMuonMedium");
  //hist("legId != 0 && muonId_globalMuonTight  == 1", "globalMuonTight" );
  //hist("legId != 0 && muonId_globalMuonPT     == 1", "globalMuonPTight");
  //hist("legId != 0 && muonId_globalMuonTest   == 1", "globalMuonTest");
}

void hist(const char* cut, const char* category)
{
  const int nbin = 100;
  const double xmin = 0.43, xmax = 0.56;

  TFile* fout = TFile::Open(Form("hist_%s_Eta_Barrel_M.root", category), "RECREATE");
  TH1F* hMassAll  = new TH1F("hMassAll" , "Vertex mass", nbin, xmin, xmax);
  TH1F* hMassPass = new TH1F("hMassPass", "Vertex mass", nbin, xmin, xmax);

  TH1F* hMassAll_Visible = new TH1F("hMassAll_Visible", "Vertex mass", nbin, xmin, xmax);
  TH1F* hMassAll_Barrel  = new TH1F("hMassAll_Barrel" , "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassAll_BarrelM  = new TH1F("hMassAll_BarrelM" , "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassAll_Overlap = new TH1F("hMassAll_Overlap", "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassAll_OverlapM = new TH1F("hMassAll_OverlapM", "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassAll_Endcap  = new TH1F("hMassAll_Endcap" , "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassAll_EndcapM  = new TH1F("hMassAll_EndcapM" , "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassAll_Forward = new TH1F("hMassAll_Forward", "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassAll_ForwardM = new TH1F("hMassAll_ForwardM", "Vertex mass", nbin, xmin, xmax);

  TH1F* hMassPass_Visible = new TH1F("hMassPass_Visible", "Vertex mass", nbin, xmin, xmax);
  TH1F* hMassPass_Barrel  = new TH1F("hMassPass_Barrel" , "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassPass_BarrelM  = new TH1F("hMassPass_BarrelM" , "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassPass_Overlap = new TH1F("hMassPass_Overlap", "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassPass_OverlapM = new TH1F("hMassPass_OverlapM", "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassPass_Endcap  = new TH1F("hMassPass_Endcap" , "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassPass_EndcapM  = new TH1F("hMassPass_EndcapM" , "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassPass_Forward = new TH1F("hMassPass_Forward", "Vertex mass", nbin, xmin, xmax);
  //TH1F* hMassPass_ForwardM = new TH1F("hMassPass_ForwardM", "Vertex mass", nbin, xmin, xmax);

  const char* cut_Visible  = "abs(track2.eta()) < 0.8                         && track2.pt() < 20";
  //const char* cut_Barrel   = "track1.eta() < 0.8      && track1.eta() > 0     && track1.pt() < 20";
  //const char* cut_BarrelM  = "track1.eta() > -0.8     && track1.eta() < 0     && track1.pt() < 20";
  //const char* cut_Overlap  = "track1.eta() >= 0.8     && track1.eta() < 1.2   && track1.pt() < 20";
  //const char* cut_OverlapM = "track1.eta() <= -0.8    && track1.eta() > -1.2  && track1.pt() < 20";
  //const char* cut_Endcap   = "track1.eta() >= 1.2     && track1.eta() < 1.6   && track1.pt() < 20";
  //const char* cut_EndcapM  = "track1.eta() <= -1.2    && track1.eta() > -1.6  && track1.pt() < 20";
  //const char* cut_Forward  = "track1.eta() >= 1.6                             && track1.pt() < 20";
  //const char* cut_ForwardM = "track1.eta() <= -1.6                            && track1.pt() < 20";

  const char* cut_Barrel  = "abs(track2.eta()) <= 0.8                             && track2.pt() < 20";
  //const char* cut_Overlap = "abs(track1.eta()) > 0.8  && abs(track1.eta()) <= 1.2 && track1.pt() < 20";
  //const char* cut_Endcap  = "abs(track1.eta()) > 1.2  && abs(track1.eta()) <= 1.6 && track1.pt() < 20";
  //const char* cut_Forward = "abs(track2.eta()) > 1.6  && abs(track2.eta()) <= 2.4 && track2.pt() < 20";

  TFile* f = TFile::Open("addCharge_20130111.root");
  TTree* tree = (TTree*) f->Get("fakeKshort/tree");

  fout->cd();
  tree->Draw("vertexMass >> hMassAll         ", "" , "goff");
  tree->Draw("vertexMass >> hMassPass        ", cut, "goff");

  tree->Draw("vertexMass >> hMassAll_Visible ", cut_Visible, "goff");
  tree->Draw("vertexMass >> hMassAll_Barrel  ", cut_Barrel , "goff");
  //tree->Draw("vertexMass >> hMassAll_BarrelM  ", cut_BarrelM , "goff");
  //tree->Draw("vertexMass >> hMassAll_Overlap ", cut_Overlap, "goff");
  //tree->Draw("vertexMass >> hMassAll_OverlapM ", cut_OverlapM, "goff");
  //tree->Draw("vertexMass >> hMassAll_Endcap  ", cut_Endcap , "goff");
  //tree->Draw("vertexMass >> hMassAll_EndcapM  ", cut_EndcapM , "goff");
  //tree->Draw("vertexMass >> hMassAll_Forward ", cut_Forward, "goff");
  //tree->Draw("vertexMass >> hMassAll_ForwardM ", cut_ForwardM, "goff");
  tree->Draw("vertexMass >> hMassPass_Visible", Form("%s && %s", cut, cut_Visible), "goff");
  tree->Draw("vertexMass >> hMassPass_Barrel ", Form("%s && %s", cut, cut_Barrel ), "goff");
  //tree->Draw("vertexMass >> hMassPass_BarrelM ", Form("%s && %s", cut, cut_BarrelM ), "goff");
  //tree->Draw("vertexMass >> hMassPass_Overlap", Form("%s && %s", cut, cut_Overlap), "goff");
  //tree->Draw("vertexMass >> hMassPass_OverlapM", Form("%s && %s", cut, cut_OverlapM), "goff");
  //tree->Draw("vertexMass >> hMassPass_Endcap ", Form("%s && %s", cut, cut_Endcap ), "goff");
  //tree->Draw("vertexMass >> hMassPass_EndcapM ", Form("%s && %s", cut, cut_EndcapM ), "goff");
  //tree->Draw("vertexMass >> hMassPass_Forward", Form("%s && %s", cut, cut_Forward), "goff");
  //tree->Draw("vertexMass >> hMassPass_ForwardM", Form("%s && %s", cut, cut_ForwardM), "goff");

  fout->Write();

}
