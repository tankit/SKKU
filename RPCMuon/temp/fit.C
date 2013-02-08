void fit()
{
  TString outStr;
  
  //const char* category = "RPCMuLoose";
  //const char* category = "RPCMuMedium";
  const char* category = "RPCMuPTight";
  //const char* category = "TMOneStationLoose";
  //const char* category = "TMOneStationTight";
  //const char* category = "TMTwoStationTest";
  //const char* category = "globalMuonMedium";
  //const char* category = "globalMuonLoose";
  //const char* category = "globalMuonTight";
  //const char* category = "globalMuonPTight";

  //const double xbins[] = {-2.5, -1.6, -1.2, -0.8, 0, 0.8, 1.2, 1.6, 2.5};
  //const double xbins[] = { -1.6, -1.2, -0.8, 0, 0.8, 1.2, 1.6};
  //const char* xlabels[] = {"ForwardM", "EndcapM", "OverlapM", "BarrelM", "Barrel", "Overlap", "Endcap", "Forward"};
  //const char* xlabels[] = { "EndcapM", "OverlapM", "BarrelM", "Barrel", "Overlap", "Endcap"};

  const double xbins[] = {  0, 0.8 }; //,1.2, 1.6, 2.5};
  const char* xlabels[] = {  "Barrel" };//, "Overlap", "Endcap"};//, "Forward"};

  TGraphAsymmErrors* grp = new TGraphAsymmErrors();
  TGraphAsymmErrors* grpAll = new TGraphAsymmErrors();
  grp->SetTitle("Fake rate in eta;Pseudorapidity |#eta|;Fake rate (%)");
  grp->SetMinimum(0);
  grpAll->SetLineColor(kRed);
  grpAll->SetFillColor(kRed);
  grpAll->SetFillStyle(3004);
  for ( int i=0; i<1; ++i )
  {
    RooRealVar r = fitRatio(category, TString(Form("%s", xlabels[i])).Data());
    grp->SetPoint(i, (xbins[i]+xbins[i+1])/2, 100*r.getVal());
    const double dx = xbins[i+1]-xbins[i];
    grp->SetPointError(i, dx/2, dx/2, -100*r.getErrorLo(), 100*r.getErrorHi());
    outStr += Form("%s : Ratio = %.3g + %.3g - %.3g (%%)\n", xlabels[i], 100*r.getVal(), 100*r.getErrorHi(), 100*r.getErrorLo());
  }
  RooRealVar rAll = fitRatio(category, "Visible");
  grpAll->SetPoint(0, xbins[0], 100*rAll.getVal());
  //grpAll->SetPointError(0, 0, xbins[3]-xbins[0], -100*rAll.getErrorLo(), 100*rAll.getErrorHi());
  outStr += Form("All : Ratio = %.3g + %.3g - %.3g (%%)\n", 100*rAll.getVal(), 100*rAll.getErrorHi(), 100*rAll.getErrorLo());

  TCanvas* c = new TCanvas("c", "c", 500, 500);
  //c->SetLogy();
  grp->SetMaximum(1); grp->SetMinimum(1e-5);
  //grp->GetXaxis()->SetRangeUser(-1.6, 1.6);
  grp->GetXaxis()->SetRangeUser(-2.5, 2.5);
  grp->GetYaxis()->SetRangeUser(-0.01, 0.03);
  grp->Draw("AP");
  //grpAll->Draw("5");

  c->Print(Form("c_%s_eta_RPC_fakeRate.png", category));
  //c->Print(Form("c_%s_piP_fakeRate.png", category));
  //c->Print(Form("c_%s_piM_fakeRate.png", category));

  cout << outStr << endl;
}

RooRealVar fitRatio(const char* category, const char* region)
{
  using namespace RooFit;
  TFile* f = TFile::Open(Form("hist_%s_Eta_Barrel.root", category));
  //TFile* f = TFile::Open(Form("hist_%s_piP.root", category));
  //TFile* f = TFile::Open(Form("hist_%s_piM.root", category));

  TH1F* hMassA = (TH1F*)f->Get(Form("hMassAll_%s" , region));
  TH1F* hMassB = (TH1F*)f->Get(Form("hMassPass_%s", region));
  const double nA = hMassA->Integral();
  const double nB = hMassB->Integral();

  RooRealVar vertexMass("vertexMass", "Kshort mass", 0.43, 0.56, "GeV/c^{2}");
  RooDataHist hDataA("hDataA", "Kshort mass", RooArgList(vertexMass), hMassA);
  RooDataHist hDataB("hDataB", "Kshort mass", RooArgList(vertexMass), hMassB);

  RooWorkspace w("w");
  w.import(vertexMass);
  w.import(hDataA);
  w.import(hDataB);

  //w.factory("Voigtian::sigA(vertexMass, m0[0.498, 0.48, 0.51], w0[1e-2, 1e-5, 1e-1], sigmaA[1e-2, 1e-3, 1e+0])");
  //w.factory("Voigtian::sigB(vertexMass, m0, w0, sigmaB[1e-2, 1e-3, 1e+0])");
  
  w.factory("BreitWigner::bw(vertexMass, m0[0.498, 0.48, 0.55], w0[1e-2, 1e-5, 1e-1])");
  w.factory("CBShape::resA(vertexMass, mrA[0, -0.01, 0], sigmaA[5e-3, 1e-3, 1e-1], alphaA[1, 0, 3], nA[0, 0, 6])");
  w.factory("CBShape::resB(vertexMass, mrB[0, -0.01, 0], sigmaB[5e-3, 1e-3, 1e-1], alphaB[1, 0, 3], nB[0, 0, 6])");

  //w.factory("BreitWigner::bw(vertexMass, m0[0.498, 0.48, 0.55], w0[1e-3, 1e-5, 1e-1])");
  //w.factory("CBShape::resA(vertexMass, mrA[0, -0.01, 0], sigmaA[5e-3, 1e-3, 1e-1], alphaA[1, 0, 3], nA[0, 0, 6])");
  //w.factory("CBShape::resB(vertexMass, mrB[0, -0.01, 0], sigmaB[2e-3, 1e-3, 1e-1], alphaB[1, 0, 3], nB[0, 0, 6])");

  w.factory("FCONV::sigA(vertexMass, bw, resA)");
  w.factory("FCONV::sigB(vertexMass, bw, resB)");
  w.factory("Chebychev::bkgA(vertexMass, {p0A[0, -10, 10]})"); //, p1A[0, -10, 10]})");
  w.factory("Chebychev::bkgB(vertexMass, {p0B[0, -10, 10]})"); //, p1B[0, -10, 10]})");
  w.factory(Form("SUM::pdfA(nSigA[%.f, 0, %.f]*sigA, nBkgA[%.f, 0, %.f]*bkgA)", 0.5*nA, 1.1*nA, 0.5*nA, 1.1*nA));
  w.factory(Form("SUM::pdfB(nSigB[%.f, 0, %.f]*sigB, nBkgB[%.f, 0, %.f]*bkgB)", 0.5*nB, 1.1*nB, 0.5*nB, 1.1*nB));
  w.factory("weight[1, 0, 1e9]");
  w.factory("index[A,B]");
  w.exportToCint();

  RooSimultaneous simPdf("simPdf", "simPdf", w::index);
  simPdf.addPdf(w::pdfA, "A");
  simPdf.addPdf(w::pdfB, "B");

  RooDataSet dataA("dataA", "Kshort mass", RooArgList(w::vertexMass, w::weight), RooFit::WeightVar("weight"));
  RooDataSet dataB("dataB", "Kshort mass", RooArgList(w::vertexMass, w::weight), RooFit::WeightVar("weight"));

  for ( int i=0; i<hDataA.numEntries(); ++i )
  {
    dataA.add(*(hDataA.get(i)), hDataA.weight());
    dataB.add(*(hDataB.get(i)), hDataB.weight());
  }

  RooDataSet dataAB("dataAB", "Kshort mass", RooArgSet(vertexMass, w::weight), Index(w::index), Import("A", dataA), Import("B", dataB), WeightVar(w::weight));

  RooFitResult* result = simPdf.fitTo(dataAB, RooFit::Extended());
  w::m0.setConstant(); w::w0.setConstant();
  w::mrA.setConstant(); w::sigmaA.setConstant(); w::alphaA.setConstant(); w::nA.setConstant();
  w::mrB.setConstant(); w::sigmaB.setConstant(); w::alphaB.setConstant(); w::nB.setConstant();
  w::p0A.setConstant(); //w::p1A.setConstant();
  w::p0B.setConstant(); //w::p1B.setConstant();
  result = simPdf.fitTo(dataAB, RooFit::Save(), RooFit::Extended(), RooFit::Minos());
  
  RooPlot* frameA = vertexMass.frame();
  dataAB.plotOn(frameA, RooFit::Cut("index==index::A"));
  simPdf.plotOn(frameA, Slice(w::index, "A"), ProjWData(w::index, dataAB), RooFit::LineColor(kBlue));
  simPdf.plotOn(frameA, Slice(w::index, "A"), ProjWData(w::index, dataAB), RooFit::LineColor(kRed), RooFit::Components("bkgA"), RooFit::LineStyle(kDashed));

  RooPlot* frameB = vertexMass.frame();
  dataAB.plotOn(frameB, RooFit::Cut("index==index::B"));
  simPdf.plotOn(frameB, Slice(w::index, "B"), ProjWData(w::index, dataAB), RooFit::LineColor(kBlue));
  simPdf.plotOn(frameB, Slice(w::index, "B"), ProjWData(w::index, dataAB), RooFit::LineColor(kRed), RooFit::Components("bkgB"), RooFit::LineStyle(kDashed));

  TCanvas* cA = new TCanvas(Form("cA_%s", region), Form("c_%s", region), 500, 500);
  frameA->Draw();
  TCanvas* cB = new TCanvas(Form("cB_%s", region), Form("c_%s", region), 500, 500);
  frameB->Draw();

  result->Print("v");

  RooRealVar ratio("ratio", "ratio", -999);

  RooArgList floatPars = result->floatParsFinal();
  const RooRealVar* nSigA = floatPars.find("nSigA");
  const RooRealVar* nSigB = floatPars.find("nSigB");
  if ( !nSigA || !nSigB )
  {
    cout << "Cannot find nSigA or nSigB" << endl;
    return;
  }
  else
  {
    const double nA = nSigA->getVal();
    const double nB = nSigB->getVal();
    const double r = nB/nA;
    const double errHi = r*TMath::Hypot(nSigA->getErrorHi()/nA, nSigB->getErrorHi()/nB);
    const double errLo = r*TMath::Hypot(nSigA->getErrorLo()/nA, nSigB->getErrorLo()/nB);

    ratio.setVal(r);
    ratio.setAsymError(-errLo, errHi);
  }


  cA->Print(Form("cA_%s_eta_RPC_%s.png", category, region));
  cB->Print(Form("cB_%s_eta_RPC_%s.png", category, region));

  //cA->Print(Form("cA_%s_piP_%s.png", category, region));
  //cB->Print(Form("cB_%s_piP_%s.png", category, region));

  //cA->Print(Form("cA_%s_piM_%s.png", category, region));
  //cB->Print(Form("cB_%s_piM_%s.png", category, region));

  return ratio;
}

