{
  // In unnamed scripts, variables not forgotten at end, so must delete them before rerunning script, so ...
  gROOT->Reset("a");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat("");
  //gStyle->SetOptStat("emr");
  //  gStyle->SetOptStat("euom");
  gStyle->SetStatFontSize(0.035);
  gStyle->SetHistFillColor(kBlue);
  gStyle->SetHistFillStyle(1001);
  gStyle->SetMarkerSize(2.0);

  gStyle->SetStatFormat("5.3f");
  gStyle->SetStatFontSize(0.04);
  gStyle->SetOptFit(0111);
  gStyle->SetStatW(0.30);
  gStyle->SetStatH(0.02);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetTitleSize(0.05, "XYZ");

  gStyle->SetLabelSize(.04,"x");
  gStyle->SetLabelSize(.04,"y");

  gStyle->SetCanvasDefH(500);
  gStyle->SetCanvasDefW(800);

  TCanvas d1("d1");

  TFile *file[7];
  
  file[1] = new TFile("out_ttbar_ultimate_inv2RcorrOff_20180731_115437/Hist.root"); // Corr 1 off
  file[2] = new TFile("out_ttbar_ultimate_tiltedOff_20180725_141945/Hist.root"); // Corr 2 off
  file[3] = new TFile("out_ttbar_ultimate_helixExpOff_20180725_142024/Hist.root"); // Corr 3 off
  file[4] = new TFile("out_ttbar_ultimate_alpha0_20180725_142215/Hist.root"); // Corr 4 off
  file[5] = new TFile("out_ttbar_ultimate_20180725_141716/Hist.root");  // All on
  file[6] = new TFile("out_ttbar_ultimate_dodgy_20180725_143909/Hist.root"); // All off
  TLegend leg(0.7,0.15,0.9,0.45);    

  /*
  file[1] = new TFile("out_muon_ultimate_inv2RcorrOff_20180731_120259/Hist.root"); // Corr 1 off
  file[2] = new TFile("out_muonpt40_ultimate_tiltedOff_20180725_144435/Hist.root"); // Corr 2 off
  file[3] = new TFile("out_muonpt40_ultimate_helixExpOff_20180725_144513/Hist.root"); // Corr 3 off
  file[4] = new TFile("out_muonpt40_ultimate_alpha0_20180725_144629/Hist.root"); // Corr 4 off
  file[5] = new TFile("out_muonpt40_ultimate_20180725_144323//Hist.root");  // All on
  file[6] = new TFile("out_muonpt40_ultimate_dodgy_20180725_144726//Hist.root"); // All off
  TLegend leg(0.2,0.6,0.4,0.9);
  */

  TString name[7] = {"", "Corr. 1 off", "Corr. 2 off", "Corr. 3 off", "Corr. 4 off", "All on", "All off"};
  unsigned int icol[7] = {0, 1, 2, 3, 6, 8, 9};

  TH1F* his;
  TH2F* his2D;
  TProfile *prof[7];
  TEfficiency *teffi1, *teffi2, *teffi3;

  bool first = true;
  for (unsigned int i = 1; i <= 6; i++) {

    //file[i]->GetObject("TMTrackProducer/KF4ParamsComb/QoverPtResVsTrueEta_KF4ParamsComb", prof[i]);
    file[i]->GetObject("TMTrackProducer/KF4ParamsComb/QoverPtResVsTrueInvPt_KF4ParamsComb", prof[i]);
    //file[i]->GetObject("TMTrackProducer/KF4ParamsComb/Z0ResVsTrueEta_KF4ParamsComb", prof[i]);

    //prof[i]->SetMaximum(0.02);

    if (prof[i] == nullptr) {
      cout<<"ERROR: Input histogram missing "<<i<<endl;
      cin.get();
      continue;
    }

    prof[i]->SetMarkerStyle(20 + i);
    prof[i]->SetMarkerColor(icol[i]);
    prof[i]->SetMinimum(0.0);
    if (first) {
      first = false;
      prof[i]->Draw("P ");
    } else {
      prof[i]->Draw("P SAME");
    }
    leg.AddEntry(prof[i], name[i], "P");
    //leg.Draw();
    d1.Draw(); d1.Update(); 
  }
 
  //prof1->SetTitle(";1/Pt (1/GeV); #phi_{0} resolution");
  
  d1.Print("plot.pdf");
  cin.get(); 

  for (unsigned int i = 1; i <= 6; i++) {
    file[i]->Close();
    delete file[i];
  }
}
