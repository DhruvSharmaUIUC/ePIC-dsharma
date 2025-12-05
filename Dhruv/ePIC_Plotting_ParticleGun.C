// Run this macro, from the Linux / Terminal command line, as
// root -l -q ePIC_Plotting.C
// The only condition is that previously ran ePIC_Analysis.C  with the identical "strang" as hardcoded in this macro (below)
// Strang must match the string in a file starting with out*, in this directory
// CR 2024-08-14/20

void ePIC_Plotting_ParticleGun(int gevnum = 10)
{
  gSystem->Exec("date");

  // define nHCal acceptance (2024-07-16) - make sure this is consistent with *_Analysis.C:
    int nlayers_nhcal = 10;
  //
  
  // Define name of MC file:
    //TString strang = Form("Muon_%d_GeV", gevnum);
    TString strang = Form("Muon_%d_GeV_Uniform", gevnum);
    //TString filetype_dir = Form("ParticleGun/Muon_%d_GeV/", gevnum);
    TString filetype_dir = Form("ParticleGun/Muon_%d_GeV_Uniform/", gevnum);
    
  cout << "Analyzed data will be of the type:\n " << strang << " .\n";
  
  // define flavor of this macro:  
  const char flavor[]="ePIC";

  //define and create, if not existing, pdf output directory:  
  TString pdfdir = TString("graphs/") + filetype_dir;
  if (!pdfdir.IsNull()) {
    // Directory does not exist: try to/Users/dhruv/Desktop/ePIC/FromCaroline make it
    gSystem->mkdir(pdfdir.Data(), kTRUE);
    cout << "Created output pdf directory:\n " << pdfdir << " .\n";
  }
    
  // define and open input file:
  TString infile_ram= TString("out_files/") + filetype_dir + "out." + strang + TString("-") + flavor + TString(".root");
  const char *infile=infile_ram.Data();
  TFile *ifile = TFile::Open(infile,"READ");
  
  cout << "Reading infile:\n " << infile << " .\n";

    // FILE 1 - muonEta histogram //
    TString name1 = TString("muonEta");
     
    TString filename1 = pdfdir + TString("/") + TString(name1) + TString(".pdf");
     
     TH1D *muonEta = (TH1D*)ifile->Get("muonEta");
     
     gStyle->SetOptStat(0);

     TCanvas *canvas1 = new TCanvas(name1, strang, 800, 600);

    muonEta->Draw("HIST");
    muonEta->Draw("E SAME");
    muonEta->Draw();
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.016);
    latex.DrawLatex(0.03, 0.91,
                    Form("env TILE_SIZE=10  N_LAYERS=10  SCINTILLATOR_THICKNESS=2.4  ABSORBER_THICKNESS=4  MOMENTUM=%d PARTICLE=mu-  NUMBER_OF_EVENTS=1000  ./run_submit.sh ",gevnum));
    
    canvas1->Print(filename1, "pdf");
     
     // End of muonEta Histogram

    // FILE 2 - nHCalRecHitsPosZ_all histogram //
    TString name2 = TString("nHCalRecHitsPosZ_all");
     
    TString filename2 = pdfdir + TString("/") + TString(name2) + TString(".pdf");
     
     TH1D *nHCalRecHitsPosZ_all = (TH1D*)ifile->Get("nHCalRecHitsPosZ_all");
     
     gStyle->SetOptStat(0);

     TCanvas *canvas2 = new TCanvas(name2, strang, 800, 600);

    nHCalRecHitsPosZ_all->Draw("HIST");
    nHCalRecHitsPosZ_all->Draw("E SAME");
    nHCalRecHitsPosZ_all->Draw();
    
    latex.SetNDC();
    latex.SetTextSize(0.016);
    latex.DrawLatex(0.03, 0.91,
                    Form("env TILE_SIZE=10  N_LAYERS=10  SCINTILLATOR_THICKNESS=2.4  ABSORBER_THICKNESS=4  MOMENTUM=%d PARTICLE=mu-  NUMBER_OF_EVENTS=1000  ./run_submit.sh ",gevnum));
    
    canvas2->Print(filename2, "pdf");
     // End of nHCalRecHitsPosZ_all Histogram

    // FILE 3 - nHCalRecHitsPosXY_all histogram //
    TString name3 = TString("nHCalRecHitsPosXY_all");
     
    TString filename3 = pdfdir + TString("/") + TString(name3) + TString(".pdf");
     
     TH2D *nHCalRecHitsPosXY_all = (TH2D*)ifile->Get("nHCalRecHitsPosXY_all");
     
     gStyle->SetOptStat(0);

     TCanvas *canvas3 = new TCanvas(name3, strang, 800, 600);

    nHCalRecHitsPosXY_all->Draw("HIST");
    nHCalRecHitsPosXY_all->Draw("E SAME");
    nHCalRecHitsPosXY_all->Draw();
    
    latex.SetNDC();
    latex.SetTextSize(0.016);
    latex.DrawLatex(0.03, 0.91,
                    Form("env TILE_SIZE=10  N_LAYERS=10  SCINTILLATOR_THICKNESS=2.4  ABSORBER_THICKNESS=4  MOMENTUM=%d PARTICLE=mu-  NUMBER_OF_EVENTS=1000  ./run_submit.sh ",gevnum));
    
    canvas3->Print(filename3, "pdf");
     
     // End of nHCalRecHitsPosXY_all Histogram
 
    // FILE 4 - nHCalRecHitsPosXYZ_all histogram //
    TString name4 = TString("nHCalRecHitsPosXYZ_all");
     
    TString filename4 = pdfdir + TString("/") + TString(name4) + TString(".pdf");
     
     TH1D *nHCalRecHitsPosXYZ_all = (TH1D*)ifile->Get("nHCalRecHitsPosXYZ_all");
     
     gStyle->SetOptStat(0);
     gStyle->SetPalette(kBird);
     TCanvas *canvas4 = new TCanvas(name4, strang, 800, 600);
    
    canvas4->SetLeftMargin(0.15);
    canvas4->SetRightMargin(0.20);   // room for colorbar
    canvas4->SetBottomMargin(0.15);
    canvas4->SetTopMargin(0.10);
    
    nHCalRecHitsPosXYZ_all->GetXaxis()->SetTitleOffset(1.8);
    nHCalRecHitsPosXYZ_all->GetYaxis()->SetTitleOffset(2.1);
    nHCalRecHitsPosXYZ_all->GetZaxis()->SetTitleOffset(1.8);
    //nHCalRecHitsPosXYZ_all->Draw("HIST");
    //nHCalRecHitsPosXYZ_all->Draw("E SAME");
    nHCalRecHitsPosXYZ_all->Draw("LEGO2Z");
    // Adjust these angles:
    // Caption at bottom
    latex.SetNDC();
    latex.SetTextSize(0.018);
    latex.DrawLatex(0.01, 0.03,
                    Form("env TILE_SIZE=10  N_LAYERS=10  SCINTILLATOR_THICKNESS=2.4  ABSORBER_THICKNESS=4  MOMENTUM=%d ",gevnum));
    latex.DrawLatex(0.01,0.01,"PARTICLE=mu-  NUMBER_OF_EVENTS=1000  ./run_submit.sh");
    canvas4->Update();  // must update first so the palette exists

    TPaletteAxis *palette =
        (TPaletteAxis*)nHCalRecHitsPosXYZ_all->GetListOfFunctions()->FindObject("palette");

    if (palette) {
        // Shift the palette to the right
        palette->SetX1NDC(0.88);  // left edge of colorbar
        palette->SetX2NDC(0.93);  // right edge of colorbar
        palette->SetY1NDC(0.15);  // bottom edge (optional)
        palette->SetY2NDC(0.90);  // top edge (optional)
    }

    canvas4->Modified();
    canvas4->Update();
    canvas4->Print(filename4, "pdf");
     
     // End of nHCalRecHitsPosXYZ_all Histogram

    // FILE 5 - nHCalRecHitsE_all histogram //
    TString name5 = TString("nHCalRecHitsE_all");
     
    TString filename5 = pdfdir + TString("/") + TString(name5) + TString(".pdf");
     
     TH1D *nHCalRecHitsE_all = (TH1D*)ifile->Get("nHCalRecHitsE_all");
     
     gStyle->SetOptStat(0);

     TCanvas *canvas5 = new TCanvas(name5, strang, 800, 600);

    nHCalRecHitsE_all->Draw("HIST");
    nHCalRecHitsE_all->Draw("E SAME");
    nHCalRecHitsE_all->GetXaxis()->SetTitle("nHCal Rec Hits x [mm]");
    nHCalRecHitsE_all->GetYaxis()->SetTitle("nHCal Rec Hits y [mm]");
    nHCalRecHitsE_all->GetZaxis()->SetTitle("nHCal Rec Hits z [mm]");
    nHCalRecHitsE_all->Draw();
    
    latex.SetNDC();
    latex.SetTextSize(0.016);
    latex.DrawLatex(0.03, 0.91,
                    Form("env TILE_SIZE=10  N_LAYERS=10  SCINTILLATOR_THICKNESS=2.4  ABSORBER_THICKNESS=4  MOMENTUM=%d PARTICLE=mu-  NUMBER_OF_EVENTS=1000  ./run_submit.sh ",gevnum));
    
    canvas5->Print(filename5, "pdf");
     
     // End of nHCalRecHitsE_all Histogram
    
    // FILE 6 - nHCalRecHitsE_all histogram //
    TString name6 = TString("nHCal_ByLayer_E");
     
    TString filename6 = pdfdir + TString("/") + TString(name6) + TString(".pdf");
     
    TH1D* nHCalLayers_Arr[nlayers_nhcal];
    for (int i = 0; i < nlayers_nhcal; i++) {
        nHCalLayers_Arr[i] = (TH1D*)ifile->Get(Form("nhCal_LayerE_%d",i));
    }
     
    gStyle->SetOptStat(0);

    TCanvas *canvas6 = new TCanvas(name6, strang, 1000, 600);
    
    canvas6->Divide(4,3);
    //gStyle->SetOptStat(1111);
    
    for (int i = 0; i < nlayers_nhcal; i++) {
        canvas6->cd(i+1);
        nHCalLayers_Arr[i]->Draw("HIST SAME");
    }
    
    canvas6->cd(nlayers_nhcal+1);
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.DrawLatex(0, 0.7, "env TILE_SIZE=10  N_LAYERS=10");
    latex.DrawLatex(0, 0.6, "SCINTILLATOR_THICKNESS=2.4");
    latex.DrawLatex(0, 0.5, "ABSORBER_THICKNESS=4");
    latex.DrawLatex(0, 0.4, Form("MOMENTUM=%d PARTICLE=mu-",gevnum));
    latex.DrawLatex(0, 0.3, "NUMBER_OF_EVENTS=1000  ./run_submit.sh");

    canvas6->Print(filename6, "pdf");
     
     // End of nHCalRecHitsE_all Histogram
    
    // FILE 7 - nHCalRecHitsE_all histogram //
    TString name7 = TString("nHCalClustersE_all");
     
    TString filename7 = pdfdir + TString("/") + TString(name7) + TString(".pdf");

    TH1D *nHCalClustersE_all = (TH1D*)ifile->Get("nHCalClustersE_all");
    
    gStyle->SetOptStat(0);

    TCanvas *canvas7 = new TCanvas(name7, strang, 800, 600);
    nHCalClustersE_all->Draw("HIST");
    
    latex.SetNDC();
    latex.SetTextSize(0.016);
    latex.DrawLatex(0.03, 0.91,
                    Form("env TILE_SIZE=10  N_LAYERS=10  SCINTILLATOR_THICKNESS=2.4  ABSORBER_THICKNESS=4  MOMENTUM=%d PARTICLE=mu-  NUMBER_OF_EVENTS=1000  ./run_submit.sh ",gevnum));

    canvas7->Print(filename7, "pdf");
     
     // End of nHCalRecHitsE_all Histogram

    // nHCalRecHitsPosZ_all_nonAssoc histogram //
    TString name8 = TString("nHCalRecHitsPosZ_all_nonAssoc");
     
    TString filename8 = pdfdir + TString("/") + TString(name8) + TString(".pdf");
     
     TH1D *nHCalRecHitsPosZ_all_nonAssoc = (TH1D*)ifile->Get("nHCalRecHitsPosZ_all_nonAssoc");
     
     gStyle->SetOptStat(0);

     TCanvas *canvas8 = new TCanvas(name8, strang, 800, 600);

    nHCalRecHitsPosZ_all_nonAssoc->Draw("HIST");
    nHCalRecHitsPosZ_all_nonAssoc->Draw("E SAME");
    nHCalRecHitsPosZ_all_nonAssoc->Draw();
    
    latex.SetNDC();
    latex.SetTextSize(0.016);
    latex.DrawLatex(0.03, 0.91,
                    Form("env TILE_SIZE=10  N_LAYERS=10  SCINTILLATOR_THICKNESS=2.4  ABSORBER_THICKNESS=4  MOMENTUM=%d PARTICLE=mu-  NUMBER_OF_EVENTS=1000  ./run_submit.sh ",gevnum));
    
    canvas8->Print(filename8, "pdf");
    
    // End of nHCalRecHitsPosZ_all_nonAssoc Histogram
    
    // FILE 9 - muonPhi histogram //
    
    TString name9 = TString("muonPhi");

    TString filename9 = pdfdir + TString("/") + TString(name9) + TString(".pdf");

     TH1D *muonPhi = (TH1D*)ifile->Get("muonPhi");

     gStyle->SetOptStat(0);

     TCanvas *canvas9 = new TCanvas(name9, strang, 800, 600);
    
    muonPhi->Draw("HIST");
    muonPhi->Draw("E SAME");
    muonPhi->Draw();

    latex.SetNDC();
    latex.SetTextSize(0.016);
    latex.DrawLatex(0.03, 0.91,
                    Form("env TILE_SIZE=10  N_LAYERS=10  SCINTILLATOR_THICKNESS=2.4  ABSORBER_THICKNESS=4  MOMENTUM=%d PARTICLE=mu-  NUMBER_OF_EVENTS=1000  ./run_submit.sh ",gevnum));
    
    canvas9->Print(filename9, "pdf");
    
     // End of muonPhi Histogram
    
    // FILE 10 - eta_v_phi_true graph //
    
    TString name10 = TString("eta_v_phi_true");
     
    TString filename10 = pdfdir + TString("/") + TString(name10) + TString(".pdf");
     
    TGraph *eta_v_phi_true = (TGraph*)ifile->Get("eta_v_phi_true");
    cout << (eta_v_phi_true == nullptr) << "\n";

    eta_v_phi_true->SetMinimum(1); //Set min val of q2 to 1
    
    //Make the points appear as larger red circles
    eta_v_phi_true->SetMarkerColor(kRed);
    eta_v_phi_true->SetMarkerStyle(20);
    eta_v_phi_true->SetMarkerSize(0.5);
    
    //Add Title and Axis Labels
    eta_v_phi_true->SetTitle("trueEta vs truePhi values");
    eta_v_phi_true->GetXaxis()->SetTitle("trueEta");
    eta_v_phi_true->GetXaxis()->CenterLabels();
    eta_v_phi_true->GetYaxis()->SetTitle("truePhi");
    eta_v_phi_true->GetYaxis()->CenterLabels();

    gStyle->SetOptStat(0);

    TCanvas *canvas10 = new TCanvas(name10, strang, 800, 600);
    canvas10->SetTopMargin(0.07);

    eta_v_phi_true->Draw("AP"); //Draw "A" - Axes, and "P" - Points(NOT lines)
    canvas10->Draw();
    
    latex.SetNDC();
    latex.SetTextSize(0.016);
    latex.DrawLatex(0.03, 0.91,
                    Form("env TILE_SIZE=10  N_LAYERS=10  SCINTILLATOR_THICKNESS=2.4  ABSORBER_THICKNESS=4  MOMENTUM=%d PARTICLE=mu-  NUMBER_OF_EVENTS=1000  ./run_submit.sh ",gevnum));
    
    canvas10->Print(filename10, "pdf");
    
    // End of eta_v_phi_true graph
    
    // FILE 11 - eta_v_phi_true_hist HISTOGRAM //
     
    // Define the name of the plot:
    TString name11 = TString("eta_v_phi_true_hist");
     
    // Define the name of the pdf file:
    TString filename11 = pdfdir + TString("/") + TString(name11) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *eta_v_phi_true_hist = (TH2F*)ifile->Get("eta_v_phi_true_hist");
    eta_v_phi_true_hist->SetMarkerColor(kWhite);
    eta_v_phi_true_hist->SetMarkerSize(2);
    eta_v_phi_true_hist->GetXaxis()->SetTitle("Eta");
    eta_v_phi_true_hist->GetYaxis()->SetTitle("Phi");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas11 = new TCanvas(name11, strang, 600, 600);
    canvas11->SetTopMargin(0.07);

     // Set color palette and axis labels
    canvas11->SetRightMargin(0.1);
    canvas11->SetLeftMargin(0.1);

    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    eta_v_phi_true_hist->Draw("COLZ");

    // Update the canvas
    canvas11->Update();
    
    canvas11->Draw();
    canvas11->Print(filename11, "pdf");

    // end of xB_q2 histogram
    
    
    
  //
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
  
} // end of macro  
//////////////////

  //TPaveText *t = new TPaveText(.05,.3,.95,.6, "NDC");                                                                                   
  //t->AddText("This line is blue"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);                                          
  //t->AddText("This line is red");  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kRed);                                           
  //t->Draw();     
