// Run this macro, from the Linux / Terminal command line, as
// root -l -q ePIC_Plotting.C
// The only condition is that previously ran ePIC_Analysis.C  with the identical "strang" as hardcoded in this macro (below)
// Strang must match the string in a file starting with out*, in this directory
// CR 2024-08-14/20

void ePIC_2Plotting()
{
  gSystem->Exec("date");
  
  // Define name of MC file:
  const char strang1[]="Sartre_Au_phi_10runs";
  const char filetype1_dir[]="Sartre/Au_phi/";
    
  const char strang2[]="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.0998.eicrecon.tree.edm4eic";
  const char filetype2_dir[]="Pythia/BeamEffects/";
  
  const char filetype_dir[]="2Graphs/";
    
  cout << "Analyzed data will be of the type:\n " << strang1 << " . and \n " << strang2 << " .\n";
  
  // define flavor of this macro:  
  const char flavor[]="ePIC";

  //define and create, if not existing, pdf output directory:  
  TString pdfdir = TString("graphs/") + filetype_dir + strang1;
    cout << pdfdir;
  if (!pdfdir.IsNull()) {
    // Directory does not exist: try to/Users/dhruv/Desktop/ePIC/FromCaroline make it
    gSystem->mkdir(pdfdir.Data(), kTRUE);
    cout << "Created output pdf directory:\n " << pdfdir << " .\n";
  }
  
  // define and open input file:
    TString name1 = TString("xBjorken");
    
    TString filename1 = pdfdir + TString("/") + TString(name1) + TString(".pdf");

    TString sartre_name = TString("out_files/") + filetype1_dir + "out." + strang1 + TString("-") + flavor + TString(".root");
    TFile *sartre_file = TFile::Open(sartre_name);
    TH1D *sartre_xb = (TH1D*)sartre_file->Get("xBjorken");
    
    TString pythia_name = TString("out_files/") + filetype2_dir + "out." + strang2 + TString("-") + flavor + TString(".root");
    TFile *pythia_file = TFile::Open(pythia_name);
    TH1D *pythia_xb = (TH1D*)pythia_file->Get("newXb");
         
    gStyle->SetOptStat(0);

    TCanvas *canvas1 = new TCanvas(name1, strang1, 800, 600);
    canvas1->SetLogx();
    canvas1->SetTopMargin(0.07);
    
    pythia_xb->SetLineColor(kBlack);
    pythia_xb->GetYaxis()->SetTitle("Counts (Pythia file)");
    pythia_xb->GetYaxis()->CenterTitle(true);
    pythia_xb->SetTitle("Sartre vs Pythia x_{Bj}");
    pythia_xb->Draw("HIST");
    
    canvas1->Update();
    
    gPad->Update();
    double rightmax = 1.05 * sartre_xb->GetMaximum();
    double scale = gPad->GetUymax() / rightmax;
    
    sartre_xb->SetLineColor(kBlue);
    sartre_xb->SetLineWidth(1);
    sartre_xb->Scale(scale);
    sartre_xb->Draw("HIST SAME");
    
    auto axis = new TGaxis(2,gPad->GetUymin(), 2, gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetLineColor(kBlue);
    axis->SetTextColor(kBlue);
    axis->SetLabelColor(kBlue);
    axis->SetLabelFont(42);
    axis->SetLabelSize(0.03);
    axis->SetTitle("Counts (Sartre File)");
    axis->SetTitleFont(42);
    axis->CenterTitle(true);
    axis->Draw();
    
    canvas1->Draw();
    canvas1->Print(filename1, "pdf");
  
  
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
  
} // end of macro  
//////////////////
