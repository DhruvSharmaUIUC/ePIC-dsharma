// Run this macro, from the Linux / Terminal command line, as
// root -l -q ePIC_Plotting.C
// The only condition is that previously ran ePIC_Analysis.C  with the identical "strang" as hardcoded in this macro (below)
// Strang must match the string in a file starting with out*, in this directory
// CR 2024-08-14/20

void ePIC_Plotting()
{
  gSystem->Exec("date");

  // define nHCal acceptance (2024-07-16) - make sure this is consistent with *_Analysis.C:
  double eta_min = -4.05;
  double eta_max = -1.2;
  //
  
  // Define name of MC file:
  const char strang[]="Sartre_Au_phi_10runs";
  const char filetype_dir[]="Sartre/Au_phi/";
    
  cout << "Analyzed data will be of the type:\n " << strang << " .\n";
  
  // define flavor of this macro:  
  const char flavor[]="ePIC";

  //define and create, if not existing, pdf output directory:  
  TString pdfdir = TString("graphs/") + filetype_dir + strang;
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
  
  // FILE 1 - generated eta //
  
  // Define the name of the plot:
  TString name1 = TString("trueEta_species");
  // Define the name of the pdf file:
  TString filename1 = pdfdir + TString("/") + TString(name1) + TString(".pdf");

  TH1F *electronEta = (TH1F*)ifile->Get("electronEta");
  TH1F *muonEta = (TH1F*)ifile->Get("muonEta");
  TH1F *protonEta = (TH1F*)ifile->Get("protonEta");
  TH1F *pionEta = (TH1F*)ifile->Get("pionEta");
  TH1F *kaonEta = (TH1F*)ifile->Get("kaonEta");
  TH1F *rho0Eta = (TH1F*)ifile->Get("rho0Eta");
  TH1F *phiEta = (TH1F*)ifile->Get("phiEta");

  gStyle->SetOptStat(0); //no stats box
  
  TCanvas *canvas = new TCanvas(name1, strang, 800, 600);
    canvas->SetTopMargin(0.07);
  electronEta->SetTitle("Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay K+K-, proton, electron distribution");
  electronEta->SetLineColor(kBlue);
  electronEta->Draw();
  muonEta->SetLineColor(kMagenta);
  muonEta->Draw("same");
  protonEta->SetLineColor(kRed);
  protonEta->Draw("same");
  pionEta->SetLineColor(kOrange+1);
  pionEta->Draw("same");
  kaonEta->SetLineColor(kCyan);
  kaonEta->Draw("same");
  rho0Eta->SetLineColor(kBlack);
  rho0Eta->Draw("same");
    
  canvas->Draw();

  auto leg = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header  
  leg->SetHeader("Generated particles", "C"); // option "C" allows to center the header
  leg->SetFillStyle(0);
  leg->AddEntry(electronEta,"electrons (#pm)","l");
  leg->AddEntry(muonEta,"muons (#pm)","l");
  leg->AddEntry(pionEta,"pions (#pm)","l");
  leg->AddEntry(kaonEta,"kaons (#pm)","l");
  leg->AddEntry(protonEta,"protons (#pm)","l");
  leg->AddEntry(rho0Eta,"#rho^{0} (770)","l");
  leg->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_1 = electronEta->GetMaximumBin();
  Double_t y_max1 = electronEta->GetBinContent(binmax_1);
  TLine *eta_min_line1= new TLine(eta_min,0.,eta_min,y_max1);  // (x1,y1,x2,y2)
  eta_min_line1->SetLineColor(kBlack);
  eta_min_line1->SetLineWidth(2);
  eta_min_line1->SetLineStyle(kDashed);
  eta_min_line1->Draw("same");
  TLine *eta_max_line1= new TLine(eta_max,0.,eta_max,y_max1);  // (x1,y1,x2,y2)
  eta_max_line1->SetLineColor(kBlack);
  eta_max_line1->SetLineWidth(2);
  eta_max_line1->SetLineStyle(kDashed);
  eta_max_line1->Draw("same");

  //
  canvas->Print(filename1, "pdf");          
  // end of generated eta

  // FILE 2 - reconstructed eta //
  
  // Define the name of the plot:
  TString name2 = TString("recEta_species");
  // Define the name of the pdf file:
  TString filename2 = pdfdir + TString("/") + TString(name2) + TString(".pdf");

  TH1F *electronRecEta = (TH1F*)ifile->Get("electronRecEta");
  TH1F *muonRecEta = (TH1F*)ifile->Get("muonRecEta");
  TH1F *protonRecEta = (TH1F*)ifile->Get("protonRecEta");
  TH1F *pionRecEta = (TH1F*)ifile->Get("pionRecEta");
  TH1F *kaonRecEta = (TH1F*)ifile->Get("kaonRecEta");
  

  gStyle->SetOptStat(0); //no stats box    

  TCanvas *canvas2 = new TCanvas(name2, strang, 800, 600);
  electronRecEta->SetTitle(strang);
  electronRecEta->SetLineColor(kBlue);
  electronRecEta->Draw();
  muonRecEta->SetLineColor(kMagenta);
  muonRecEta->Draw("same");
  protonRecEta->SetLineColor(kRed);
  protonRecEta->Draw("same");
  pionRecEta->SetLineColor(kOrange+1);
  pionRecEta->Draw("same");
  kaonRecEta->SetLineColor(kCyan);
  kaonRecEta->Draw("same");

  canvas2->Draw();

  auto leg2 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header                                                                         
  leg2->SetHeader("Reconstructed particles", "C"); // option "C" allows to center the header
  leg2->SetFillStyle(0);
  leg2->AddEntry(electronRecEta,"electrons (#pm)","l");
  leg2->AddEntry(muonRecEta,"muons (#pm)","l");
  leg2->AddEntry(pionRecEta,"pions (#pm)","l");
  leg2->AddEntry(kaonRecEta,"kaons (#pm)","l");
  leg2->AddEntry(protonRecEta,"protons (#pm)","l");
  leg2->Draw();

  // add vertical lines for nHCal acceptance                                                                                              
  Int_t binmax_2 = electronRecEta->GetMaximumBin();
  Double_t y_max2 = electronRecEta->GetBinContent(binmax_2);
  TLine *eta_min_line2= new TLine(eta_min,0.,eta_min,y_max2);  // (x1,y1,x2,y2)
  eta_min_line2->SetLineColor(kBlack);
  eta_min_line2->SetLineWidth(2);
  eta_min_line2->SetLineStyle(kDashed);
  eta_min_line2->Draw("same");
  TLine *eta_max_line2= new TLine(eta_max,0.,eta_max,y_max2);  // (x1,y1,x2,y2)
  eta_max_line2->SetLineColor(kBlack);
  eta_max_line2->SetLineWidth(2);
  eta_max_line2->SetLineStyle(kDashed);
  eta_max_line2->Draw("same");
  
  //
  canvas2->Print(filename2, "pdf");
  // end of reconstructed eta  

  // FILE 3 - gen & rec eta: //

  // Define the name of the plot:
  TString name3 = TString("gen-recEta_species");
  // Define the name of the pdf file:
  TString filename3 = pdfdir + TString("/") + TString(name3) + TString(".pdf");

  TCanvas *canvas3 = new TCanvas(name3, strang, 1200, 600);

  //devide the canvas into several pads (x,y):
  canvas3->Divide(3,2);

  gStyle->SetOptStat(1111); //now we want some stat boxes  
  
  // go to the first one and do your thing
  canvas3->cd(1);
  electronEta->SetTitle(strang);
  electronEta->SetLineColor(kBlue);
  electronEta->Draw();
  electronRecEta->SetLineStyle(2);
  electronRecEta->Draw("same");

  auto leg31 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
  leg31->SetBorderSize(0);
  leg31->SetTextSize(0.05);
  leg31->SetFillStyle(0);
  leg31->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg31->AddEntry(electronEta,"electrons (#pm) gen","l");
  leg31->AddEntry(electronRecEta,"electrons (#pm) rec","l");
  leg31->Draw();

  // add vertical lines for nHCal acceptance (use existing lines for generated eta)
  eta_min_line1->Draw("same");
  eta_max_line1->Draw("same");
  
  // pad 2
  canvas3->cd(2);
  muonEta->SetTitle(strang);
  muonEta->SetLineColor(kMagenta);
  muonEta->Draw();
  muonRecEta->SetLineStyle(2);
  muonRecEta->Draw("same");

  auto leg32 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header                                                                       
  leg32->SetBorderSize(0);
  leg32->SetTextSize(0.05);
  leg32->SetFillStyle(0);
  leg32->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg32->AddEntry(muonEta,"muons (#pm) gen","l");
  leg32->AddEntry(muonRecEta,"muons (#pm) rec","l");
  leg32->Draw();

// add vertical lines for nHCal acceptance
  Int_t binmax_32 = muonEta->GetMaximumBin();
  Double_t y_max32 = muonEta->GetBinContent(binmax_32);
  TLine *eta_min_line32= new TLine(eta_min,0.,eta_min,y_max32);
  eta_min_line32->SetLineColor(kBlack);
  eta_min_line32->SetLineWidth(2);
  eta_min_line32->SetLineStyle(kDashed);
  eta_min_line32->Draw("same");
  TLine *eta_max_line32= new TLine(eta_max,0.,eta_max,y_max32);
  eta_max_line32->SetLineColor(kBlack);
  eta_max_line32->SetLineWidth(2);
  eta_max_line32->SetLineStyle(kDashed);
  eta_max_line32->Draw("same");

  // pad 3
  canvas3->cd(3);
  pionEta->SetTitle(strang);
  pionEta->SetLineColor(kOrange+1);
  pionEta->Draw();
  pionRecEta->SetLineStyle(2);
  pionRecEta->Draw("same");

  auto leg33 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
                                                                                                                                          
  leg33->SetBorderSize(0);
  leg33->SetTextSize(0.05);
  leg33->SetFillStyle(0);
  leg33->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header                                          
  leg33->AddEntry(pionEta,"pions (#pm) gen","l");
  leg33->AddEntry(pionRecEta,"pions (#pm) rec","l");
  leg33->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_33 = pionEta->GetMaximumBin();
  Double_t y_max33 = pionEta->GetBinContent(binmax_33);
  TLine *eta_min_line33= new TLine(eta_min,0.,eta_min,y_max33);  
  eta_min_line33->SetLineColor(kBlack);
  eta_min_line33->SetLineWidth(2);
  eta_min_line33->SetLineStyle(kDashed);
  eta_min_line33->Draw("same");
  TLine *eta_max_line33= new TLine(eta_max,0.,eta_max,y_max33);  
  eta_max_line33->SetLineColor(kBlack);
  eta_max_line33->SetLineWidth(2);
  eta_max_line33->SetLineStyle(kDashed);
  eta_max_line33->Draw("same");
  
  // pad 4
  canvas3->cd(4);
  protonEta->SetTitle(strang);
  protonEta->SetLineColor(kRed);
  protonEta->Draw();
  protonRecEta->SetLineStyle(2);
  protonRecEta->Draw("same");

  auto leg34 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
					      
  leg34->SetBorderSize(0);
  leg34->SetTextSize(0.05);
  leg34->SetFillStyle(0);
  leg34->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg34->AddEntry(protonEta,"protons (#pm) gen","l");
  leg34->AddEntry(protonRecEta,"protons (#pm) rec","l");
  leg34->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_34 = protonEta->GetMaximumBin();
  Double_t y_max34 = protonEta->GetBinContent(binmax_34);
  TLine *eta_min_line34= new TLine(eta_min,0.,eta_min,y_max34);
  eta_min_line34->SetLineColor(kBlack);
  eta_min_line34->SetLineWidth(2);
  eta_min_line34->SetLineStyle(kDashed);
  eta_min_line34->Draw("same");
  TLine *eta_max_line34= new TLine(eta_max,0.,eta_max,y_max34);
  eta_max_line34->SetLineColor(kBlack);
  eta_max_line34->SetLineWidth(2);
  eta_max_line34->SetLineStyle(kDashed);
  eta_max_line34->Draw("same");

  // pad 5
  canvas3->cd(5);
  kaonEta->SetTitle(strang);
  kaonEta->SetLineColor(kCyan);
  kaonEta->Draw();
  kaonRecEta->SetLineStyle(2);
  kaonRecEta->Draw("same");

  auto leg35 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
  						  
  leg35->SetBorderSize(0);
  leg35->SetTextSize(0.05);
  leg35->SetFillStyle(0);
  leg35->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg35->AddEntry(kaonEta,"kaons (#pm) gen","l");
  leg35->AddEntry(kaonRecEta,"kaons (#pm) rec","l");
  leg35->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_35 = kaonEta->GetMaximumBin();
  Double_t y_max35 = kaonEta->GetBinContent(binmax_35);
  TLine *eta_min_line35= new TLine(eta_min,0.,eta_min,y_max35);
  eta_min_line35->SetLineColor(kBlack);
  eta_min_line35->SetLineWidth(2);
  eta_min_line35->SetLineStyle(kDashed);
  eta_min_line35->Draw("same");
  TLine *eta_max_line35= new TLine(eta_max,0.,eta_max,y_max35);
  eta_max_line35->SetLineColor(kBlack);
  eta_max_line35->SetLineWidth(2);
  eta_max_line35->SetLineStyle(kDashed);
  eta_max_line35->Draw("same");
  
  // draw canvas and print to pdf
  canvas3->Draw();
  canvas3->Print(filename3, "pdf");
  // end of gen & rec eta       

  // FILE 4 - eta decay rho0 to pi+ pi- //                                                                                                           
  // Define the name of the plot:
  TString name4 = TString("Eta_decay_rho0");
  // Define the name of the pdf file:
  TString filename4 = pdfdir + TString("/") + TString(name4) + TString(".pdf");

  TH1F *pipmfromrho0Eta = (TH1F*)ifile->Get("pipmfromrho0Eta");
  TH1F *pipmfromrho0RecEta = (TH1F*)ifile->Get("pipmfromrho0RecEta");
  // exists already: 
  //TH1F *rho0Eta = (TH1F*)ifile->Get("rho0Eta");

  gStyle->SetOptStat(0); //no stats box                 

  TCanvas *canvas4 = new TCanvas(name4, strang, 800, 600);
  canvas4->SetTopMargin(0.07);
  rho0Eta->SetLineColor(kBlack);
  rho0Eta->SetLineStyle(1);
  rho0Eta->Draw();
  pipmfromrho0RecEta->SetLineStyle(2);
  pipmfromrho0RecEta->SetTitle(strang);
  pipmfromrho0RecEta->SetLineColor(kRed);
  pipmfromrho0RecEta->Draw("same");
  
  canvas4->Draw();

  auto leg4 = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetTextSize(0.05);
  leg4->SetHeader("generated #rho^{0}(770) and decay pions", "C"); 
  leg4->AddEntry(pipmfromrho0RecEta,"reco pions (#pm)","l");
  leg4->AddEntry(rho0Eta,"gen #rho^{0}","l");
  leg4->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_4 = pipmfromrho0Eta->GetMaximumBin();
  Double_t y_max4 = 0.6*pipmfromrho0Eta->GetBinContent(binmax_4);
  TLine *eta_min_line4= new TLine(eta_min,0.,eta_min,y_max4); 
  eta_min_line4->SetLineColor(kBlack);
  eta_min_line4->SetLineWidth(2);
  eta_min_line4->SetLineStyle(kDashed);
  eta_min_line4->Draw("same");
  TLine *eta_max_line4= new TLine(eta_max,0.,eta_max,y_max4);
  eta_max_line4->SetLineColor(kBlack);
  eta_max_line4->SetLineWidth(2);
  eta_max_line4->SetLineStyle(kDashed);
  eta_max_line4->Draw("same");
  //                                                                                                                                        
  canvas4->Print(filename4, "pdf");
  // end of decay rho0 to pipi eta   

   // FILE 5 - eta decay phi to K+ K- //                                                                                                           
  // Define the name of the plot:
  TString name5 = TString("Eta_decay_phi");
  // Define the name of the pdf file:
  TString filename5 = pdfdir + TString("/") + TString(name5) + TString(".pdf");

  TH1F *kpmfromphiEta = (TH1F*)ifile->Get("kpmfromphiEta");
  TH1F *kpmfromphiRecEta = (TH1F*)ifile->Get("kpmfromphiRecEta");
    

  gStyle->SetOptStat(0); //no stats box                 

  TCanvas *canvas5 = new TCanvas(name5, strang, 800, 400);
    canvas5->SetTopMargin(0.07);

  // meager channels:
  //phiEta->SetLineColor(kBlack);
  //phiEta->SetLineStyle(1);
  //phiEta->Draw();
  //kpmfromphiRecEta->SetLineStyle(2);
  //kpmfromphiRecEta->SetTitle(strang);
  //kpmfromphiRecEta->SetLineColor(kRed);
  //kpmfromphiRecEta->Draw("same");
  // fat channels:
  kpmfromphiRecEta->SetLineStyle(2);
  kpmfromphiRecEta->SetTitle("Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', Gen vs Reco decay K+K-");
  kpmfromphiRecEta->SetLineColor(kRed);
  kpmfromphiRecEta->Draw();
    kaonEta->SetLineStyle(1);
    kaonEta->SetLineColor(kBlack);
    kaonEta->SetLineStyle(1);
  kaonEta->Draw("same");
  /*phiEta->SetLineStyle(1);
  phiEta->SetLineColor(kBlack);
  phiEta->SetLineStyle(1);
  phiEta->Draw("same");
  canvas5->Draw();*/

  auto leg5 = new TLegend(0.6,0.7,0.78,0.8); //x1,y1,x2,y2,header
  leg5->SetFillStyle(0);
  leg5->AddEntry(kpmfromphiRecEta,"reco kaons (#pm)","l");
  leg5->AddEntry(phiEta,"gen kaons","l");
  leg5->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_5 = kpmfromphiEta->GetMaximumBin();
  Double_t y_max5 = 0.9*kpmfromphiEta->GetBinContent(binmax_5);
  TLine *eta_min_line5= new TLine(eta_min,0.,eta_min,y_max5);
  eta_min_line5->SetLineColor(kBlack);
  eta_min_line5->SetLineWidth(2);
  eta_min_line5->SetLineStyle(kDashed);
  eta_min_line5->Draw("same");
  TLine *eta_max_line5= new TLine(eta_max,0.,eta_max,y_max5);
  eta_max_line5->SetLineColor(kBlack);
  eta_max_line5->SetLineWidth(2);
  eta_max_line5->SetLineStyle(kDashed);
  eta_max_line5->Draw("same");
  TLine *LFHCAL_eta_min_line= new TLine(1.18,0.,1.18,y_max5);  // (x1,y1,x2,y2)
    LFHCAL_eta_min_line->SetLineColor(kBlack);
    LFHCAL_eta_min_line->SetLineWidth(2);
    LFHCAL_eta_min_line->SetLineStyle(kDashed);
    LFHCAL_eta_min_line->Draw("same");
  TLine *LFHCAL_eta_max_line= new TLine(4.2,0.,4.2,y_max5);  // (x1,y1,x2,y2)
    LFHCAL_eta_max_line->SetLineColor(kBlack);
    LFHCAL_eta_max_line->SetLineWidth(2);
    LFHCAL_eta_max_line->SetLineStyle(kDashed);
    LFHCAL_eta_max_line->Draw("same");
  //
  canvas5->Print(filename5, "pdf");
  // end of decay phi to KK eta   
  
    
    
   // FILE 6 - kaonOccurrence HISTOGRAM //
    
   // Define the name of the plot:
   TString name6 = TString("kaonOccurrence");
    
   // Define the name of the pdf file:
   TString filename6 = pdfdir + TString("/") + TString(name6) + TString(".pdf");
    
   //Get the histogram data from ePIC_Analysis.C
    TH2D *kaonOccurrence = (TH2D*)ifile->Get("kaonOccurrence");
    kaonOccurrence->SetMinimum(-0.001);
    kaonOccurrence->SetMarkerColor(kWhite);
    kaonOccurrence->SetMarkerSize(2);
    //Set bin labels
    kaonOccurrence->GetXaxis()->SetBinLabel(1,"nHCal");
    kaonOccurrence->GetXaxis()->SetBinLabel(2,"Barrel");
    kaonOccurrence->GetXaxis()->SetBinLabel(3,"LFHCAL");
    kaonOccurrence->GetXaxis()->SetBinLabel(4,"Any HCal");
    
    kaonOccurrence->GetYaxis()->SetBinLabel(1,"nHCal");
    kaonOccurrence->GetYaxis()->SetBinLabel(2,"Barrel");
    kaonOccurrence->GetYaxis()->SetBinLabel(3,"LFHCAL");
    kaonOccurrence->GetYaxis()->SetBinLabel(4,"Any HCal");
    
    kaonOccurrence->GetZaxis()->SetTitle("geom acc. (%)");
    kaonOccurrence->GetZaxis()->SetTitleOffset(1.3);
    
    gStyle->SetOptStat(0);
    TCanvas *canvas6 = new TCanvas(name6, strang, 600, 600);
    canvas6->SetTopMargin(0.07);

    // Set color palette and axis labels
    canvas6->SetRightMargin(0.15);
    canvas6->SetLeftMargin(0.15);

    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    kaonOccurrence->Draw("COLZ TEXT");

   // Update the canvas
   canvas6->Update();
   
   canvas6->Draw();
   canvas6->Print(filename6, "pdf");

   // end of kaon occurrence histogram
    
   // FILE 7 - xBjorken histogram //
   TString name7 = TString("xBjorken");
    
   TString filename7 = pdfdir + TString("/") + TString(name7) + TString(".pdf");
    
    TH1D *xBjorken = (TH1D*)ifile->Get("xBjorken");
    
    gStyle->SetOptStat(0);

    TCanvas *canvas7 = new TCanvas(name7, strang, 800, 600);
    canvas7->SetLogx();
    canvas7->SetTopMargin(0.07);
    
    xBjorken->Draw("HIST");
    xBjorken->Draw("E SAME");
    canvas7->Draw();
    canvas7->Print(filename7, "pdf");
    
    // End of xBjorken Histogram
    
// FILE 8 - q^2 histogram //
    TString name8 = TString("q2");
     
    TString filename8 = pdfdir + TString("/") + TString(name8) + TString(".pdf");
     
    TH1D *q2 = (TH1D*)ifile->Get("q2");
     
    gStyle->SetOptStat(0);

    TCanvas *canvas8 = new TCanvas(name8, strang, 800, 600);
    canvas8->SetLogx();
    canvas8->SetTopMargin(0.07);
    q2->Draw("HIST");
    q2->Draw("E SAME");
    canvas8->Draw();
    canvas8->Print(filename8, "pdf");
     
// End of FILE 8 - q^2 Histogram //
    
// FILE 9 - x_B vs q^2 Scatterplot  //
    
    TString name9 = TString("xB_v_q2");
     
    TString filename9 = pdfdir + TString("/") + TString(name9) + TString(".pdf");
     
    TGraph *xB_v_q2 = (TGraph*)ifile->Get("xB_v_q2");
    
    xB_v_q2->SetMinimum(1); //Set min val of q2 to 1
    
    //Make the points appear as larger red circles
    xB_v_q2->SetMarkerColor(kRed);
    xB_v_q2->SetMarkerStyle(20);
    xB_v_q2->SetMarkerSize(0.5);
    
    //Add Title and Axis Labels
    xB_v_q2->SetTitle("x_{Bj} vs Q^{2} values");
    xB_v_q2->GetXaxis()->SetTitle("log_{10}x_{Bj}");
    xB_v_q2->GetXaxis()->CenterLabels();
    xB_v_q2->GetYaxis()->SetTitle("log_{10}Q^{2}");
    xB_v_q2->GetYaxis()->CenterLabels();

    gStyle->SetOptStat(0);

    TCanvas *canvas9 = new TCanvas(name9, strang, 800, 600);
    canvas9->SetLogx(); //Set log scales
    canvas9->SetLogy();
    canvas9->SetTopMargin(0.07);

    xB_v_q2->Draw("AP"); //Draw "A" - Axes, and "P" - Points(NOT lines)
    canvas9->Draw();
    canvas9->Print(filename9, "pdf");
     
// End of FILE 9 - xB_vs_q2 scatterplot //
    
// FILE 10 - xB_v_percent plot //
    
    TString name10 = TString("xB_v_percent_0");
     
    TString filename10 = pdfdir + TString("/") + TString(name10) + TString(".pdf");
     
    TH1D *xB_v_percent_0 = (TH1D*)ifile->Get("xB_v_percent_0");
    TH1D *xB_v_percent_1 = (TH1D*)ifile->Get("xB_v_percent_1");
    TH1D *xB_v_percent_2 = (TH1D*)ifile->Get("xB_v_percent_2");
    TH1D *xB_v_percent_all = (TH1D*)ifile->Get("xB_v_percent_all");
    
    gStyle->SetOptStat(0);

    TCanvas *canvas10 = new TCanvas(name10, strang, 600, 600);
    canvas10->SetLogx();
    canvas10->SetTopMargin(0.07);

    xB_v_percent_0->SetMarkerStyle(8);
    xB_v_percent_0->SetMarkerColor(kRed);
    xB_v_percent_0->SetLineColor(kRed);
    xB_v_percent_0->SetLineWidth(3);
    xB_v_percent_0->SetMaximum(1.4);
    xB_v_percent_0->SetMinimum(-0.001);
    xB_v_percent_0->GetXaxis()->SetTitleOffset(1.2);
    
    xB_v_percent_1->SetMarkerStyle(8);
    xB_v_percent_1->SetMarkerColor(kBlue);
    xB_v_percent_1->SetLineColor(kBlue);
    xB_v_percent_1->SetLineWidth(3);
    
    xB_v_percent_2->SetMarkerStyle(8);
    xB_v_percent_2->SetMarkerColor(kMagenta);
    xB_v_percent_2->SetLineColor(kMagenta);
    xB_v_percent_2->SetLineWidth(3);
    
    xB_v_percent_all->SetMarkerStyle(8);
    xB_v_percent_all->SetMarkerColor(kBlack);
    xB_v_percent_all->SetLineColor(kBlack);
    xB_v_percent_all->SetLineWidth(3);
    
    gStyle->SetEndErrorSize(0);
    
    xB_v_percent_0->Draw("E1X0 P");
    xB_v_percent_1->Draw("E1X0 P SAME");
    xB_v_percent_2->Draw("E1X0 P SAME");
    xB_v_percent_all->Draw("E1X0 P SAME");
    
    auto xB_percent_leg = new TLegend(0.48,0.75,0.78,0.88);
    xB_percent_leg->SetFillStyle(0);
    xB_percent_leg->AddEntry(xB_v_percent_0, "0-K in nHCAL", "P");
    xB_percent_leg->AddEntry(xB_v_percent_1, "1-K in nHCAL", "P");
    xB_percent_leg->AddEntry(xB_v_percent_2, "2-K in nHCAL", "P");
    xB_percent_leg->AddEntry(xB_v_percent_all, "All Calo", "P");
    xB_percent_leg->Draw();
    
    TLine *hundred_line = new TLine(0.,1.,1,1.);  // (x1,y1,x2,y2)
    hundred_line->SetLineColor(kBlack);
    hundred_line->SetLineWidth(2);
    hundred_line->SetLineStyle(kDashed);
    hundred_line->Draw("same");
    
    canvas10->Draw();
    canvas10->Print(filename10, "pdf");
     
// End of xB_v_percent plot
    
    // FILE 11 - xB_q2 HISTOGRAM //
     
    // Define the name of the plot:
    TString name11 = TString("xB_q2_hist");
     
    // Define the name of the pdf file:
    TString filename11 = pdfdir + TString("/") + TString(name11) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *xB_q2_hist = (TH2F*)ifile->Get("xB_q2_hist");
    xB_q2_hist->SetMarkerColor(kWhite);
    xB_q2_hist->SetMarkerSize(2);
    xB_q2_hist->GetXaxis()->SetTitle("x_{Bj}");
    xB_q2_hist->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas11 = new TCanvas(name11, strang, 600, 600);
    canvas11->SetLogx();
    canvas11->SetLogy();
    canvas11->SetTopMargin(0.07);

     // Set color palette and axis labels
    canvas11->SetRightMargin(0.1);
    canvas11->SetLeftMargin(0.1);

    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    xB_q2_hist->Draw("COLZ");

    // Update the canvas
    canvas11->Update();
    
    canvas11->Draw();
    canvas11->Print(filename11, "pdf");

    // end of xB_q2 histogram
    
    // start of all eta histogram
    
    TString name12 = TString("all_eta");
     
    TString filename12 = pdfdir + TString("/") + TString(name12) + TString(".pdf");
     
    TH1D *all_eta = (TH1D*)ifile->Get("all_eta");
     
    gStyle->SetOptStat(0);
    //gStyle->SetTitleOffset(1.2, "y");

    TCanvas *canvas12 = new TCanvas(name12, strang, 800, 600);
    canvas12->SetTopMargin(0.07);
    
    all_eta->Draw("HIST");
    all_eta->Draw("E SAME");
    
    TLine *nHCAL_min_line= new TLine(-4.05,0.,-4.05,750);  // (x1,y1,x2,y2)
    nHCAL_min_line->SetLineColor(kBlack);
    nHCAL_min_line->SetLineWidth(2);
    nHCAL_min_line->SetLineStyle(kDashed);
    nHCAL_min_line->Draw("same");
    TLine *nHCAL_max_line= new TLine(-1.2,0.,-1.2,750);  // (x1,y1,x2,y2)
    nHCAL_max_line->SetLineColor(kBlack);
    nHCAL_max_line->SetLineWidth(2);
    nHCAL_max_line->SetLineStyle(kDashed);
    nHCAL_max_line->Draw("same");
    TLine *LFHCAL_min_line= new TLine(1.18,0.,1.18,750);  // (x1,y1,x2,y2)
    LFHCAL_min_line->SetLineColor(kBlack);
    LFHCAL_min_line->SetLineWidth(2);
    LFHCAL_min_line->SetLineStyle(kDashed);
    LFHCAL_min_line->Draw("same");
    TLine *LFHCAL_max_line= new TLine(4.2,0.,4.2,750);  // (x1,y1,x2,y2)
    LFHCAL_max_line->SetLineColor(kBlack);
    LFHCAL_max_line->SetLineWidth(2);
    LFHCAL_max_line->SetLineStyle(kDashed);
    LFHCAL_max_line->Draw("same");
    
    canvas12->Draw();
    canvas12->Print(filename12, "pdf");
    
    // end of all eta histogram
    /*
    // start of New q2 Histogram
    TString name14 = TString("newQ2");
     
    TString filename14 = pdfdir + TString("/") + TString(name14) + TString(".pdf");
     
    TH1D *newQ2 = (TH1D*)ifile->Get("newQ2");
     
    gStyle->SetOptStat(0);

    TCanvas *canvas14 = new TCanvas(name14, strang, 800, 600);
    canvas14->SetLogx();
    canvas14->SetTopMargin(0.07);
    newQ2->Draw("HIST");
    //newQ2->Draw("E SAME");
    canvas14->Draw();
    canvas14->Print(filename14, "pdf");
    // end of new q2 Histogram
    
    // start of New Xb histogram
    TString name15 = TString("newXb");
     
    TString filename15 = pdfdir + TString("/") + TString(name15) + TString(".pdf");
     
    TH1D *newXb = (TH1D*)ifile->Get("newXb");
     
    gStyle->SetOptStat(0);

    TCanvas *canvas15 = new TCanvas(name15, strang, 800, 600);
    canvas15->SetLogx();
    canvas15->SetTopMargin(0.07);
    newXb->Draw("HIST");
    //newXb->Draw("E SAME");
    canvas15->Draw();
    canvas15->Print(filename15, "pdf");
    // end of new Xb histogram
    */
    
    // start of xb by num in nHCal histogram
    
    TString name16 = TString("Xb_num_in_nHCal");
    
    TString filename16 = pdfdir + TString("/") + TString(name16) + TString(".pdf");
    
    TH1D *xBjorken_0 = (TH1D*)ifile->Get("xBjorken0");
    TH1D *xBjorken_1 = (TH1D*)ifile->Get("xBjorken1");
    TH1D *xBjorken_2 = (TH1D*)ifile->Get("xBjorken2");
    
    gStyle->SetOptStat(0);
    
    TCanvas *canvas16 = new TCanvas(name16, strang, 800, 600);
    canvas16->SetLogx();
    canvas16->SetTopMargin(0.07);
    xBjorken->SetLineColor(kBlack);
    xBjorken->SetLineStyle(kDashed);
    xBjorken->Draw("HIST");
    
    xBjorken_0->SetLineColor(kRed);
    xBjorken_0->Draw("HIST SAME");
    
    xBjorken_1->SetLineColor(kBlue);
    xBjorken_1->Draw("HIST SAME");
    
    xBjorken_2->SetLineColor(kMagenta);
    xBjorken_2->Draw("HIST SAME");
    
    auto xB_num_in_nHCal_leg = new TLegend(0.65,0.75,0.85,0.9);
    xB_num_in_nHCal_leg->SetFillStyle(0);
    xB_num_in_nHCal_leg->AddEntry(xBjorken, "All events", "l");
    xB_num_in_nHCal_leg->AddEntry(xBjorken_0, "0-K in nHCAL", "l");
    xB_num_in_nHCal_leg->AddEntry(xBjorken_1, "1-K in nHCAL", "l");
    xB_num_in_nHCal_leg->AddEntry(xBjorken_2, "2-K in nHCAL", "l");
    xB_num_in_nHCal_leg->Draw();
    
    canvas16->Draw();
    canvas16->Print(filename16, "pdf");
    
    // End of xb by num in nHCal Histogram
    
    //start of q2 by num in nhcal histogram
    
    TString name17 = TString("Q2_num_in_nHCal");
    
    TString filename17 = pdfdir + TString("/") + TString(name17) + TString(".pdf");
    
    TH1D *q2_0 = (TH1D*)ifile->Get("q2_0");
    TH1D *q2_1 = (TH1D*)ifile->Get("q2_1");
    TH1D *q2_2 = (TH1D*)ifile->Get("q2_2");
    
    gStyle->SetOptStat(0);
    
    TCanvas *canvas17 = new TCanvas(name17, strang, 800, 600);
    canvas17->SetLogx();
    canvas17->SetTopMargin(0.07);
    q2->SetLineColor(kBlack);
    q2->SetLineStyle(kDashed);
    q2->Draw("HIST");
    
    q2_0->SetLineColor(kRed);
    q2_0->Draw("HIST SAME");
    
    q2_1->SetLineColor(kBlue);
    q2_1->Draw("HIST SAME");
    
    q2_2->SetLineColor(kMagenta);
    q2_2->Draw("HIST SAME");
    
    auto q2_num_in_nHCal_leg = new TLegend(0.65,0.75,0.85,0.9);
    q2_num_in_nHCal_leg->SetFillStyle(0);
    q2_num_in_nHCal_leg->AddEntry(q2, "All events", "l");
    q2_num_in_nHCal_leg->AddEntry(q2_0, "0-K in nHCAL", "l");
    q2_num_in_nHCal_leg->AddEntry(q2_1, "1-K in nHCAL", "l");
    q2_num_in_nHCal_leg->AddEntry(q2_2, "2-K in nHCAL", "l");
    q2_num_in_nHCal_leg->Draw();
    
    canvas17->Draw();
    canvas17->Print(filename17, "pdf");
    
    // end of q2 by num in nhcal histogram
    
    // FILE 18 - eta vs pT Scatterplot  //
        
        TString name18 = TString("eta_v_pT");
        TString filename18 = pdfdir + TString("/") + TString(name18) + TString(".pdf");
        TGraph *eta_v_pT = (TGraph*)ifile->Get("eta_v_pT");
        
        eta_v_pT->GetXaxis()->SetLimits(-4.5,4.5);
        eta_v_pT->SetMaximum(4);
        eta_v_pT->SetMarkerColor(kRed);
        eta_v_pT->SetMarkerStyle(20);
        eta_v_pT->SetMarkerSize(0.3);
        
        eta_v_pT->SetTitle("eta vs pT values");
        eta_v_pT->GetXaxis()->SetTitle("eta");
        eta_v_pT->GetYaxis()->SetTitle("pT");

        gStyle->SetOptStat(0);

        TCanvas *canvas18 = new TCanvas(name18, strang, 800, 700);
        canvas18->SetTopMargin(0.07);
        eta_v_pT->Draw("AP"); //Draw "A" - Axes, and "P" - Points(NOT lines)
        canvas18->Draw();
        canvas18->Print(filename18, "pdf");
         
    // End of FILE 18 - eta vs pT scatterplot //
    
    // FILE 19 - eta vs xB Scatterplot  //
        
        TString name19 = TString("eta_v_xB");
         
        TString filename19 = pdfdir + TString("/") + TString(name19) + TString(".pdf");
         
        TGraph *eta_v_xB = (TGraph*)ifile->Get("eta_v_xB");
        
        eta_v_xB->GetXaxis()->SetLimits(-4.5,4.5);
        eta_v_xB->GetYaxis()->SetLimits(0.,0.1);
        
        //Make the points appear as larger red circles
        eta_v_xB->SetMarkerColor(kRed);
        eta_v_xB->SetMarkerStyle(20);
        eta_v_xB->SetMarkerSize(0.3);
        
        //Add Title and Axis Labels
        eta_v_xB->SetTitle("eta vs xB values");
        eta_v_xB->GetXaxis()->SetTitle("eta");
        eta_v_xB->GetYaxis()->SetTitle("xB");

        gStyle->SetOptStat(0);

        TCanvas *canvas19 = new TCanvas(name19, strang, 800, 600);
        canvas19->SetTopMargin(0.07);

        eta_v_xB->Draw("AP"); //Draw "A" - Axes, and "P" - Points(NOT lines)
        canvas19->Draw();
        //canvas18->SetLogy();
        canvas19->Print(filename19, "pdf");
         
    // End of FILE 19 - eta vs xB scatterplot //
    
    // FILE 20 - eta vs Q2 Scatterplot  //
        
        TString name20 = TString("eta_v_q2");
         
        TString filename20 = pdfdir + TString("/") + TString(name20) + TString(".pdf");
         
        TGraph *eta_v_q2 = (TGraph*)ifile->Get("eta_v_q2");
        
    eta_v_q2->GetXaxis()->SetLimits(-4.5,4.5);
    eta_v_q2->GetYaxis()->SetLimits(0.,10);
        
        //Make the points appear as larger red circles
    eta_v_q2->SetMarkerColor(kRed);
    eta_v_q2->SetMarkerStyle(20);
    eta_v_q2->SetMarkerSize(0.3);
        
        //Add Title and Axis Labels
    eta_v_q2->SetTitle("eta vs Q2 values");
    eta_v_q2->GetXaxis()->SetTitle("eta");
    eta_v_q2->GetYaxis()->SetTitle("Q2");

        gStyle->SetOptStat(0);

        TCanvas *canvas20 = new TCanvas(name20, strang, 800, 600);
        canvas20->SetTopMargin(0.07);

    eta_v_q2->Draw("AP"); //Draw "A" - Axes, and "P" - Points(NOT lines)
        canvas20->Draw();
        //canvas18->SetLogy();
        canvas20->Print(filename20, "pdf");
         
    // End of FILE 20 - eta vs Q2 scatterplot //
    
    // FILE 21 - eta_v_pT HISTOGRAM //
     
    // Define the name of the plot:
    TString name21 = TString("eta_v_pT_hist");
     
    // Define the name of the pdf file:
    TString filename21 = pdfdir + TString("/") + TString(name21) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *eta_v_pT_hist = (TH2F*)ifile->Get("eta_v_pT_hist");
    eta_v_pT_hist->SetMarkerColor(kWhite);
    eta_v_pT_hist->SetMarkerSize(2);
    eta_v_pT_hist->GetXaxis()->SetTitle("eta");
    eta_v_pT_hist->GetYaxis()->SetTitle("pT");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas21 = new TCanvas(name21, strang, 700, 600);
    canvas21->SetLogy();
    canvas21->SetTopMargin(0.07);

    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    eta_v_pT_hist->Draw("COLZ");

    // Update the canvas
    canvas21->Update();
    canvas21->Draw();
    canvas21->Print(filename21, "pdf");

    // end of FILE 21 - eta_v_pT histogram
    
    // FILE 22 - eta_v_xB HISTOGRAM //
     
    // Define the name of the plot:
    TString name22 = TString("eta_v_xB_hist");
     
    // Define the name of the pdf file:
    TString filename22 = pdfdir + TString("/") + TString(name22) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *eta_v_xB_hist = (TH2F*)ifile->Get("eta_v_xB_hist");
    eta_v_xB_hist->SetMarkerColor(kWhite);
    eta_v_xB_hist->SetMarkerSize(2);
    eta_v_xB_hist->GetXaxis()->SetTitle("eta");
    eta_v_xB_hist->GetYaxis()->SetTitle("x_{Bj}");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas22 = new TCanvas(name22, strang, 700, 600);
    canvas22->SetTopMargin(0.07);

    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    eta_v_xB_hist->Draw("COLZ");

    // Update the canvas
    canvas22->Update();
    canvas22->SetLogy();
    canvas22->Draw();
    canvas22->Print(filename22, "pdf");

    // end of FILE 22 - eta_v_xB histogram
    
    // FILE 23 - eta_v_q2 HISTOGRAM //
     
    // Define the name of the plot:
    TString name23 = TString("eta_v_q2_hist");
     
    // Define the name of the pdf file:
    TString filename23 = pdfdir + TString("/") + TString(name23) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *eta_v_q2_hist = (TH2F*)ifile->Get("eta_v_q2_hist");
    eta_v_q2_hist->SetMarkerColor(kWhite);
    eta_v_q2_hist->SetMarkerSize(2);
    eta_v_q2_hist->GetXaxis()->SetTitle("eta");
    eta_v_q2_hist->GetYaxis()->SetTitle("Q^{2}");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas23 = new TCanvas(name23, strang, 700, 600);
    canvas23->SetTopMargin(0.07);

    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    eta_v_q2_hist->Draw("COLZ");

    // Update the canvas
    canvas23->Update();
    canvas23->SetLogy();
    canvas23->Draw();
    canvas23->Print(filename23, "pdf");

    // end of FILE 23 - eta_v_q2 histogram
    
    // start of FILE 24 - Calculated X_Bjorken Histogram
    
    TString name24 = TString("calcXb");
     
    TString filename24 = pdfdir + TString("/") + TString(name24) + TString(".pdf");
     
    TH1D *reco_calcXb = (TH1D*)ifile->Get("reco_calcXb");
    TH1D *gen_calcXb = (TH1D*)ifile->Get("gen_calcXb");
     
    gStyle->SetOptStat(0);

    TCanvas *canvas24 = new TCanvas(name24, strang, 800, 600);
    canvas24->SetLogx();
    canvas24->SetTopMargin(0.07);
    
    xBjorken->SetLineColor(kBlack);
    xBjorken->SetLineStyle(kDashed);
    xBjorken->Draw("HIST");
    
    reco_calcXb->SetLineColor(kBlue);
    reco_calcXb->Draw("HIST SAME");
    
    gen_calcXb->SetLineColor(kRed);
    gen_calcXb->Draw("HIST SAME");
    canvas24->Draw();
    
    auto leg24 = new TLegend(0.65,0.8,0.85,0.9); //x1,y1,x2,y2,header
    leg24->SetBorderSize(0);
    leg24->SetFillStyle(0);
    leg24->SetTextSize(0.03);
    leg24->AddEntry(xBjorken,"ROOT x_{Bj}","l");
    leg24->AddEntry(reco_calcXb,"reco calc x_{Bj}","l");
    leg24->AddEntry(gen_calcXb,"gen calc x_{Bj}","l");
    leg24->Draw();
    
    canvas24->Print(filename24, "pdf");

    // end of FILE 24 - Calculated X_Bjorken Histogram
    
    // start of FILE 25 - Calculated Q-squared Histogram
    
    TString name25 = TString("calcQ2");
     
    TString filename25 = pdfdir + TString("/") + TString(name25) + TString(".pdf");
     
    TH1D *reco_calcQ2 = (TH1D*)ifile->Get("reco_calcQ2");
    TH1D *gen_calcQ2 = (TH1D*)ifile->Get("gen_calcQ2");
     
    gStyle->SetOptStat(0);

    TCanvas *canvas25 = new TCanvas(name25, strang, 800, 600);
    canvas25->SetLogx();
    canvas25->SetTopMargin(0.07);
    
    q2->SetLineColor(kBlack);
    q2->SetLineStyle(kDashed);
    q2->Draw("HIST");
    
    reco_calcQ2->SetLineColor(kBlue);
    reco_calcQ2->Draw("HIST SAME");
    
    gen_calcQ2->SetLineColor(kRed);
    gen_calcQ2->Draw("HIST SAME");
    canvas25->Draw();
    
    auto leg25 = new TLegend(0.65,0.8,0.85,0.9); //x1,y1,x2,y2,header
    leg25->SetBorderSize(0);
    leg25->SetFillStyle(0);
    leg25->SetTextSize(0.03);
    leg25->AddEntry(q2,"ROOT Q^{2}","l");
    leg25->AddEntry(reco_calcQ2,"reco calc Q^{2}","l");
    leg25->AddEntry(gen_calcQ2,"gen calc Q^{2}","l");
    leg25->Draw();
    
    canvas25->Print(filename25, "pdf");

    // end of FILE 25 - Calculated Q-squared Histogram
    
    // start of FILE 26 - Reconstructed Calculated xBj vs Q2 Histogram
    
    TString name26 = TString("reco_xB_q2_hist");
     
    // Define the name of the pdf file:
    TString filename26 = pdfdir + TString("/") + TString(name26) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *reco_xB_q2_hist = (TH2F*)ifile->Get("reco_xB_q2_hist");
    reco_xB_q2_hist->SetMarkerColor(kWhite);
    reco_xB_q2_hist->SetMarkerSize(2);
    reco_xB_q2_hist->GetXaxis()->SetTitle("x_{Bj}");
    reco_xB_q2_hist->GetYaxis()->SetTitle("Q^{2}");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas26 = new TCanvas(name26, strang, 600, 600);
    canvas26->SetLogx();
    canvas26->SetLogy();
    canvas26->SetTopMargin(0.07);

    
    
    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    reco_xB_q2_hist->Draw("COLZ");

    // Update the canvas
    canvas26->Update();
    
    canvas26->Draw();
    canvas26->Print(filename26, "pdf");

    // end of FILE 26 - Reconstructed Calculated xBj vs Q2 Histogram
    
    // start of FILE 27 - Xbj % Error Histogram
    
    TString name27 = TString("calcXbj_Percent_Error");
     
    TString filename27 = pdfdir + TString("/") + TString(name27) + TString(".pdf");
     
    TH1D *reco_calcXb_PercentError = (TH1D*)ifile->Get("reco_calcXb_PercentError");
    TH1D *gen_calcXb_PercentError = (TH1D*)ifile->Get("gen_calcXb_PercentError");
    
    gStyle->SetOptStat(0);

    TCanvas *canvas27 = new TCanvas(name27, strang, 800, 600);
    canvas27->SetTopMargin(0.07);
    
    reco_calcXb_PercentError->SetLineColor(kBlue);
    gen_calcXb_PercentError->SetLineColor(kRed);
    
    gen_calcXb_PercentError->Draw("HIST");
    reco_calcXb_PercentError->Draw("HIST SAME");
    canvas27->Draw();
    
    auto leg27 = new TLegend(0.65,0.8,0.85,0.9); //x1,y1,x2,y2,header
    leg27->SetBorderSize(0);
    leg27->SetFillStyle(0);
    leg27->SetTextSize(0.03);
    leg27->AddEntry(reco_calcXb_PercentError,"reco calc x_{Bj}","l");
    leg27->AddEntry(gen_calcXb_PercentError,"gen calc x_{Bj}","l");
    leg27->Draw();
    
    canvas27->Print(filename27, "pdf");

    // end of FILE 27 - Xbj % Error Histogram
    
    // start of FILE 28 - Q2 % Error Histogram
    
    TString name28 = TString("calcQ2_Percent_Error");
     
    TString filename28 = pdfdir + TString("/") + TString(name28) + TString(".pdf");
     
    TH1D *reco_calcQ2_PercentError = (TH1D*)ifile->Get("reco_calcQ2_PercentError");
    TH1D *gen_calcQ2_PercentError = (TH1D*)ifile->Get("gen_calcQ2_PercentError");
    
    gStyle->SetOptStat(0);

    TCanvas *canvas28 = new TCanvas(name28, strang, 800, 600);
    canvas28->SetTopMargin(0.07);
    
    reco_calcQ2_PercentError->SetLineColor(kBlue);
    gen_calcQ2_PercentError->SetLineColor(kRed);
    
    gen_calcQ2_PercentError->Draw("HIST");
    reco_calcQ2_PercentError->Draw("HIST SAME");
    
    auto leg28 = new TLegend(0.65,0.8,0.85,0.9); //x1,y1,x2,y2,header
    leg28->SetBorderSize(0);
    leg28->SetFillStyle(0);
    leg28->SetTextSize(0.03);
    leg28->AddEntry(reco_calcXb_PercentError,"reco calc Q^{2}","l");
    leg28->AddEntry(gen_calcXb_PercentError,"gen calc Q^{2}","l");
    leg28->Draw();
    
    canvas28->Draw();
    
    canvas28->Print(filename28, "pdf");

    // end of FILE 28 - Q2 % Error Histogram
    
    // start of FILE 29 - Generated Calculated xBj vs Q2 Histogram
    
    TString name29 = TString("gen_xB_q2_hist");
     
    // Define the name of the pdf file:
    TString filename29 = pdfdir + TString("/") + TString(name29) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *gen_xB_q2_hist = (TH2F*)ifile->Get("gen_xB_q2_hist");
    gen_xB_q2_hist->SetMarkerColor(kWhite);
    gen_xB_q2_hist->SetMarkerSize(2);
    gen_xB_q2_hist->GetXaxis()->SetTitle("x_{Bj}");
    gen_xB_q2_hist->GetYaxis()->SetTitle("Q^{2}");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas29 = new TCanvas(name29, strang, 600, 600);
    canvas29->SetLogx();
    canvas29->SetLogy();
    canvas29->SetTopMargin(0.07);

    
    
    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    gen_xB_q2_hist->Draw("COLZ");

    // Update the canvas
    canvas29->Update();
    
    canvas29->Draw();
    canvas29->Print(filename29, "pdf");

    // end of FILE 29 - Generated Calculated xBj vs Q2 Histogram
    
    // start of FILE 30 - Generated vs Reconstructed Calculated xBj Histogram
    
    TString name30 = TString("gen_reco_xB_hist");
     
    // Define the name of the pdf file:
    TString filename30 = pdfdir + TString("/") + TString(name30) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *gen_reco_xB_hist = (TH2F*)ifile->Get("gen_reco_xB_hist");
    gen_reco_xB_hist->SetMarkerColor(kWhite);
    gen_reco_xB_hist->SetMarkerSize(2);
    gen_reco_xB_hist->GetXaxis()->SetTitle("Generated Calc x_{Bj}");
    gen_reco_xB_hist->GetYaxis()->SetTitle("Reconstructed Calc x_{Bj}");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas30 = new TCanvas(name30, strang, 600, 600);
    canvas30->SetLogx();
    canvas30->SetLogy();
    canvas30->SetTopMargin(0.07);
    canvas30->SetLeftMargin(0.13);
    
    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    gen_reco_xB_hist->Draw("COLZ");

    // Update the canvas
    canvas30->Update();
    
    canvas30->Draw();
    canvas30->Print(filename30, "pdf");

    // end of FILE 30 - Generated vs Reconstructed Calculated xBj Histogram
    
    // start of FILE 31 - Generated vs Reconstructed Calculated xBj Histogram
    
    TString name31 = TString("gen_reco_q2_hist");
     
    // Define the name of the pdf file:
    TString filename31 = pdfdir + TString("/") + TString(name31) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *gen_reco_q2_hist = (TH2F*)ifile->Get("gen_reco_q2_hist");
    gen_reco_q2_hist->SetMarkerColor(kWhite);
    gen_reco_q2_hist->SetMarkerSize(2);
    gen_reco_q2_hist->GetXaxis()->SetTitle("Generated Calc Q^{2}");
    gen_reco_q2_hist->GetYaxis()->SetTitle("Reconstructed Calc Q^{2}");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas31 = new TCanvas(name31, strang, 600, 600);
    canvas31->SetLogx();
    canvas31->SetLogy();
    canvas31->SetTopMargin(0.07);
    
    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    gen_reco_q2_hist->Draw("COLZ");

    // Update the canvas
    canvas31->Update();
    
    canvas31->Draw();
    canvas31->Print(filename31, "pdf");

    // end of FILE 31 - Generated vs Reconstructed Calculated xBj Histogram
    
    // start of FILE 32 - Generated Calc vs ROOT xBj Histogram
    
    TString name32 = TString("gen_root_xB_hist");
     
    // Define the name of the pdf file:
    TString filename32 = pdfdir + TString("/") + TString(name32) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *gen_root_xB_hist = (TH2F*)ifile->Get("gen_root_xB_hist");
    gen_root_xB_hist->SetMarkerColor(kWhite);
    gen_root_xB_hist->SetMarkerSize(2);
    gen_root_xB_hist->GetXaxis()->SetTitle("Generated Calc x_{Bj}");
    gen_root_xB_hist->GetYaxis()->SetTitle("ROOT Prodived x_{Bj}");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas32 = new TCanvas(name32, strang, 600, 600);
    canvas32->SetLogx();
    canvas32->SetLogy();
    canvas32->SetTopMargin(0.07);
    
    canvas32->SetRightMargin(0.13);
    canvas32->SetLeftMargin(0.12);
    
    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    gen_root_xB_hist->Draw("COLZ");

    // Update the canvas
    canvas32->Update();
    
    canvas32->Draw();
    canvas32->Print(filename32, "pdf");

    // end of FILE 32 - Generated Calc vs ROOT xBj Histogram
    
    // start of FILE 33 - Reconstructed Calc vs ROOT xBj Histogram
    
    TString name33 = TString("reco_root_xB_hist");
     
    // Define the name of the pdf file:
    TString filename33 = pdfdir + TString("/") + TString(name33) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *reco_root_xB_hist = (TH2F*)ifile->Get("reco_root_xB_hist");
    reco_root_xB_hist->SetMarkerColor(kWhite);
    reco_root_xB_hist->SetMarkerSize(2);
    reco_root_xB_hist->GetXaxis()->SetTitle("Reconstructed Calc x_{Bj}");
    reco_root_xB_hist->GetYaxis()->SetTitle("ROOT Provided x_{Bj}");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas33 = new TCanvas(name33, strang, 600, 600);
    canvas33->SetLogx();
    canvas33->SetLogy();
    canvas33->SetTopMargin(0.07);
    
    canvas33->SetRightMargin(0.1);
    canvas33->SetLeftMargin(0.15);
    
    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    reco_root_xB_hist->Draw("COLZ");

    // Update the canvas
    canvas33->Update();
    
    canvas33->Draw();
    canvas33->Print(filename33, "pdf");

    // end of FILE 33 - Reconstructed Calc vs ROOT xBj Histogram
    
    // start of FILE 34 - ALL GEN/RECO/ROOT DIFFS
    
    TString name34 = TString("gen_reco_xBj_q2_diffs");
    // Define the name of the pdf file:
    TString filename34 = pdfdir + TString("/") + TString(name34) + TString(".pdf");

    TCanvas *canvas34 = new TCanvas(name34, strang, 1000, 600);
    
    TH1F *xBj_gen_reco_calc_diff = (TH1F*)ifile->Get("xBj_gen_reco_calc_diff");
    TH1F *q2_gen_reco_calc_diff = (TH1F*)ifile->Get("q2_gen_reco_calc_diff");
    TH1F *reco_calcXb_Diff = (TH1F*)ifile->Get("reco_calcXb_Diff");
    TH1F *reco_calcQ2_Diff = (TH1F*)ifile->Get("reco_calcQ2_Diff");
    TH1F *gen_calcXb_Diff = (TH1F*)ifile->Get("gen_calcXb_Diff");
    TH1F *gen_calcQ2_Diff = (TH1F*)ifile->Get("gen_calcQ2_Diff");
    
    //devide the canvas into several pads (x,y):
    canvas34->Divide(2,2);
    
    // go to the first one and do your thing
    canvas34->cd(1);
    xBj_gen_reco_calc_diff->Draw();
    
    // pad 2
    canvas34->cd(2);
    q2_gen_reco_calc_diff->Draw();

    // pad 3
    canvas34->cd(3);
    reco_calcXb_Diff->SetLineColor(kBlue);
    reco_calcXb_Diff->SetMaximum(7500);
    reco_calcXb_Diff->Draw();
    gen_calcXb_Diff->SetLineColor(kRed);
    gen_calcXb_Diff->Draw("SAME");

    auto leg_xb_diffs = new TLegend(0.65,0.8,0.85,0.9); //x1,y1,x2,y2,header
    leg_xb_diffs->SetBorderSize(0);
    leg_xb_diffs->SetTextSize(0.05);
    leg_xb_diffs->SetFillStyle(0);
    leg_xb_diffs->AddEntry(reco_calcXb_Diff,"reco calc x_{Bj}","l");
    leg_xb_diffs->AddEntry(gen_calcXb_Diff,"gen calc x_{Bj}","l");
    leg_xb_diffs->Draw();
    
    // pad 4
    canvas34->cd(4);
    reco_calcQ2_Diff->SetLineColor(kBlue);
    reco_calcQ2_Diff->SetMaximum(6000);
    reco_calcQ2_Diff->Draw();
    gen_calcQ2_Diff->SetLineColor(kRed);
    gen_calcQ2_Diff->Draw("SAME");

    auto leg_q2_diffs = new TLegend(0.65,0.8,0.85,0.9); //x1,y1,x2,y2,header
    leg_q2_diffs->SetBorderSize(0);
    leg_q2_diffs->SetTextSize(0.05);
    leg_q2_diffs->SetFillStyle(0);
    leg_q2_diffs->AddEntry(reco_calcQ2_Diff,"reco calc x_{Bj}","l");
    leg_q2_diffs->AddEntry(gen_calcQ2_Diff,"gen calc x_{Bj}","l");
    leg_q2_diffs->Draw();
    
    // draw canvas and print to pdf
    canvas34->Draw();
    canvas34->Print(filename34, "pdf");
    
    // end of FILE 34 - ALL GEN/RECO/ROOT DIFFS

    // FILE 35 - K1_K2_eta_hist HISTOGRAM //
     
    // Define the name of the plot:
    TString name35 = TString("K1_K2_eta_hist");
     
    // Define the name of the pdf file:
    TString filename35 = pdfdir + TString("/") + TString(name35) + TString(".pdf");
     
    //Get the histogram data from ePIC_Analysis.C
    TH2F *K1_K2_eta_hist = (TH2F*)ifile->Get("K1_K2_eta_hist");
    K1_K2_eta_hist->SetMarkerColor(kWhite);
    K1_K2_eta_hist->SetMarkerSize(2);
    K1_K2_eta_hist->GetXaxis()->SetTitle("#eta_{K1}");
    K1_K2_eta_hist->GetYaxis()->SetTitle("#eta_{K2}");
    
     //gStyle->SetOptStat(0);
    TCanvas *canvas35 = new TCanvas(name35, strang, 700, 600);
    canvas35->SetTopMargin(0.07);

    gStyle->SetPalette(kBird);// Choose color palette
    gStyle->SetNumberContours(256);

    K1_K2_eta_hist->Draw("COLZ");

    // Update the canvas
    canvas35->Update();
    canvas35->Draw();
    canvas35->Print(filename35, "pdf");

    // end of FILE 35 - K1_K2_eta_hist HISTOGRAM
        
  //
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
  
} // end of macro  
//////////////////

  //TPaveText *t = new TPaveText(.05,.3,.95,.6, "NDC");                                                                                   
  //t->AddText("This line is blue"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);                                          
  //t->AddText("This line is red");  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kRed);                                           
  //t->Draw();     
