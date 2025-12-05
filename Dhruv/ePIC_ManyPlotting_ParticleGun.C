// Run this macro, from the Linux / Terminal command line, as
// root -l -q ePIC_Plotting.C
// The only condition is that previously ran ePIC_Analysis.C  with the identical "strang" as hardcoded in this macro (below)
// Strang must match the string in a file starting with out*, in this directory
// CR 2024-08-14/20

void ePIC_ManyPlotting_ParticleGun()
{
  gSystem->Exec("date");
  
    int nlayers_nhcal = 10;
    std::vector<int> momentum_vals = {1, 3, 5, 7, 10};
    
    const char flavor[]="ePIC";

    TString filetype_dir = TString("ParticleGun/");
    TString pdfdir = TString("graphs/ParticleGun/ManyPlots/");
    
    //=======================================================
    
    std::vector<TFile*> ifiles;
    
    for (int i = 0; i < momentum_vals.size(); i++) {
        std::cout << "MOMENTUM:" << momentum_vals[i] << "\n";
        TString strang = Form("Muon_%d_GeV", momentum_vals[i]);
        std::cout << strang<< " AND " << filetype_dir << "\n";
        TString infile_ram = TString("out_files/") + filetype_dir + strang + "/" + "out." + strang + TString("-") + flavor + TString(".root");
        ifiles.push_back(TFile::Open(infile_ram.Data(), "READ"));
        std::cout << infile_ram << " THE FILENAME\n";
    }
    
    // Create canvas
    TCanvas *c = new TCanvas("c", "Mean Energy per Layer", 900, 700);

    // Prepare multigraph
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Mean Hit Energy vs Layer;Layer Index;Mean Energy");

    // Color palette for 5 momenta
    int colors[5] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7};

    // Create one TGraph per momentum file
    for (int i = 0; i < ifiles.size(); i++) {

        std::vector<double> layer_idx(nlayers_nhcal);
        std::vector<double> layer_mean(nlayers_nhcal);
        std::vector<double> layer_err(nlayers_nhcal);

        for (int L = 0; L < nlayers_nhcal; L++) {
            TString histname = Form("nhCal_LayerE_%d", L);
            TH1F *h = (TH1F*) ifiles[i]->Get(histname);

            if (!h) {
                std::cout << "Missing hist: " << histname
                          << " in file " << i << std::endl;
                layer_mean[L] = 0;
                layer_idx[L] = L;
                continue;
            }
            layer_idx[L] = L;
            layer_mean[L] = h->GetMean();
            
            int n = h->GetEntries();
            cout << "SIZE IS " << n << "\n";
            if (n > 0)
                layer_err[L] = 1.0 / std::sqrt(n);
            else
                layer_err[L] = 0.0;
        }

        // Create the line for this momentum
        TGraphErrors *gr = new TGraphErrors(nlayers_nhcal,
                                layer_idx.data(),
                                layer_mean.data(), nullptr, layer_err.data());

        gr->SetLineColor(colors[i]);
        gr->SetMarkerColor(colors[i]);
        gr->SetMarkerStyle(20);
        gr->SetLineWidth(2);
        gr->SetTitle(Form("p = %d GeV", momentum_vals[i]));

        mg->Add(gr, "LP");
    }

    // Draw result
    mg->Draw("A");

    c->BuildLegend();
    c->SetGrid();
    c->Print(pdfdir + TString("nHCal_byLayer_Momentum_E.pdf"));
    // Save PDF
    
    //=======================================================
    
    // Create canvas
    TCanvas *c2 = new TCanvas("c2", "Mean Energy per Layer", 900, 700);

    // Prepare multigraph
    TMultiGraph *mg2 = new TMultiGraph();
    mg2->SetTitle("Energy Deposit Per Layer By Momentum");

    // Color palette for 5 momenta
    int colors2[10] = {kRed, kOrange-3, kYellow+1, kSpring, kGreen+2, kAzure+10, kBlue, kViolet-6, kMagenta, kBlack};
    // Create one TGraph per layer
    for (int i = 0; i < nlayers_nhcal; i++) {

        std::vector<double> momentum_idx(5);
        std::vector<double> layer_mean(nlayers_nhcal);
        std::vector<double> mom_err(nlayers_nhcal);

        for (int L = 0; L < momentum_vals.size(); L++) {
            TString histname = Form("nhCal_LayerE_%d", L);
            TH1F *h = (TH1F*) ifiles[L]->Get(Form("nhCal_LayerE_%d", i));

            if (!h) {
                std::cout << "Missing hist: " << histname
                          << " in file " << L << std::endl;
                layer_mean[L] = 0;
                momentum_idx[L] = momentum_vals[L];
                continue;
            }
            
            momentum_idx[L] = momentum_vals[L];
            layer_mean[L] = h->GetMean();
            
            int n = h->GetEntries();
            cout << "SIZE IS " << n << "\n";
            if (n > 0)
                mom_err[L] = 1.0 / std::sqrt(n);
            else
                mom_err[L] = 0.0;
        }

        // Create the line for this momentum
        TGraphErrors *gr2 = new TGraphErrors(momentum_vals.size(),
                                momentum_idx.data(),
                                layer_mean.data(), nullptr, mom_err.data());

        gr2->SetLineColor(colors2[i]);
        gr2->SetMarkerColor(colors2[i]);
        gr2->SetMarkerStyle(20);
        gr2->SetLineWidth(2);
        gr2->SetTitle(Form("Layer %i", i));

        mg2->Add(gr2, "LP");
    }

    // Draw result
    mg2->Draw("A");

    c2->BuildLegend();
    c2->SetGrid();
    c2->Print(pdfdir + TString("nHCal_byMomentum_Layer_E.pdf"));
    // Save PDF
    
    
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
  
} // end of macro  
//////////////////
