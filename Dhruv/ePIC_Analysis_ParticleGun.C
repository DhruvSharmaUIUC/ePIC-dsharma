// Run this macro, from the Linux / Terminal command line, as
// root -l -q ePIC_Analysis.C/Users/dhruv/Desktop/ePIC/Dhruv/ePIC_Plotting.C
// The only condition is that there is a subdirectory called "data" that contains the input MC file
// (provided by Caroline on Box / copied from SDCC)
// The name of this input MC file (variable "strang") is hardcoded in this macro and must match with the file name
// CR 2024-08-14/20

//Defines constants and collection of Calorimeters/acceptances
const char* nHCal_name = "nHCal";
const char* Barrel_name = "Barrel";
const char* LFHCAL_name = "LFHCAL";
const char *cals[] = {nHCal_name, Barrel_name, LFHCAL_name};
double cal_limits[sizeof(cals)][2] = {{-4.05, -1.2}, {-1.2, 1.18}, {1.18, 4.2}};
const double elec_BeamE = 18; // electron and proton beam energies (GeV)
const double prot_BeamE = 110;
const double elec_mom_init = -17.9999999927; // initial momentum magnitudes of electron and proton
const double prot_mom_init = 109.996000636;

//Array containing the number of phi mesons with 0, 1, or 2 kaons for each calorimeter
float calNums[sizeof(cals)][3];

//Array containing the K1/K2 crossover within different calorimeters and an "All"
float calMatrix[4][4];

void ConvertToLogBins2D(TH2F*& hist, int numBinsX, double minX, double maxX, int numBinsY, double minY, double maxY) {
    // Check if minX and minY are positive
    if (minX <= 0 || minY <= 0) {
        std::cerr << hist->GetName() << " Error: Minimum values for logarithmic bins must be positive." << std::endl;
        return;
    }

    // Prepare logarithmic bin edges for x-axis
    std::vector<float> binEdgesX(numBinsX + 1);
        double logMinX = TMath::Log10(minX);
        double logMaxX = TMath::Log10(maxX);
        double logBinWidthX = (logMaxX - logMinX) / numBinsX;
        for (int i = 0; i <= numBinsX; ++i) {
            binEdgesX[i] = TMath::Power(10, logMinX + i * logBinWidthX);
        }

    // Prepare logarithmic bin edges for y-axis
    std::vector<float> binEdgesY(numBinsY + 1);
    double logMinY = TMath::Log10(minY);
    double logMaxY = TMath::Log10(maxY);
    double logBinWidthY = (logMaxY - logMinY) / numBinsY;
    for (int i = 0; i <= numBinsY; ++i) {
        binEdgesY[i] = TMath::Power(10, logMinY + i * logBinWidthY);
    }
    
    // Create a new 2D histogram with logarithmic binning on both axes
    TH2F* logHist2D = new TH2F(hist->GetName(), hist->GetTitle(), numBinsX, binEdgesX.data(), numBinsY, binEdgesY.data());

    delete hist;      // Delete the old histogram to free memory
    hist = logHist2D; // Update the pointer to point to the new histogram
}

void ConvertToLogYBins2D(TH2F*& histY, int numBinsX, double minX, double maxX, int numBinsY, double minY, double maxY) {
    if (minY <= 0) {
        std::cerr << histY->Class_Name() << "Error: Minimum values for logarithmic bins must be positive." << std::endl;
        return;
    }
    
    std::vector<float> binEdgesX(numBinsX + 1);
    double binWidthX = (maxX - minX) / numBinsX;
    for (int i = 0; i <= numBinsX; ++i) {
        binEdgesX[i] = minX + i * binWidthX;
    }
    
    std::vector<float> binEdgesY(numBinsY + 1);
    double logMinY = TMath::Log10(minY);
    double logMaxY = TMath::Log10(maxY);
    double logBinWidthY = (logMaxY - logMinY) / numBinsY;
    for (int i = 0; i <= numBinsY; ++i) {
        binEdgesY[i] = TMath::Power(10, logMinY + i * logBinWidthY);
    }
    
    TH2F* logYHist2D = new TH2F(histY->GetName(), histY->GetTitle(), numBinsX, binEdgesX.data(), numBinsY, binEdgesY.data());
    
    delete histY;
    histY = logYHist2D;
}

void ConvertToLogBins(TH1F*& hist, int numBins, double min, double max) {
    // Check if min is positive
    if (min <= 0) {
        std::cerr << "Error: Log bins can only take positive values" << std::endl;
        return;
    }

    std::vector<double> binEdges(numBins + 1);

    // Calculate log bin edges
    double logMin = TMath::Log10(min);
    double logMax = TMath::Log10(max);
    double logBinWidth = (logMax - logMin) / numBins;

    for (int i = 0; i <= numBins; ++i) {
        binEdges[i] = TMath::Power(10, logMin + i * logBinWidth);
    }

    // Create a new histogram with logarithmic binning
    TH1F* logHist = new TH1F("logHist", hist->GetTitle(), numBins, binEdges.data());

    // Transfer contents from the original histogram to the new one
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double x = hist->GetBinCenter(i);
        double y = hist->GetBinContent(i);
        //double error = hist->GetBinError(i);
        
        // Find the corresponding bin in the logarithmic histogram
        int bin = logHist->FindBin(x);
        if (bin >= 1 && bin <= numBins) {
            logHist->SetBinContent(bin, y + logHist->GetBinContent(bin));
            //logHist->SetBinError(bin, TMath::Sqrt(TMath::Power(logHist->GetBinError(bin), 2) + TMath::Power(error, 2)));
        }
    }

    // Replace the original histogram with the logarithmic one
    hist->SetBins(numBins, binEdges.data());
    for (int i = 1; i <= numBins; ++i) {
        hist->SetBinContent(i, logHist->GetBinContent(i));
        //hist->SetBinError(i, logHist->GetBinError(i));
    }

    // Clean up
    delete logHist;
}

double calcPortionErr(double part, double total) { //Manually calculate errors for ratio values part/total
    return (part/total) * TMath::Sqrt(TMath::Power(TMath::Sqrt(part) / part, 2) +
                                      TMath::Power(TMath::Sqrt(total) / total, 2));
}
//Accepts a string cal_name containing the desired calorimeter name and returns true if an eta particle_eta is within acceptance for the given calorimeter
bool in_Cal_Tolerance(const char* cal_name, float particle_eta) {
    if (cal_name == nHCal_name) {
        return (particle_eta >= -4.05 && particle_eta < -1.2);
    } else if (cal_name == Barrel_name) {
        return (particle_eta >= -1.2 && particle_eta < 1.18);
    } else if (cal_name == LFHCAL_name) {
        return (particle_eta >= 1.18 && particle_eta < 4.2);
    }
    return false;
}

//Accepts two floats, one for each kaon's eta and fills calNums according to if zero, one, or two kaons are within tolerance of each calorimeter.
void kaons_in_Cal(float p1_eta, float p2_eta) {
    for (int i = 0; i < sizeof(cals) / sizeof(cals[0]); i++) {
        int particles_in_cal = 0;
        
        if (in_Cal_Tolerance(cals[i], p1_eta)) {
            particles_in_cal++;
        }
        if (in_Cal_Tolerance(cals[i], p2_eta)) {
            particles_in_cal++;
        }
        
        calNums[i][particles_in_cal]++;
    }
}

//Accepts two floats representing each kaon's eta and fills calMatrix with the intersection of phi decay, showing how many phi mesons had decay in a given combination of calorimeters, or fell into any calorimeter acceptance
void fill_Cal_Arr(float p1_eta, float p2_eta) {
    int p1_Cal_index = -1;
    int p2_Cal_index = -1;
    const char *graph_cals[] = {LFHCAL_name, Barrel_name, nHCal_name};
    
    for (int i = 0; i < sizeof(graph_cals); i++) {
        if (in_Cal_Tolerance(graph_cals[i], p1_eta)) {
            p1_Cal_index = i;
        }
        if (in_Cal_Tolerance(graph_cals[i], p2_eta)) {
            p2_Cal_index = i;
        }
    }
    if (p1_Cal_index != -1 && p2_Cal_index != -1) {
        calMatrix[p1_Cal_index+1][3-(p2_Cal_index+1)]++;
        calMatrix[0][3-(p2_Cal_index+1)]++;
        calMatrix[p1_Cal_index+1][3]++;
        calMatrix[0][3]++;
    }
}

//Create a class phiDecay representing each decay and contain attributes of each decay in a single unit for access throughout the program
class phiDecay {
public:
    float eta1, eta2, k1_pT, k2_pT, genEta1, genEta2;
    double x_b, q2, reco_calcXb, reco_calcQ2, gen_calcXb, gen_calcQ2;
    double reco_elec_mom_x_f, reco_elec_mom_y_f, reco_elec_mom_z_f, reco_elec_energy_f;
    double gen_elec_mom_x_f, gen_elec_mom_y_f, gen_elec_mom_z_f, gen_elec_energy_f;
    int num_in_nHCal;
    int reco_parts;
    bool k1_reco = false, k2_reco = false, reco_electron = false;
    
    phiDecay(float e1, float e2, double xb_val, double q2_val) {
        eta1 = e1;
        eta2 = e2;
        x_b = xb_val;
        q2 = q2_val;
        
        
        int particles_in_cal = 0;
        
        if (in_Cal_Tolerance("nHCal", eta1)) {
            particles_in_cal++; }
        if (in_Cal_Tolerance("nHCal", eta2)) {
            particles_in_cal++; }
        
        num_in_nHCal = particles_in_cal;
    }
    
    phiDecay() {
        
    }
    
    void set_num_in_nHCal() {
        int particles_in_cal = 0;
        
        if (in_Cal_Tolerance("nHCal", eta1)) {
            particles_in_cal++; }
        if (in_Cal_Tolerance("nHCal", eta2)) {
            particles_in_cal++; }
        
        num_in_nHCal = particles_in_cal;
    }
    
    void calcRecoXbQ2() {
        reco_calcQ2 = -(((elec_BeamE - reco_elec_energy_f) * (elec_BeamE - reco_elec_energy_f)) - (reco_elec_mom_x_f * reco_elec_mom_x_f) - (reco_elec_mom_y_f * reco_elec_mom_y_f) - ((elec_mom_init - reco_elec_mom_z_f) * (elec_mom_init - reco_elec_mom_z_f)));
        
        reco_calcXb = reco_calcQ2/(2 * ((prot_BeamE * (elec_BeamE - reco_elec_energy_f)) - prot_mom_init * (elec_mom_init - reco_elec_mom_z_f)));
    }
    
    void calcGenXbQ2() {
        gen_elec_energy_f = sqrt((0.000511 * 0.000511) + (gen_elec_mom_x_f * gen_elec_mom_x_f) + (gen_elec_mom_y_f * gen_elec_mom_y_f) + (gen_elec_mom_z_f * gen_elec_mom_z_f));
        
        gen_calcQ2 = -(((elec_BeamE - gen_elec_energy_f) * (elec_BeamE - gen_elec_energy_f)) - (gen_elec_mom_x_f * gen_elec_mom_x_f) - (gen_elec_mom_y_f * gen_elec_mom_y_f) - ((elec_mom_init - gen_elec_mom_z_f) * (elec_mom_init - gen_elec_mom_z_f)));
        
        gen_calcXb = gen_calcQ2/(2 * ((prot_BeamE * (elec_BeamE - gen_elec_energy_f)) - prot_mom_init * (elec_mom_init - gen_elec_mom_z_f)));
    }
};

//Array to contain all decays
std::vector<phiDecay> decays;

void ePIC_Analysis_ParticleGun(int gevnum = 10, char particle_name[] = "Muon"){

  gSystem->Exec("date");
  //
    int nlayers_nhcal=10; //currently hardcoded
    
    const double hx_min_nhcal = -2800.;
    const double hx_max_nhcal = 2800.;
    const double hy_min_nhcal = -2800.;
    const double hy_max_nhcal = 2800.;
    const double hz_min_nhcal = -4500.;
    const double hz_max_nhcal = -3950.;
    
    const double z_min_nhcal = -395.; // (2024-12-03) start of nHCal in z-direction [cm]
    const double z_thickness_nhcal = 45.; // nHCal thickness in z [cm]
    const double z_max_nhcal = z_min_nhcal - z_thickness_nhcal;
    
    std::map<double,int> layer_zs;
    
    for (int i = 0; i < nlayers_nhcal; i++) {
        layer_zs.insert({10 * (z_min_nhcal - z_thickness_nhcal + 0.25) + (i * 45), i+1});
    }
    

  ////////////////////////////////////////////////////
  //// String definitions - modify here as needed ////
  /////////////////////////////////////////////////////
    cout << "particle name is " << particle_name << "\n";
  // Define input directory if reading locally:
  const char indir[]="data";
    //TString filetype_dir = Form("ParticleGun/%s_%d_GeV/", particle_name, gevnum);
    TString filetype_dir = Form("ParticleGun/%s_%d_GeV_Uniform/", particle_name, gevnum);
  cout << "Input directory is: " << indir << "/" << filetype_dir << " \n";
    
  // Define name of local input MC file:
  //const char strang[]="sartre_bnonsat_Au_phi_ab_eAu_1.0000.eicrecon.tree.edm4eic"; // "strang = data string"
    TString strang = Form("%s_%d_GeV", particle_name,gevnum);  //const char strang[]="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.0998.eicrecon.tree.edm4eic"; // "strang = data string"

    //
    //TString runlist_ram = Form("local_runlists/%s_%d_GeV_runlist.txt", particle_name, gevnum);  // MULTI FILE RUNNING
    TString runlist_ram = Form("local_runlists/%s_%d_GeV_Uniform_runlist.txt", particle_name, gevnum);
    const char *runlist=runlist_ram.Data(); // MULTI FILE RUNNING
    cout << "RUNLIST = " << runlist_ram.Data() << "\n";
  // define flavor of this macro:
  const char flavor[]="ePIC";
  
  ////////////////////////////////////
  //// end of string definitions ////
  ///////////////////////////////////

  // If reading from a locally copied file:
  TString infile_ram=indir + TString("/") + filetype_dir + strang + TString(".root");  // remains in RAM
  const char *infile=infile_ram.Data(); // points to a valid address in RAM
  cout << "Analyzed MC file will be: " << infile << " \n";
  

  //TString outfile_ram= TString("out_files/") + filetype_dir + "out." + strang + TString("-") + flavor + TString(".root");
    TString outfile_ram= TString("out_files/") + filetype_dir + "out." + strang + TString("_Uniform-") + flavor + TString(".root");
  const char *outfile=outfile_ram.Data();
  TFile *ofile = TFile::Open(outfile,"RECREATE"); // RECREATE overwrites an existing file of the same name
  
  TChain *mychain = new TChain("events");

  // if reading a single file:
  //mychain->Add(infile); //SINGLE FILE RUNNING
    
  std::ifstream in(runlist); //MULTI FILE RUNNING
  std::string file("");  //MULTI FILE RUNNING
  while (in >> file) mychain->Add(file.data());  //MULTI FILE RUNNING
  ///////////////////////////////////////
  //// end of automated definitions ////
  //////////////////////////////////////
  // Initialize reader
  TTreeReader tree_reader(mychain);
    
  //Get Event-Level Information
  //Get x_bjorken And q^2 Information
  TTreeReaderArray<float> partXb(tree_reader, "InclusiveKinematicsTruth.x");
  TTreeReaderArray<float> partQ2(tree_reader, "InclusiveKinematicsTruth.Q2");
    
  // Get Particle Information
  TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
  TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
  TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
  TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
  //TTreeReaderArray<float> partEnergy(tree_reader, "MCParticles.energy");
  TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");

  // Get Reconstructed Track Information
  TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
  TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
  TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
  TTreeReaderArray<float> trackPartEnergy(tree_reader, "ReconstructedChargedParticles.energy");
    
  // Get Associations Between MCParticles and ReconstructedChargedParticles
  TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
  TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");

  // Get parent and daugther information
  TTreeReaderArray<int> parents_index(tree_reader, "_MCParticles_parents.index");
  TTreeReaderArray<unsigned int> parents_begin(tree_reader, "MCParticles.parents_begin");
  TTreeReaderArray<unsigned int> parents_end(tree_reader, "MCParticles.parents_end");
  
  TTreeReaderArray<int> daughters_index(tree_reader, "_MCParticles_daughters.index");
  TTreeReaderArray<unsigned int> daughters_begin(tree_reader, "MCParticles.daughters_begin");
  TTreeReaderArray<unsigned int> daughters_end(tree_reader, "MCParticles.daughters_end");
  
    TTreeReaderArray<unsigned int> nHCalClustersnHits(tree_reader, "HcalEndcapNClusters.nhits");
    TTreeReaderArray<float> nHCalClustersE(tree_reader, "HcalEndcapNClusters.energy");
    TTreeReaderArray<float> nHCalClustersPosX(tree_reader, "HcalEndcapNClusters.position.x");
    TTreeReaderArray<float> nHCalClustersPosY(tree_reader, "HcalEndcapNClusters.position.y");
    TTreeReaderArray<float> nHCalClustersPosZ(tree_reader, "HcalEndcapNClusters.position.z");
    
    TTreeReaderArray<float> nHCalRecHitsE(tree_reader, "HcalEndcapNRecHits.energy");
    TTreeReaderArray<float> nHCalRecHitsPosX(tree_reader, "HcalEndcapNRecHits.position.x");
    TTreeReaderArray<float> nHCalRecHitsPosY(tree_reader, "HcalEndcapNRecHits.position.y");
    TTreeReaderArray<float> nHCalRecHitsPosZ(tree_reader, "HcalEndcapNRecHits.position.z");
    
    TTreeReaderArray<unsigned int> recoAssocClusters_nHCal(tree_reader, "HcalEndcapNClusterAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssocClusters_nHCal(tree_reader, "HcalEndcapNClusterAssociations.simID");

    TTreeReaderArray<unsigned int> nHCalClustershits_begin(tree_reader, "HcalEndcapNClusters.hits_begin");
    TTreeReaderArray<unsigned int> nHCalClustershits_end(tree_reader, "HcalEndcapNClusters.hits_end");
      
  // Define Histograms
    
    TH1D *nHCalRecHitsPosZ_all_nonAssoc = new TH1D("nHCalRecHitsPosZ_all_nonAssoc", "nHCalRecHits_all_nonAssoc Z; nHCalRecHits.position.z [mm]", 30,hz_min_nhcal,hz_max_nhcal);
    
    TH1D *nHCalRecHitsPosZ_all = new TH1D("nHCalRecHitsPosZ_all","nHCalRecHits_all Z; nHCalRecHits.position.z [mm]", 30,hz_min_nhcal,hz_max_nhcal);
    TH2D *nHCalRecHitsPosXY_all = new TH2D("nHCalRecHitsPosXY_all","nHCalRecHits_all XY; nHCalRecHits.position.x [mm]; nHCalRecHits.position.y [mm]", 50,hx_min_nhcal,hx_max_nhcal,50,hy_min_nhcal,hy_max_nhcal);
    TH3D *nHCalRecHitsPosXYZ_all = new TH3D("nHCalRecHitsPosXYZ_all","nHCalRecHits_all XYZ; nHCalRecHits.position.x [mm]; nHCalRecHits.position.y [mm]; nHCalRecHits.position.z [mm]", 30,hx_min_nhcal,hx_max_nhcal, 30,hy_min_nhcal,hy_max_nhcal, 30,hz_min_nhcal,hz_max_nhcal);
    TH1D *nHCalRecHitsE_all = new TH1D("nHCalRecHitsE_all","nHCalRecHits_all Energy; nHCalRecHits.energy [GeV]", 100,0.,12.);
    
    TH1D *nHCalClustersE_all = new TH1D("nHCalClustersE_all", "nHCalClusters_E_all; Cluster Energy [UNITS]; Counts", 100, 0, 6);
    
    TH1F* nHCalLayers_Arr[nlayers_nhcal]; // Energy By Layer

    for (int L = 0; L < nlayers_nhcal; L++) {
        nHCalLayers_Arr[L] = new TH1F(Form("nhCal_LayerE_%d", L),
                              Form("Hit Energies in Layer %d", L),
                              100, 0, 2);   // adjust binning as you want
    }
    
    
  //Creates Kaon Occurrence Histogram
  TH2D *kaonOccurrence = new TH2D("kaonOccurrence", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', Kaon HCal Acceptance (%);K_{1};K_{2}", 4, 0, 4, 4, 0, 4); // nbinsx, xlow, xup, nbinsy, ylow, yup
  
  //xB vs q2 histogram
  TH2F *xB_q2_hist = new TH2F("xB_q2_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay x_{Bj} vs Q^{2} Occurrence", 40, 1e-4, 0.02, 40, 1, 10); // nbinsx, xlow, xup, nbinsy, ylow, yup
  ConvertToLogBins2D(xB_q2_hist, 40, 1e-4, 0.02, 40, 1, 10);
    
  //K1 vs K2 eta histogram
  TH2F *K1_K2_eta_hist = new TH2F("K1_K2_eta_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay K_{1} vs K_{2} Pseudorapodity(#eta)", 60, -4.05, 4.2, 60, -4.05, 4.2); // nbinsx, xlow, xup, nbinsy, ylow, yup
    
  //Histogram displaying all eta of all kaons of run
  TH1F *all_eta = new TH1F("all_eta", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', K+K- eta Distribution; #eta", 100, -5, 5);
    
  //X-bjorken Histogram
  TH1F *xBjorken = new TH1F("xBjorken", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay x_{Bj} distribution; x_{Bj}", 75, 1e-5, 0.075);
  ConvertToLogBins(xBjorken, 75, 1e-5, 0.075);
    
    TH1F *reco_calcXb = new TH1F("reco_calcXb", "Calculated vs Reconstructed x_{Bj}", 75, 1e-5, 0.075);
    ConvertToLogBins(reco_calcXb, 75, 1e-5, 0.075);
    TH1F *reco_calcXb_Diff = new TH1F("reco_calcXb_Diff", "Reco Calculated vs ROOT x_{Bj} Difference (reco-root)", 75, -0.05, 0.05);
    TH1F *reco_calcXb_PercentError = new TH1F("reco_calcXb_PercentError", "ROOT vs Calculated x_{Bj} % Error (100(calc-reco)/reco)", 75, -300, 300);
    
    TH1F *reco_calcQ2 = new TH1F("reco_calcQ2", "Calculated vs Reconstructed Q^{2}", 50, 0.1, 10);
    ConvertToLogBins(reco_calcQ2, 50, 0.1, 10);
    TH1F *reco_calcQ2_Diff = new TH1F("reco_calcQ2_Diff", "Reco Calculated vs ROOT Q^{2} Difference (reco-root)", 75, -4, 4);
    TH1F *reco_calcQ2_PercentError = new TH1F("reco_calcQ2_PercentError", "ROOT vs Calculated Q^{2} % Error (100(calc-reco)/reco)", 75, -100, 100);
    
    TH1F *gen_calcXb = new TH1F("gen_calcXb", "Calculated vs Reconstructed x_{Bj}", 75, 1e-5, 0.075);
    ConvertToLogBins(gen_calcXb, 75, 1e-5, 0.075);
    TH1F *gen_calcXb_Diff = new TH1F("gen_calcXb_Diff", "Gen Calculated vs ROOT x_{Bj} Difference (gen-root)", 75, -0.05, 0.05);
    TH1F *gen_calcXb_PercentError = new TH1F("gen_calcXb_PercentError", "ROOT vs Calculated x_{Bj} % Error (100(calc-reco)/reco)", 75, -300, 300);
    
    TH1F *gen_calcQ2 = new TH1F("gen_calcQ2", "Calculated vs Reconstructed Q^{2}", 50, 0.1, 10);
    ConvertToLogBins(gen_calcQ2, 50, 0.1, 10);
    TH1F *gen_calcQ2_Diff = new TH1F("gen_calcQ2_Diff", "Gen Calculated vs ROOTd Q^{2} Difference (gen-root)", 75, -4, 4);
    TH1F *gen_calcQ2_PercentError = new TH1F("gen_calcQ2_PercentError", "ROOT vs Calculated Q^{2} % Error (100(calc-reco)/reco)", 75, -100, 100);
    
    TH1F *xBj_gen_reco_calc_diff = new TH1F("xBj_gen_reco_calc_diff", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay x_{Bj} Difference(Gen - Reco); x_{Bj}", 75, -0.05, 0.05);
    TH1F *q2_gen_reco_calc_diff = new TH1F("q2_gen_reco_calc_diff", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay x_{Bj} Difference(Gen - Reco); x_{Bj}", 75, -4, 4);
    
    TH1F *xBjorken_0 = new TH1F("xBjorken0", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay x_{Bj} distribution, 0 K+K- in nHCal; x_{Bj}", 75, 1e-4, 0.075);
    ConvertToLogBins(xBjorken_0, 75, 1e-4, 0.075);
    TH1F *xBjorken_1 = new TH1F("xBjorken1", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay x_{Bj} distribution, 1 K+K- in nHCal; x_{Bj}", 75, 1e-4, 0.075);
    ConvertToLogBins(xBjorken_1, 75, 1e-4, 0.075);
    TH1F *xBjorken_2 = new TH1F("xBjorken2", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay x_{Bj} distribution, 2 K+K- in nHCal; x_{Bj}", 75, 1e-4, 0.075);
    ConvertToLogBins(xBjorken_2, 75, 1e-4, 0.075);

  //q^2 Histogram
  TH1F *q2 = new TH1F("q2", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay Q^{2} distribution; Q^{2} (GeV^{2})", 50, 0.1,10);
    ConvertToLogBins(q2, 50, 0.1, 10);
    
    TH1F *q2_0 = new TH1F("q2_0", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay Q^{2} distribution, 0 K+K- in nHCal; Q^{2} (GeV^{2})", 50, 1,10);
      ConvertToLogBins(q2_0, 50, 1, 10);
    TH1F *q2_1 = new TH1F("q2_1", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay Q^{2} distribution, 1 K+K- in nHCal; Q^{2} (GeV^{2})", 50, 1,10);
      ConvertToLogBins(q2_1, 50, 1, 10);
    TH1F *q2_2 = new TH1F("q2_2", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay Q^{2} distribution, 2 K+K- in nHCal; Q^{2} (GeV^{2})", 50, 1,10);
      ConvertToLogBins(q2_2, 50, 1, 10);
    
     TH1F *newQ2 = new TH1F("newQ2", "Pythia: 18 x 275 Beam Effects, Q^{2} distribution; Q^{2}", 50, 1, 10);
    ConvertToLogBins(newQ2, 50, 1, 10);
    
    TH1F *newXb = new TH1F("newXb", "Pythia: 18 x 275 Beam Effects, x_{Bj} distribution; x_{Bj}", 75, 1e-5, 2);
    ConvertToLogBins(newXb, 75, 1e-5, 2);
    
  //xB_v_q2 Graph
  TGraph *xB_v_q2 = new TGraph();
    TH2F *reco_xB_q2_hist = new TH2F("reco_xB_q2_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay calculated x_{Bj} vs Q^{2}", 40, 1e-4, 0.02, 40, 0.4, 10);
    ConvertToLogBins2D(reco_xB_q2_hist, 40, 1e-4, 0.02, 40, 0.4, 10);
    TH2F *gen_xB_q2_hist = new TH2F("gen_xB_q2_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay calculated x_{Bj} vs Q^{2}", 40, 1e-4, 0.02, 40, 1, 10);
    ConvertToLogBins2D(gen_xB_q2_hist, 40, 1e-4, 0.02, 40, 1, 10);
    
    TH2F *gen_reco_xB_hist = new TH2F("gen_reco_xB_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay gen(x) v reco(y) calculated x_{Bj}", 40, 1e-4, 0.02, 40, 1e-4, 0.02);
    ConvertToLogBins2D(gen_reco_xB_hist, 40, 1e-4, 0.02, 40, 1e-4, 0.02);
    TH2F *gen_reco_q2_hist = new TH2F("gen_reco_q2_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay gen(x) v reco(y) calculated Q^{2}", 40, 1, 10, 40, 0.4, 10);
    ConvertToLogBins2D(gen_reco_q2_hist, 40, 1, 10, 40, 0.4, 10);
    
    TH2F *gen_root_xB_hist = new TH2F("gen_root_xB_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay gen calc(x) vs root(y) x_{Bj}", 40, 1e-4, 0.02, 40, 1e-4, 0.02);
    ConvertToLogBins2D(gen_root_xB_hist, 40, 1e-4, 0.02, 40, 1e-4, 0.02);
    TH2F *reco_root_xB_hist = new TH2F("reco_root_xB_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', #phi meson decay reco calc(x) vs root(y) x_{Bj}", 40, 1e-4, 0.02, 40, 1e-4, 0.02);
    ConvertToLogBins2D(reco_root_xB_hist, 40, 1e-4, 0.02, 40, 1e-4, 0.02);
    
    TGraph *eta_v_pT = new TGraph();
    TH2F *eta_v_pT_hist = new TH2F("eta_v_pT_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', K+K- pseudorapidity vs pT", 55, -4.5, 4.5, 55, 0.1, 4); // nbinsx, xlow, xup, nbinsy, ylow, yup
   ConvertToLogYBins2D(eta_v_pT_hist, 55, -4.5, 4.5, 55, 0.1, 4);
    
    TGraph *eta_v_xB = new TGraph();
    TH2F *eta_v_xB_hist = new TH2F("eta_v_xB_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', K+K- pseudorapidity vs x_{Bj}", 55, -4.5, 4.5, 55, 0.0001, 0.07); // nbinsx, xlow, xup, nbinsy, ylow, yup
    ConvertToLogYBins2D(eta_v_xB_hist, 55, -4.5, 4.5, 55, 0.0001, 0.07);

    TGraph *eta_v_q2 = new TGraph();
    TH2F *eta_v_q2_hist = new TH2F("eta_v_q2_hist", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', K+K- pseudorapidity vs Q^{2}", 55, -4.5, 4.5, 55, 1, 20); // nbinsx, xlow, xup, nbinsy, ylow, yup
    ConvertToLogYBins2D(eta_v_q2_hist, 55, -4.5, 4.5, 55, 1, 20);

    
  //xB vs percentage Histograms for 0, 1, 2 kaons in nHCal
  int xB_percent_nBins = 60;
  TH1F *xB_v_percent_0 = new TH1F("xB_v_percent_0", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', x_{Bj} vs percent ocurrence (%);x_{Bj}; geom. acc(%)", xB_percent_nBins, 0.00001, 1);
    ConvertToLogBins(xB_v_percent_0, xB_percent_nBins, 0.00001, 1);
  TH1F *xB_v_percent_1 = new TH1F("xB_v_percent_1", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', x_{Bj} vs percent ocurrence (%);x_{Bj}; geom. acc(%)", xB_percent_nBins, 0.00001, 1);
    ConvertToLogBins(xB_v_percent_1, xB_percent_nBins, 0.00001, 1);
  TH1F *xB_v_percent_2 = new TH1F("xB_v_percent_2", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', x_{Bj} vs percent ocurrence (%);x_{Bj}; geom. acc(%)", xB_percent_nBins, 0.00001, 1);
    ConvertToLogBins(xB_v_percent_2, xB_percent_nBins, 0.00001, 1);
  TH1F *xB_v_percent_all = new TH1F("xB_v_percent_all", "Sartre: e + Au#rightarrow e\'+#phi(KK)+Au\', x_{Bj} vs percent ocurrence (%);x_{Bj}; geom. acc(%)", xB_percent_nBins, 0.00001, 1);
    ConvertToLogBins(xB_v_percent_all, xB_percent_nBins, 0.00001, 1);
  TH1F *xB_v_percent_denom = new TH1F("xB_v_percent_denom", "Sartre: ", xB_percent_nBins, 0.00001, 1);
    ConvertToLogBins(xB_v_percent_denom, xB_percent_nBins, 0.00001, 1);
    
  //generatorStatus
  TH1D *generatorStatus = new TH1D("generatorStatus","Status of generated particles, all; generatorStatus",101,0,100);
  
  // eta (pseudorapidity)
  TH1D *partEta = new TH1D("partEta","Eta of thrown particles; #eta",120,-6.,6.);
  TH1D *recEta = new TH1D("recEta","Eta of reconstructed tracks that have matching thrown particle; #eta",120,-6.,6.);

  TH1D *electronEta = new TH1D("electronEta", "Eta of thrown e; #eta",120,-6.,6.);
  TH1D *electronRecEta = new TH1D("electronRecEta","Eta of reco e;#eta",120,-6.,6.);
  
  TH1D *protonEta = new TH1D("protonEta","Eta of thrown p;#eta",120,-6.,6.);
  TH1D *protonRecEta = new TH1D("protonRecEta","Eta of reco p;#eta",120,-6.,6.);
  
  TH1D *muonEta = new TH1D("muonEta","Eta of thrown #mu;#eta",120,-6.,6.);
  TH1D *muonRecEta = new TH1D("muonRecEta","Eta of reco #mu;#eta",120,-6.,6.);
  
  TH1D *muonPhi = new TH1D("muonPhi","Phi of thrown #mu;#phi",120,-6.,6.);
  TH1D *muonRecPhi = new TH1D("muonRecPhi","Phi of reco #mu;#phi",120,-6.,6.);
    
  TGraph *eta_v_phi_true = new TGraph();
  TH2F *eta_v_phi_true_hist = new TH2F("eta_v_phi_true_hist", "Eta vs Phi of thrown muon", 55, -6, 6, 55, -6, 6);
    
  TGraph *eta_v_phi_reco = new TGraph();

  TH1D *pionEta = new TH1D("pionEta","Eta of thrown #pi;#eta",120,-6.,6.);
  TH1D *pionRecEta = new TH1D("pionRecEta","Eta of reco #pi;#eta",120,-6.,6.);

  TH1D *pi0Eta = new TH1D("pi0Eta","Eta of thrown #pi;#eta",120,-6.,6.);
  
  TH1D *kaonEta = new TH1D("kaonEta","Eta of thrown K;#eta",120,-6.,6.);
  TH1D *kaonRecEta = new TH1D("kaonRecEta","Eta of reco K;#eta",120,-6.,6.);
  
  TH1D *rho0Eta = new TH1D("rho0Eta","Eta of thrown #rho^{0};#eta",120,-6.,6.);
  TH1D *pipmfromrho0Eta = new TH1D("pipmfromrho0Eta","generated #eta of pi^{#pm} from rho^{0} decay;#eta",120,-6.,6.);
  TH1D *pipmfromrho0RecEta = new TH1D("pipmfromrho0RecEta","reconstructed #eta of pi^{#pm} from rho^{0} decay;#eta",120,-6.,6.);
  
  // phi(1020)
  TH1D *phiEta = new TH1D("phiEta","Eta of thrown #phi(1020);#eta",120,-6.,6.);
  TH1D *kpmfromphiEta = new TH1D("kpmfromphiEta","generated #eta of K^{#pm} from #phi(1020) decay;#eta",120,-6.,6.);
  TH1D *kpmfromphiRecEta = new TH1D("kpmfromphiRecEta","reconstructed #eta of K^{#pm} from #phi(1020) decay;#eta",120,-6.,6.);
  
  // momentum
  TH1D *partMom = new TH1D("partMom","Mom of thrown particles; P [GeV]",150,0.,150.);
  TH1D *recP = new TH1D("recP","Momentum of reconstructed tracks; P [GeV]",150,0.,150.);
  
  // phi
  TH1D *partPhi = new TH1D("partPhi","Phi of thrown charged particles; #phi [rad]",150,-3.2,3.2);
  TH1D *recPhi = new TH1D("recPhi","Phi of reconstructed tracks; #phi [rad]",150,-3.2,3.2);
    

  // count events:
  int ievgen = 0;
  // count generated particles:
  int ngen_electrons = 0; //+-
  int ngen_protons = 0; //+-
  int ngen_muons = 0; //+-
  int ngen_pions = 0; //+-
  int ngen_pi0=0;
  int ngen_kaons = 0; //+-
  int ngen_rho0 = 0;
  int ngen_rhop = 0;
  int ngen_phi = 0;
  int ngen_omega = 0;
  int ngen_jpsi = 0;
  int ngen_upsilon = 0;
  int ngen_d0 = 0;
  int ngen_b0 = 0;
  // count reconstructed particles (+-):
  int nrec_electrons = 0;
  int nrec_protons = 0;
  int nrec_muons = 0;
  int nrec_pions = 0;
  int nrec_kaons = 0;
  // count number of decays:
  int ndecay_pi0_gg = 0;
  int ndecay_rho0_pp = 0;
  int ndecay_rho0_mumu = 0;
  int ndecay_rho0_ee = 0;
  int ndecay_phi_kk = 0;
  // count number of decay particles (reco level) in nHCal acceptance:
  int ndecay_rho0_pionpm_nHCal = 0;
  int ndecay_phi_kaonpm_nHCal = 0;
  int ndecay_phi_kaon_pair_nHCal = 0;
  int ndecay_phi_kaon_single_nHCal = 0;
  int ndecay_phi_kaon_none_nHCal = 0;
  // tag decays on generated particle level:
  int is_rho0decay_pp = 0; // 0 or 1 for a given generated particle
  int is_rho0decay_mumu = 0;
  int is_rho0decay_ee = 0;
  int is_phidecay_kk = 0;
    int numTimes = 0;
  int nreco_phi = 0; //total # of phi decay with both kaons reconstructed

  cout << "+ Ready to run over events... \n";
    int inB = 0;

  while(tree_reader.Next()) { // Loop over events

      //newQ2->Fill(partQ2[0]); //PYTHIA STUFF
      //newXb->Fill(partXb[0]);

    ievgen++;
      
    int kaons_in_nHCal = 0; //set the number of kaons within tolerance to 0 and reset for each particle
    float k1_eta;
    float k2_eta;
    //cout << "+ Entering event #: " << ievgen << " \n";
    
    //cout << "Event #: " << ievgen << ", " << partGenStat.GetSize() << " gen particles, " << parents_index.GetSize() << " parent particles, " << daughters_index.GetSize() << " daughter particles \n";   // parent_index and daughter_index must be of the same length since they are in the same tree (is that what pushback does?)

    for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over generated particles
      {
    //cout << "++ Entering generated particle #: " << i << " \n";
    
    
    int pdg = TMath::Abs(partPdg[i]);
    TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);
    float trueEta = trueMom.PseudoRapidity();
    float truePhi = trueMom.Phi();
    float trueTheta = trueMom.Theta();
    
          
      for (unsigned int h = 0; h < nHCalRecHitsPosZ.GetSize(); h++) {
          nHCalRecHitsPosZ_all_nonAssoc->Fill(nHCalRecHitsPosZ[h]);
      }
    generatorStatus->Fill(partGenStat[i]);
    // reset particle decays for each new generated particle:
    is_rho0decay_pp = 0;
    is_rho0decay_mumu = 0;
    is_rho0decay_ee = 0;
    is_phidecay_kk = 0;
          
    int reco_kk = 0;
    
    int i_parents_begin = parents_begin[i];
        int i_parents_end = parents_end[i];
    int i_parents = parents_end[i] - parents_begin[i];

    int i_daughters_begin = daughters_begin[i];
        int i_daughters_end = daughters_end[i] - 1;
    int i_daughters = daughters_end[i] - daughters_begin[i];

    //Consider only selected generated particles:
    if( (partGenStat[i] == 1) || (partGenStat[i] == 2 || (partGenStat[i] == 4)) /*|| (partGenStat[i] = 4 */) // Select only stable or decay particles
    {
      //cout << "Ev#: " << ievgen << ", P-index: " << i <<", PDG: " << partPdg[i] << ", GenStatus:" << partGenStat[i] << ", i_parents: "<< i_parents<<", i_daughters: " << i_daughters << ", pb: " << parents_index[i_parents_begin] << ", pe: " << parents_index[i_parents_end] <<  ", db: " << daughters_index[i_daughters_begin] << ", de: " << daughters_index[i_daughters_end] << " \n";
      // rho0 decays
      if( partPdg[i] == 113 )
        {
          //cout << "Event " << ievgen << " with gen rho0 #: " << partPdg[i] << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";

          // count the 2-body rho0 decays:
          if( i_daughters == 2 )
        {
          if( (partPdg[daughters_index[i_daughters_begin]] == 211 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -211 ) || (partPdg[daughters_index[i_daughters_begin]] == -211 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 211) )
            {
              //cout << "-> Event " << ievgen << " found rho0 decayed into pi+ pi-: "  << partPdg[daughters_index[i_daughters_begin]] << " and " << partPdg[daughters_index[i_daughters_begin]+1] << "  \n";
              
              // count it:
              ndecay_rho0_pp++;
              // tag it:
              is_rho0decay_pp = 1;
              
              TVector3 trueMom_rho0_pi1(partMomX[daughters_index[i_daughters_begin]],partMomY[daughters_index[i_daughters_begin]],partMomZ[daughters_index[i_daughters_begin]]);
              TVector3 trueMom_rho0_pi2(partMomX[daughters_index[i_daughters_begin]+1],partMomY[daughters_index[i_daughters_begin]+1],partMomZ[daughters_index[i_daughters_begin]+1]);
              TVector3 trueMom_rho0_pi12=trueMom_rho0_pi1 + trueMom_rho0_pi2;
              
              float trueEta_rho0_pi1 = trueMom_rho0_pi1.PseudoRapidity();
              float trueEta_rho0_pi2 = trueMom_rho0_pi2.PseudoRapidity();
              pipmfromrho0Eta->Fill(trueEta_rho0_pi1);
              pipmfromrho0Eta->Fill(trueEta_rho0_pi2);

              //          cout << "--> Event " << ievgen << " rho0 decay to 2pi: generated rho0 eta: " << trueEta << ", pi1: " << trueEta_rho0_pi1 << ", pi2: " << trueEta_rho0_pi2 << "  \n";
              //  cout << "            trueMomrho0 X: " << trueMom.X() << ", trueMomrho0 Y: " << trueMom.Y() <<", trueMomrho0 Z: " << trueMom.Z() << "  \n";
              // cout << "            trueMom_rho0_pi12 X: " << trueMom_rho0_pi12.X() << ", trueMom_rho0_pi12 Y: " << trueMompi12.Y() <<", trueMom_rho0_pi12 Z: " << trueMom_rho0_pi12.Z() << "  \n";
            
            } // end of rho0 to pi+pi- decays
          else if( (partPdg[daughters_index[i_daughters_begin]] == 13 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -13 ) || (partPdg[daughters_index[i_daughters_begin]] == -13 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 13) )
            {
              ndecay_rho0_mumu++;
            } // end of rho0 to mumu decays
          else if( (partPdg[daughters_index[i_daughters_begin]] == 11 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -11 ) || (partPdg[daughters_index[i_daughters_begin]] == -11 ) && ( partPdg[daughters_index[i_daughters_begin]] == 11) )
            {
              ndecay_rho0_ee++;
            } // end of ee decays
        } // end of 2-body decays
        } // end of rho0 decays
      // pi0 decays
      else if( partPdg[i] == 111 )
        {
          //cout << "Event with gen pi0 #: " << ievgen << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";
          if( (partPdg[daughters_index[i_daughters_begin]] == 22) && (partPdg[daughters_index[i_daughters_begin]+1] == 22) )
        {
          ndecay_pi0_gg++;
        } // end of gg decays
        }// end of pi0 decay

      // phi(1020) decays:
            if( partPdg[i] == 333 )
        {
          //cout << "Event " << ievgen << " with gen phi(1020): " << partPdg[i] << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";

          // count the 2-body phi(1020) decays:
          if( i_daughters == 2 )
        {
          if( (partPdg[daughters_index[i_daughters_begin]] == 321 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -321 ) || (partPdg[daughters_index[i_daughters_begin]] == -321 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 321) )
            {
              //cout << "-> Event " << ievgen << " found phi(1020) decayed into K+ K-: "  << partPdg[daughters_index[i_daughters_begin]] << " and " << partPdg[daughters_index[i_daughters_begin]+1] << "  \n";
              
              // count it:
              ndecay_phi_kk++;
              // tag it:
              is_phidecay_kk = 1;
              decays.push_back(phiDecay());

              TVector3 trueMom_phi_k1(partMomX[daughters_index[i_daughters_begin]],partMomY[daughters_index[i_daughters_begin]],partMomZ[daughters_index[i_daughters_begin]]);
              TVector3 trueMom_phi_k2(partMomX[daughters_index[i_daughters_begin]+1],partMomY[daughters_index[i_daughters_begin]+1],partMomZ[daughters_index[i_daughters_begin]+1]);
              TVector3 trueMom_phi_k12=trueMom_phi_k1 + trueMom_phi_k2;
              
              float trueEta_phi_k1 = trueMom_phi_k1.PseudoRapidity();
              float trueEta_phi_k2 = trueMom_phi_k2.PseudoRapidity();
              kpmfromphiEta->Fill(trueEta_phi_k1);
              kpmfromphiEta->Fill(trueEta_phi_k2);
                decays[decays.size() - 1].genEta1 = trueEta_phi_k1;
                decays[decays.size() - 1].genEta2 = trueEta_phi_k2;

              //cout << "--> Event " << ievgen << " phi(1020) decay to 2pi: generated phi eta: " << trueEta << ", pi1: " << trueEta_phi_k1 << ", pi2: " << trueEta_phi_k2 << "  \n";
              //cout << "            trueMomphi X: " << trueMom.X() << ", trueMomphi Y: " << trueMom.Y() <<", trueMomphi Z: " << trueMom.Z() << "  \n";
              //cout << "            trueMom_phi_k12 X: " << trueMom_phi_k12.X() << ", trueMom_phi_k12 Y: " << trueMom_phi_k12.Y() <<", trueMom_phi_k12 Z: " << trueMom_phi_k12.Z() << "  \n";
            
            } // end of phi to pi+pi- decays
        } // end of 2-body decays
        } // end of phi meson decays
      // count any generated particles:
      if( pdg == 11){
        ngen_electrons++;
        electronEta->Fill(trueEta);
          if (i_parents == 1 && partGenStat[i] == 1) {
              if (!decays.empty()) {
                  decays[decays.size()-1].gen_elec_mom_x_f = partMomX[i];
                  decays[decays.size()-1].gen_elec_mom_y_f = partMomY[i];
                  decays[decays.size()-1].gen_elec_mom_z_f = partMomZ[i];
                  //decays[decays.size()-1].reco_elec_energy_f = partEnergy[i];
              }
              /*cout << "the scattered electron has Q2 " << partQ2[0] << " and Xb " << partXb[0] << " with energy " << CPartEnergy << " and momentum x " << trackMomX[recoAssoc[j]] << " y " << trackMomY[recoAssoc[j]] << " z " << trackMomZ[recoAssoc[j]] << "\n";*/
              }
      }// electrons
      else if( pdg == 13){
        ngen_muons++;
        muonEta->Fill(trueEta);
        muonPhi->Fill(truePhi);
        eta_v_phi_true->AddPoint(trueEta, truePhi);
          eta_v_phi_true_hist->Fill(trueEta,truePhi);
      }// muons
      else if( pdg == 211){
        ngen_pions++;
        pionEta->Fill(trueEta);
      }//pions_pm
      else if( pdg == 111){
        ngen_pi0++;
        pi0Eta->Fill(trueEta);
      }//pions_pm
      else if( pdg == 321 ){
        ngen_kaons++;
        kaonEta->Fill(trueEta);
      } // kaons_pm
      else if( pdg == 113){
        ngen_rho0++;
        rho0Eta->Fill(trueEta);
      } // rho(770)
      else if( pdg == 443){
        ngen_jpsi++;
      } // J/Psi(1S)
      else if( pdg == 2212){
        ngen_protons++;
        protonEta->Fill(trueEta);
      }// protons
      else if( pdg == 213){
        ngen_rhop++;
      }// rhop
      else if( pdg == 333){
        ngen_phi++;
        phiEta->Fill(trueEta);
      }// phi(1020)
      else if( pdg == 223){
        ngen_omega++;
      }// omega(982)
      else if( pdg == 553){
        ngen_upsilon++;
      }// Upsilon(1S)
      else if( pdg == 421){
        ngen_d0++;
      }// D0
      else if( pdg == 511){
        ngen_b0++;
      }// B0
        
        
      //Fill all true eta:
      partEta->Fill(trueEta);
        
      // Fill all true momentum:
      partMom->Fill(trueMom.Mag());
        
      // Fill all true phi:
      partPhi->Fill(truePhi);
      //hit level
        //cout << partGenStat.GetSize() << "\n";
        for (unsigned int k = 0; k < simuAssocClusters_nHCal.GetSize(); k++) {
            //nHCalRecHitsPosZ_all_nonAssoc->Fill(nHCalRecHitsPosZ[k]); // NON ASSOCIATED
            if (simuAssocClusters_nHCal[k] == i) {
                numTimes++;
                //cout << "MCParticle: " << /*i*/numTimes << ", PDG: " << partPdg[i] << ", matching cluster ID in the nHCal: " << k << ", cluster energy: " << nHCalClustersE[recoAssocClusters_nHCal[k]] << ", cluster position X: " << nHCalClustersPosX[recoAssocClusters_nHCal[k]] <<  ", cluster position Y: " << nHCalClustersPosY[recoAssocClusters_nHCal[k]] <<  ", cluster position Z: " << nHCalClustersPosZ[recoAssocClusters_nHCal[k]] << ", hits begin: " << nHCalClustershits_begin[recoAssocClusters_nHCal[k]] << ", hits end: " << nHCalClustershits_end[recoAssocClusters_nHCal[k]] << " \n";
                nHCalClustersE_all->Fill(nHCalClustersE[k]);
                
                for(unsigned int h=nHCalClustershits_begin[recoAssocClusters_nHCal[k]]; h<nHCalClustershits_end[recoAssocClusters_nHCal[k]]; h++){
                    //cout << "asso hit: " << h << " with energy: " << nHCalRecHitsE[h] << " \n";
                    // Fill hit level info:
                    nHCalRecHitsPosZ_all->Fill(nHCalRecHitsPosZ[h]);
                    //cout << "pos Z:" << nHCalRecHitsPosZ[h] << " and layer " << layer_zs[nHCalRecHitsPosZ[h]] << " with energy " << nHCalRecHitsE[h] << "\n";
                    
                    nHCalLayers_Arr[layer_zs[nHCalRecHitsPosZ[h]] - 1]->Fill(nHCalRecHitsE[h]);
                    ///cout << "FILLED AT" << nHCalRecHitsPosZ[h] << "\n";
                    nHCalRecHitsPosXY_all->Fill(nHCalRecHitsPosX[h], nHCalRecHitsPosY[h]);
                    nHCalRecHitsPosXYZ_all->Fill(nHCalRecHitsPosX[h], nHCalRecHitsPosY[h], nHCalRecHitsPosZ[h]);
                    nHCalRecHitsE_all->Fill(nHCalRecHitsE[h]);
                }
            }
        }
        
      // Loop over associations to find matching ReconstructedChargedParticle
      for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
        {
          //cout << "*** Event " << ievgen << ", generated particle " << i << ", simID " << j << " \n";
            
          if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
        {
          //cout << "***** Event " << ievgen << ", found association index: " << simuAssoc[j] << " \n";
          
          TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle
          
          float CPartMom = recMom.Mag();
          float CPartEta = recMom.PseudoRapidity();
          float CPartPhi = recMom.Phi();
          float CPartTheta = recMom.Theta();
          float CPartEnergy = trackPartEnergy[recoAssoc[j]];
            
          recEta->Fill(CPartEta);
          recPhi->Fill(CPartPhi);
          recP->Fill(recMom.Mag());
          //cout << "Particle is pdg: " << pdg << " .\n";
          
          if(pdg == 11){
            nrec_electrons++;
            electronRecEta->Fill(CPartEta);
              if (i_parents == 1 && partGenStat[i] == 1) {
                  if (!decays.empty()) {
                      decays[decays.size()-1].reco_electron = true;
                      decays[decays.size()-1].reco_elec_mom_x_f = trackMomX[recoAssoc[j]];
                      decays[decays.size()-1].reco_elec_mom_y_f = trackMomY[recoAssoc[j]];
                      decays[decays.size()-1].reco_elec_mom_z_f = trackMomZ[recoAssoc[j]];
                      decays[decays.size()-1].reco_elec_energy_f = CPartEnergy;
                  }
                  /*cout << "the scattered electron has Q2 " << partQ2[0] << " and Xb " << partXb[0] << " with energy " << CPartEnergy << " and momentum x " << trackMomX[recoAssoc[j]] << " y " << trackMomY[recoAssoc[j]] << " z " << trackMomZ[recoAssoc[j]] << "\n";*/
                  }
          }// electrons
          else if( pdg == 13){
            nrec_muons++;
            muonRecEta->Fill(CPartEta);
            muonRecPhi->Fill(CPartPhi);
            eta_v_phi_reco->AddPoint(CPartEta, CPartPhi);
          }// muons
          else if( pdg == 211){
            nrec_pions++;
            pionRecEta->Fill(CPartEta);
          }//pions
          else if( pdg == 321){
            nrec_kaons++;
            kaonRecEta->Fill(CPartEta);
          }//pions
          else if( pdg == 2212){
            nrec_protons++;
            protonRecEta->Fill(CPartEta);
          }// protons
        } // end of matched association gen to rec
          
          // Match the decay particles to their recos:
          // rho0 to pi+ pi-
          if( is_rho0decay_pp )
        {
          if( simuAssoc[j] == daughters_index[i_daughters_begin] ) // get the reco decay pi1 of the gen rho0, by accessing the correct MCParticle index
            {
              TVector3 recMom_rho0_pi1(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
              float recEta_rho0_pi1 = recMom_rho0_pi1.PseudoRapidity();
              float recPhi_rho0_pi1 = recMom_rho0_pi1.Phi();
              pipmfromrho0RecEta->Fill(recEta_rho0_pi1);

              // count the decay pions (reco level) that are within the nHCal acceptance, here pion1:
              if(in_Cal_Tolerance(nHCal_name, recEta_rho0_pi1))
            {
              ndecay_rho0_pionpm_nHCal++;
            }
              //cout << "---> Event " << ievgen << " rho0 decay, reco index rho0: " << j << " \n";
              //cout << "          reco daughter-1 eta: " << recEta_rho0_pi1 << ", reco index daughter-1: " << daughters_index[i_daughters_begin] << " \n";
            }// end of rho0 decay pi1
          else if( simuAssoc[j] == daughters_index[i_daughters_begin]+1 ) // get the reco decay pi2 of the gen rho0
            {
              TVector3 recMom_rho0_pi2(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
              float recEta_rho0_pi2 = recMom_rho0_pi2.PseudoRapidity();
              float recPhi_rho0_pi2 = recMom_rho0_pi2.Phi();
              pipmfromrho0RecEta->Fill(recEta_rho0_pi2);

              // count the decay pions (reco level) that are within the nHCal acceptance, here pion2:
              if(in_Cal_Tolerance(nHCal_name, recEta_rho0_pi2))
            {
              ndecay_rho0_pionpm_nHCal++;
            }
              
              //cout << "          reco daughter-2 eta: " << recEta_rho0_pi2  << ", reco index daughter-2: " << daughters_index[i_daughters_begin]+1 << " \n\n";
            }// end of rho0 decay pi2
        } // end of rho0 decay into pp
          // phi to K+ K-
          else if( is_phidecay_kk )
        {
            //Create a phiDecay object for this decay and add it to the array
            decays[decays.size()-1].x_b = partXb[0];
            decays[decays.size()-1].q2 = partQ2[0];
            
          if( simuAssoc[j] == daughters_index[i_daughters_begin] ) // get the reco decay pi1 of the gen rho0, by accessing the correct MCParticle index
            {
                reco_kk++;
                decays[decays.size() - 1].k1_reco = true;
                TVector3 recMom_phi_k1(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                float recEta_phi_k1 = recMom_phi_k1.PseudoRapidity();
                float recPhi_phi_k1 = recMom_phi_k1.Phi();
              kpmfromphiRecEta->Fill(recEta_phi_k1);
                decays[decays.size()-1].eta1 = recEta_phi_k1;
                decays[decays.size()-1].k1_pT = sqrt((trackMomX[recoAssoc[j]] * trackMomX[recoAssoc[j]])
                                                     + (trackMomY[recoAssoc[j]] * trackMomY[recoAssoc[j]]));
                k1_eta = recEta_phi_k1;
                if (in_Cal_Tolerance(cals[1], k1_eta)) {
                    inB++;
                }
                all_eta->Fill(k1_eta);
                
                //Check that the kaon is within a given HCal, in this case nHCal
                if(in_Cal_Tolerance(nHCal_name, recEta_phi_k1))
              {
                ndecay_phi_kaonpm_nHCal++;
                  kaons_in_nHCal++;
              }
              //cout << "---> Event " << ievgen << " phi(1020) decay, reco index phi(1020): " << j << " \n";
              //cout << "          reco daughter-1 eta: " << recEta_phi_k1 << ", reco index daughter-1: " << daughters_index[i_daughters_begin] << " \n";
            }// end of phi(1020) decay K1
          else if( simuAssoc[j] == daughters_index[i_daughters_begin]+1 ) // get the reco decay pi2 of the gen rho0
            {
                reco_kk++;
                decays[decays.size() - 1].k2_reco = true;
                TVector3 recMom_phi_k2(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                float recEta_phi_k2 = recMom_phi_k2.PseudoRapidity();
                float recPhi_phi_k2 = recMom_phi_k2.Phi();
              kpmfromphiRecEta->Fill(recEta_phi_k2);
                decays[decays.size()-1].eta2 = recEta_phi_k2;
                decays[decays.size()-1].k2_pT = sqrt((trackMomX[recoAssoc[j]] * trackMomX[recoAssoc[j]]) + (trackMomY[recoAssoc[j]] * trackMomY[recoAssoc[j]]));
                k2_eta = recEta_phi_k2;
                if (in_Cal_Tolerance(cals[1], k2_eta)) {
                    inB++;
                }
                all_eta->Fill(k2_eta);
                
                //Check that the kaon is within a given HCal, in this case nHCal
                if(in_Cal_Tolerance(nHCal_name, recEta_phi_k2))
              {
                ndecay_phi_kaonpm_nHCal++;
                  kaons_in_nHCal++;
              }
              //cout << "          reco daughter-2 eta: " << recEta_phi_k2  << ", reco index daughter-2: " << daughters_index[i_daughters_begin]+1 << " \n\n";
            }// end of phi(1020) decay K2
        
            /*if (recEta_phi_k1 == recEta_phi_k2) {
                cout << "true\n";
            } else {
                cout << "false\n";
            }*/
            decays[decays.size()-1].reco_parts = reco_kk;
            if (reco_kk == 2) {
                nreco_phi++;
            }
            
            //Check if only one of the two kaons was within the nHCal acceptance
        } // end of phi(1020) decay into KK
          
        }// End loop over associations
    } // End stable or decay particles condition
      } // End loop over thrown particles, within that event

    //for(unsigned int k=0; k<trackMomX.GetSize(); k++){ // Loop over all reconstructed tracks, thrown or not

    //  TVector3 recMom(trackMomX[k], trackMomY[k], trackMomZ[k]);
    
    //} // End loop over all reconstructed tracks, within that event

    // now go to next event
    
  } // End loop over events
    int numElec = 0;
    for (phiDecay decay : decays) {
        if (decay.reco_parts == 2 && decay.reco_electron) {
            kaons_in_Cal(decay.eta1, decay.eta2);
            fill_Cal_Arr(decay.genEta1, decay.genEta2);
            decay.calcRecoXbQ2();
            decay.calcGenXbQ2();
            
            //Fill Graphs/Histograms
            xBjorken->Fill(decay.x_b);
            q2->Fill(decay.q2);
            reco_calcXb->Fill(decay.reco_calcXb);
            reco_calcXb_PercentError->Fill(100*(decay.reco_calcXb-decay.x_b)/(decay.x_b));
            reco_calcQ2->Fill(decay.reco_calcQ2);
            reco_calcQ2_PercentError->Fill(100*(decay.reco_calcQ2-decay.q2)/(decay.q2));
            gen_calcXb->Fill(decay.gen_calcXb);
            gen_calcXb_PercentError->Fill(100*(decay.gen_calcXb-decay.x_b)/(decay.x_b));
            gen_calcQ2->Fill(decay.gen_calcQ2);
            gen_calcQ2_PercentError->Fill(100*(decay.gen_calcQ2-decay.q2)/(decay.q2));
            xBj_gen_reco_calc_diff->Fill(decay.gen_calcXb-decay.reco_calcXb);
            q2_gen_reco_calc_diff->Fill(decay.gen_calcQ2-decay.reco_calcQ2);
            reco_calcXb_Diff->Fill(decay.reco_calcXb - decay.x_b);
            reco_calcQ2_Diff->Fill(decay.reco_calcQ2 - decay.q2);
            gen_calcXb_Diff->Fill(decay.gen_calcXb - decay.x_b);
            gen_calcQ2_Diff->Fill(decay.gen_calcQ2 - decay.q2);
            
            decay.set_num_in_nHCal();
            
            if (decay.num_in_nHCal == 0) {
                xBjorken_0->Fill(decay.x_b);
                q2_0->Fill(decay.q2);
            } else if (decay.num_in_nHCal == 1) {
                xBjorken_1->Fill(decay.x_b);
                q2_1->Fill(decay.q2);
            } else if (decay.num_in_nHCal == 2) {
                xBjorken_2->Fill(decay.x_b);
                q2_2->Fill(decay.q2);
            }
            
            xB_v_q2->AddPoint(decay.x_b, decay.q2);
            
            eta_v_pT->AddPoint(decay.eta1, decay.k1_pT);
            eta_v_xB->AddPoint(decay.genEta1, decay.x_b);
            eta_v_q2->AddPoint(decay.eta1, decay.q2);
            
            eta_v_pT->AddPoint(decay.eta2, decay.k2_pT);
            eta_v_xB->AddPoint(decay.genEta2, decay.x_b);
            eta_v_q2->AddPoint(decay.eta2, decay.q2);
            
            eta_v_pT_hist->Fill(decay.eta1, decay.k1_pT);
            eta_v_pT_hist->Fill(decay.eta2, decay.k2_pT);

            eta_v_xB_hist->Fill(decay.eta1, decay.x_b);
            eta_v_xB_hist->Fill(decay.eta2, decay.x_b);
            
            eta_v_q2_hist->Fill(decay.eta1, decay.q2);
            eta_v_q2_hist->Fill(decay.eta2, decay.q2);
            
            xB_q2_hist->Fill(decay.x_b, decay.q2);
            K1_K2_eta_hist->Fill(decay.eta1, decay.eta2);
            reco_xB_q2_hist->Fill(decay.reco_calcXb, decay.reco_calcQ2);
            gen_xB_q2_hist->Fill(decay.gen_calcXb, decay.gen_calcQ2);
            gen_reco_xB_hist->Fill(decay.gen_calcXb, decay.reco_calcXb);
            gen_reco_q2_hist->Fill(decay.gen_calcQ2, decay.reco_calcQ2);
            gen_root_xB_hist->Fill(decay.gen_calcXb, decay.x_b);
            reco_root_xB_hist->Fill(decay.reco_calcXb, decay.x_b);
        }
        /*if (decay.reco_parts == 2 && decay.reco_electron == false) {
            cout << "Event Number " << numElec << " has no reco elec\n";
        }*/
        numElec++;
    }
    
    eta_v_phi_true->Write("eta_v_phi_true");
    eta_v_phi_true_hist->Write("eta_v_phi_true_hist");
    
    eta_v_phi_reco->Write("eta_v_phi_reco");

  //
  
  cout << "Number of generated events: " << ievgen << " \n\n";
  cout << "Number of generated electrons +-: " << ngen_electrons << " \n";
  cout << "Number of generated protons +-: " << ngen_protons << " \n";
  cout << "Number of generated muons +-: " << ngen_muons << " \n";
  cout << "Number of generated pions +-: " << ngen_pions << " \n";
  cout << "Number of generated pi0: " << ngen_pi0 << ", of which decay into 2 gamma: " << ndecay_pi0_gg <<  " \n";
  cout << "Number of generated kaons +-: " << ngen_kaons << " \n";
  cout << "Number of generated rho0: " << ngen_rho0 << ", of which decay into pi+ pi-: " << ndecay_rho0_pp << ", into mu+ mu-: " << ndecay_rho0_mumu << ", into e+ e-: " << ndecay_rho0_ee << " \n";
  cout << "Number of generated omega: " << ngen_omega << " \n";
  cout << "Number of generated J/Psi: " << ngen_jpsi << " \n";
  cout << "Number of generated Upsilon: " << ngen_upsilon << " \n";
  cout << "Number of generated D0: " << ngen_d0 << " \n";
  cout << "Number of generated B0: " << ngen_b0 << " \n\n";
  cout << "Number of reconstructed electrons +-: " << nrec_electrons << " \n";
  cout << "Number of reconstructed protons +-: " << nrec_protons << " \n";
  cout << "Number of reconstructed muons +-: " << nrec_muons << " \n";
  cout << "Number of reconstructed pions +-: " << nrec_pions << " \n";
  cout << "Number of reconstructed kaons +-: " << nrec_kaons << " \n";
    
  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file

  cout << "Output histograms written in: " << outfile << " \n";
  
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
}

