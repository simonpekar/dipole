#include "constants.h"

/* random generator in )0,1) */
TRandom * sim = new TRandom;

void BinLogX(TH1 * h)
{
    
    TAxis *axis = h->GetXaxis();
    int bins = axis->GetNbins();
    
    Axis_t from = axis->GetXmin();
    Axis_t to = axis->GetXmax();
    Axis_t width = (to - from) / bins;
    Axis_t *new_bins = new Axis_t[bins + 1];
    
    for (int i = 0; i <= bins; i++) {
        new_bins[i] = TMath::Power(10, from + i * width);
    }
    
    axis->Set(bins, new_bins);
    delete new_bins;
    
}

void analyse() {
    
    UInt_t ancestor;
    Double_t y,size;
    
    UInt_t random_leaves[k_leaves];
    UInt_t r;
    UInt_t num;
    
    
    TF1 * fit_Bessel = new TF1("fitBessel", "[0]*TMath::BesselI0(2.*sqrt(-2.*[1]*log(x)))");
    TF1 * P_n = new TF1("P_n", "1./x * [3]*[2]/([0]*[0]) * exp(-log(x)*log(x)/(4*[1]))");
    TF1 * n_bar = new TF1("n_bar", "1./x * exp([0]*4.*log(2.)) * exp(-log(x*x)*log(x*x)/([0]*56.*[1])) * 1./sqrt([0]*14.*TMath::Pi()*[1])");
    
    n_bar->FixParameter(0, y_max);
    n_bar->FixParameter(1, zeta_3);
    
    fit_Bessel->SetParameter(0,N);
    fit_Bessel->FixParameter(1,y_max);
    
    P_n->FixParameter(0, r_s);
    P_n->FixParameter(1, y_max);
    P_n->FixParameter(2, N);
    
    //TH1 * h_size = new TH1D("h_size", "Size of leaves", 10000, 0, 2);
    TH1 * h_rap = new TH1D("h_rap", TString::Format("Rapidity of ancestors (N = %d, y_max = %.12g, delta = %.12g, k = %d, m = %d)", N, y_max, delta, k_leaves, m_factor), 50, 0, y_max);
    TH1 * h_anc = new TH1D("h_anc", TString::Format("Size of ancestors (N = %d, y_max = %.12g, delta = %.12g, k = %d, m = %d)", N, y_max, delta, k_leaves, m_factor), 50, -2, 2);
    TH1 * h_fluct = new TH1I("h_fluct", TString::Format("Fluctuations of multiplicity (N = %d, y_max = %.12g, delta = %.12g)", N, y_max, delta), 20, 0, 10*n_bar->Eval(delta));
    
    BinLogX(h_anc);
    
    TCanvas *c1 = new TCanvas("Modèle des dipoles");
    
    TFile f("tree.root","read");
    
    for (UInt_t n = 0; n < N; ++n) {
        
        printf("%.12g %%\n",100.*n/N);
        
        TTree * t;
        TTree * l;
        
        f.GetObject(TString::Format("t%d",n),t);
        f.GetObject(TString::Format("l%d",n),l);
        
        t->SetBranchAddress("size",&size);
        t->SetBranchAddress("rapidity",&y);
        t->SetBranchAddress("ancestor",&ancestor);
        
        l->SetBranchAddress("size",&size);
        l->SetBranchAddress("rapidity",&y);
        l->SetBranchAddress("ancestor",&ancestor);
        
        if (compute_ancestors) {
            
            if (l->GetEntries() >= k_leaves && l->GetEntries() >= (int)m_factor*n_bar->Eval(delta)) {
                
                for (UInt_t k = 0; k < k_leaves; k++) {
                    
                    random_leaves[k] = l->GetEntries();
                    
                    do {
                        
                        r = (int)l->GetEntries()*(1-sim->Rndm());
                        
                    } while (std::find(random_leaves, random_leaves+k, r) != random_leaves+k);
                    
                    random_leaves[k] = r;
                    
                }
                
                for (UInt_t k = 0; k < k_leaves - 1; k++) {
                    
                    while (random_leaves[k] != random_leaves[k+1]) {
                        
                        if (random_leaves[k] < random_leaves[k+1]) {
                            
                            t->GetEvent(random_leaves[k+1]);
                            random_leaves[k+1] = ancestor;
                            
                        } else {
                            
                            t->GetEvent(random_leaves[k]);
                            random_leaves[k] = ancestor;
                            
                        }
                        
                    }
                    
                }
                
                h_anc->Fill(size);
                h_rap->Fill(y);
                
            }
            
        }
        
        if (compute_fluctuations) {
            
            num = l->GetEntries(TString::Format("size >= %.12g", r_s));
            h_fluct->Fill(num);
            
        }
        
    }
    
    
    c1->Divide(1,2);
    
    /*
    //c1->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    h_size = h_size->GetCumulative(kFALSE);
    h_size->Draw("E1");
    gStyle->SetOptFit(1011);
    h_size->Fit("fitBessel", "", "", delta, 0.1);
    */
    
    c1->cd(1);
    h_anc->GetXaxis()->SetTitle("Radius r");
    h_anc->GetYaxis()->SetTitle("Number of ancestors of radius r");
    gPad->SetLogx();
    //gPad->SetLogy();
    h_anc->Draw();
    //ancestors->Draw("size",NULL);
    
    
    c1->cd(2);
    h_rap->GetXaxis()->SetTitle("Rapidity y");
    h_rap->GetYaxis()->SetTitle("Number of ancestors of rapidity y");
    //gPad->SetLogx();
    //gPad->SetLogy();
    h_rap->Draw();
    //ancestors->Draw("size",NULL);
    
    
    /*
    //c1->cd(4);
    gPad->SetLogy();
    h_fluct->Draw("E1");
    h_fluct->Fit("P_n", "", "", n_bar->Eval(delta), 10*n_bar->Eval(delta));
    */
}
