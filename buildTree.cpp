/* constants */
const Double_t zeta_3 = 1.2020569031595942854;

/* general parameters */
const Double_t delta = 0.01; /* UV cutoff (between 10^-6 and 10^-2) */
const Double_t R = 2.; /* IR cutoff for the Gaussian form (typ 2.) */
const Double_t y_max = 3.; /* maximum rapidity in evolution (typ 3.) */
const Int_t n_real = 200; /* number of independant events (typ 10000) */

/* fluctuations calculation parameters */
const Double_t r_s = 0.05; /* minimum size to consider (>~ delta) */

/* common ancestors parameters */
const Int_t k_leaves = 4; /* compute the common ancestor of k random leaves (typ 2-10) */
const Int_t m_factor = 2; /* consider only events with a multiplicity greater than mfactor * mean multiplicity (typ 2) */


class Dipole {
public:
    Double_t px;
    Double_t py;
    Double_t vx;
    Double_t vy;
    Double_t y;
    
    Dipole(Double_t px, Double_t py, Double_t vx, Double_t vy, Double_t y);
    ~Dipole();
    
    Dipole * split();
    Dipole * invertWithRespectTo(Dipole * p);
    Double_t size();
    
private:
    static TF1 * f1;
    static TF1 * f2;
    static TF1 * g;
    static TF1 * t1;
    static TF1 * t2;
    static TF1 * r;
    
};

/* random generator in )0,1) */
TRandom * sim = new TRandom;

/* generating function for the creation of leaves (r < 1/2) */
TF1 * Dipole::f1 = new TF1("f(r) pour r<=1/2","TMath::Pi()/(x*(1.-x**2))");

/* generating function for the creation of leaves (r >= 1/2) */
TF1 * Dipole::f2 = new TF1("f(r) pour r>1/2","2./(x*abs(1.-x**2))*atan(abs(1.-x)/(1.+x)*sqrt((x+0.5)/(x-0.5)))");

/* rejecting function */
TF1 * Dipole::g = new TF1("g(r)","5.*TMath::Pi()/(3.*x*(1.+x**2))");

/* generating function for the angle of leaves (r < 1/2) */
TF1 * Dipole::t1 = new TF1("theta(r) pour r<=1/2","2.*atan((1.-[0])/(1.+[0])*1./(tan(x*TMath::Pi()/2.)))");

/* generating function for the angle of leaves (r >= 1/2) */
TF1 * Dipole::t2 = new TF1("theta(r) pour r>1/2","2.*atan(abs(1.-[0])/(1.+[0])*1./(tan(x*atan(abs(1.-[0])/(1.+[0])*sqrt(([0]+0.5)/([0]-0.5))))))");

/* generating function for the radius of the leaves */
TF1 * Dipole::r = new TF1("r(U(0,1))","1./sqrt(pow((1.+1./[0]**2),1.-x)-1.)");


Dipole * Dipole::invertWithRespectTo(Dipole * p) {
    return new Dipole(this->px+this->vx,this->py+this->vy,p->vx-this->vx,p->vy-this->vy,this->y);
}

Dipole::Dipole(Double_t px, Double_t py, Double_t vx, Double_t vy, Double_t y) {
    
    this->px = px;
    this->py = py;
    this->vx = vx;
    this->vy = vy;
    this->y = y;
    
}

Double_t Dipole::size() {
    
    return sqrt(vx*vx+vy*vy);
    
}

Dipole * Dipole::split() {
    
    Double_t lambda;
    
    r->FixParameter(0,delta/this->size());
    
    if (delta/this->size() <= 0.5)
        lambda = log(1./3.) + log((this->size()*this->size())/(delta*delta) - 1.) + (2./TMath::Pi())*Dipole::f2->Integral(0.5,TMath::Infinity(),0.001);
    else
        lambda = (2./TMath::Pi())*Dipole::f2->Integral(delta/this->size(),TMath::Infinity(),0.001);
    
    /* TODO : recompute lambda with IR cutoff */
    
    Double_t radius, ratio, theta;
    
    int quadrant;
    
    Bool_t gen = false;
    
    while(!gen) {
        
        radius = r->Eval(1-sim->Rndm());
        
        if (radius <= 0.5) ratio = (f1->Eval(radius))/(g->Eval(radius));
        else ratio = (f2->Eval(radius))/(g->Eval(radius));
        
        if (sim->Rndm() <= ratio) {
            
            /* comment the follwing line for IR cutoff */
            gen = true;
            
            if (radius <= .5) {
                t1->FixParameter(0,radius);
                theta = t1->Eval(1-sim->Rndm());
            } else {
                t2->FixParameter(0,radius);
                theta = t2->Eval(1-sim->Rndm());
            }
            
            //if (sim->Rndm() <= exp(-(radius*radius+1-radius*radius-2*radius*cos(theta))/(2*R*R))) gen = true;
            
            Double_t phi = atan2(this->vy,this->vx);
            radius *= this->size();
            
            Double_t vxd, vyd;
            
            quadrant = (int)4*(1-sim->Rndm());
            
            switch (quadrant) {
                case 3:
                    vxd = this->vx-radius*cos(phi-theta);
                    vyd = this->vy-radius*sin(phi-theta);
                    break;
                case 2:
                    vxd = this->vx-radius*cos(phi+theta);
                    vyd = this->vy-radius*sin(phi+theta);
                    break;
                case 1:
                    vxd = radius*cos(phi-theta);
                    vyd = radius*sin(phi-theta);
                    break;
                default:
                    vxd = radius*cos(phi+theta);
                    vyd = radius*sin(phi+theta);
                    break;
            }
            
            Double_t rapidity = this->y - log(1.-sim->Rndm())/lambda;
            
            return new Dipole(this->px,this->py,vxd,vyd,rapidity);
            
        }
        
    }
    
    return NULL;
    
}

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

TF1 * fit_Bessel = new TF1("fitBessel", "[0]*TMath::BesselI0(2.*sqrt(-2.*[1]*log(x)))");
TF1 * P_n = new TF1("P_n", "1./x * 1./([0]*[0]) * exp(-log(x)*log(x)/(4*[1]))", 3, 50);
TF1 * n_bar = new TF1("n_bar", "1./x * exp([0]*4.*log(2.)) * exp(-log(x*x)*log(x*x)/([0]*56.*[1])) * 1./sqrt([0]*14.*TMath::Pi()*[1])");

//TH1 * h_size = new TH1D("h_size", "Size of leaves", 10000, -3, 2);
//TH1 * h_rap = new TH1D("h_rap", "Rapidity of leaves", 10000, 0, ymax);
TH1 * h_anc = new TH1D("h_anc", TString::Format("Size of ancestors (N = %d, y_max = %.12g, delta = %.12g, k = %d)", n_real, y_max, delta, k_leaves), 50, -2, 2);
TH1 * h_fluct = new TH1I("h_fluct", TString::Format("Fluctuations of multiplicity (N = %d, y_max = %.12g, delta = %.12g)", n_real, y_max, delta), 100, 0, 100);


void buildTree() {
    
    Dipole * base = new Dipole(0.,0.,1.,0.,0.);
    
    Double_t px,py,vx,vy,y;
    Int_t ancestor;
    Double_t size;
    
    //TTree * leaves = new TTree("leaves", "Tree containing leaves");
    //leaves->Branch("size",&size,"size/D");
    //leaves->Branch("rapidity",&y,"rapidity/D");
    //leaves->Branch("ancestor",&ancestor,"ancestor/I");
    
    //TTree * ancestors = new TTree("ancestors", "Tree containing ancestors");
    //ancestors->Branch("size",&size,"size/D");
    //ancestors->Branch("rapidity",&y,"rapidity/D");
    
    Int_t from, to;
    
    Bool_t complete;
    Dipole * p, * d;
    
    Int_t result[k_leaves];
    Int_t r;
    
    Int_t num;
    
    n_bar->FixParameter(0, y_max);
    n_bar->FixParameter(1, zeta_3);
    
    fit_Bessel->SetParameter(0,n_real);
    fit_Bessel->FixParameter(1,y_max);
    
    //BinLogX(h_size);
    BinLogX(h_anc);
    
    TCanvas *c1 = new TCanvas("Modèle des dipoles");
    
    //TFile f("tree.root","recreate");
    
    for (Int_t n = 0; n < n_real; n++) {
        
        printf("%f %%\n",100.*n/n_real);
        
        TTree * t = new TTree("t","Binary tree of evolution");
        t->Branch("px",&px,"px/D");
        t->Branch("py",&py,"py/D");
        t->Branch("vx",&vx,"vx/D");
        t->Branch("vy",&vy,"vy/D");
        t->Branch("size",&size,"size/D");
        t->Branch("rapidity",&y,"rapidity/D");
        t->Branch("ancestor",&ancestor,"ancestor/I");
        
        t->SetBranchAddress("px",&px);
        t->SetBranchAddress("py",&py);
        t->SetBranchAddress("vx",&vx);
        t->SetBranchAddress("vy",&vy);
        t->SetBranchAddress("size",&size);
        t->SetBranchAddress("rapidity",&y);
        t->SetBranchAddress("ancestor",&ancestor);
        
        px = base->px;
        py = base->py;
        vx = base->vx;
        vy = base->vy;
        y = base->y;
        size = base->size();
        ancestor = 0;
        
        t->Fill();
        
        TTree * l = new TTree("l", "Tree containing leaves");
        l->Branch("px",&px,"px/D");
        l->Branch("py",&py,"py/D");
        l->Branch("vx",&vx,"vx/D");
        l->Branch("vy",&vy,"vy/D");
        l->Branch("size",&size,"size/D");
        l->Branch("rapidity",&y,"rapidity/D");
        l->Branch("ancestor",&ancestor,"ancestor/I");
        
        complete = false;
        from = 0;
        to = 0;
        
        while (!complete) {
            
            complete = true;
            
            to = t->GetEntries();
            
            for (Int_t i = from; i < to; i++) {
                
                complete = false;
                
                t->GetEvent(i);
                
                p = new Dipole(px,py,vx,vy,y);
                d = p->split();
                
                if (d->y > y_max) {
                    
                    l->Fill();
                    //leaves->Fill();
                    //h_size->Fill(size);
                    //h_rap->Fill(y);
                    
                } else {
                    
                    px = d->px;
                    py = d->py;
                    vx = d->vx;
                    vy = d->vy;
                    y = d->y;
                    size = d->size();
                    ancestor = i;
                    t->Fill();
                    
                    d = d->invertWithRespectTo(p);
                    
                    px = d->px;
                    py = d->py;
                    vx = d->vx;
                    vy = d->vy;
                    size = d->size();
                    t->Fill();
                    
                }
                
            }
            
            from = to;
            
        }
        
        //printf("n_bar = %f, n = %lli\n",n_bar->Eval(delta),l->GetEntries());
        
        if (l->GetEntries() >= k_leaves && l->GetEntries() >= m_factor*n_bar->Eval(delta)) {
            
            for(Int_t i = 0; i < k_leaves; i++) {
                
                result[i] = l->GetEntries();
                
                do {
                    
                    r = (int)l->GetEntries()*(1-sim->Rndm());
                    
                } while (std::find(result, result+i, r) != result+i);
                
                result[i] = r;
            }
            
            for (Int_t i = 0; i < k_leaves - 1; i++) {
                
                while (result[i] != result[i+1]) {
                    
                    if (result[i] < result[i+1]) {
                        
                        t->GetEvent(result[i+1]);
                        result[i+1] = ancestor;
                        
                    } else {
                        
                        t->GetEvent(result[i]);
                        result[i] = ancestor;
                        
                    }
                    
                }
                
            }
            
            //ancestors->Fill();
            h_anc->Fill(size);
            
            //printf("%f %i\n",size,l->GetEntries());
            
        }
        
        num = l->GetEntries(TString::Format("size >= %.12g", r_s));
        h_fluct->Fill(num);
        
        delete t;
        delete l;
        
        
    }
    
    //t.Write();
    //f.Close();
    
    //c1->Divide(2,2);
    
    //c1->cd(1);
    //gPad->SetLogx();
    //gPad->SetLogy();
    //h_size = h_size->GetCumulative(kFALSE);
    //h_size->Draw();
    //gStyle->SetOptFit(1011);
    //h_size->Fit("fitBessel", "", "", delta, 0.1);
    
    //c1->cd(3);
    h_anc->GetXaxis()->SetTitle("Radius r");
    h_anc->GetYaxis()->SetTitle("Number of ancestors of radius r");
    gPad->SetLogx();
    //gPad->SetLogy();
    h_anc->Draw();
    //ancestors->Draw("size",NULL);
    
    
    //c1->cd(4);
    //gPad->SetLogx();
    //gPad->SetLogy();
    //h_fluct->Draw();
    //P_n->FixParameter(0, r_s);
    //P_n->FixParameter(1, y_max);
    //h_fluct->Fit("P_n", "", "", 6, 100);
    
}
