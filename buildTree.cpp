/* constants */
const Double_t zeta_3 = 1.2020569031595942854;

/* general parameters */
const Double_t delta = 0.001; /* UV cutoff (between 10^-6 and 10^-2) */
const Double_t y_max = 2.; /* maximum rapidity in evolution (typ 3.) */
const UInt_t N = 1000; /* number of independant events (typ 10000) */
const Double_t precision = 1./sqrt(N); /* for integration */

/* IR cutoff parameters */
const Bool_t IR = false; /* activate IR cutoff */
const Double_t R = 2.; /* IR cutoff for the Gaussian form (typ 2.) */

/* fluctuations calculation parameters */
const Bool_t compute_fluctuations = true; /* compute fluctutations */
const Double_t r_s = 0.005; /* minimum size to consider (~ delta) */

/* common ancestors parameters */
const Bool_t compute_ancestors = false; /* compute common ancestors ? */
const UInt_t k_leaves = 2; /* compute the common ancestor of k random leaves (typ 2-10) */
const UInt_t m_factor = 5; /* consider only events with a multiplicity greater than mfactor * mean multiplicity (typ 2) */


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
    Double_t lambda();
    
private:
    static TF1 * f1;
    static TF1 * f2;
    static TF1 * g;
    static TF1 * t1;
    static TF1 * t2;
    static TF1 * r;
    static TF2 * l;
    static Double_t l1(Double_t *x,Double_t *parm);
    static TF1 * lr;
};

/* random generator in )0,1) */
TRandom * sim = new TRandom;

/* generating function for the creation of leaves (r < 1/2) */
TF1 * Dipole::f1 = new TF1("f(r) pour r<=1/2","TMath::Pi()/(x*(1.-x**2))");

/* generating function for the creation of leaves (r >= 1/2) */
TF1 * Dipole::f2 = new TF1("f(r) pour r>1/2","2./(x*abs(1.-x**2))*atan(abs(1.-x)/(1.+x)*sqrt((x+.5)/(x-.5)))");

/* rejecting function */
TF1 * Dipole::g = new TF1("g(r)","5.*TMath::Pi()/(3.*x*(1.+x**2))");

/* generating function for the angle of leaves (r < 1/2) */
TF1 * Dipole::t1 = new TF1("theta(r) pour r<=1/2","2.*atan((1.-[0])/(1.+[0])*1./(tan(x*TMath::Pi()/2.)))");

/* generating function for the angle of leaves (r >= 1/2) */
TF1 * Dipole::t2 = new TF1("theta(r) pour r>1/2","2.*atan(abs(1.-[0])/(1.+[0])*1./(tan(x*atan(abs(1.-[0])/(1.+[0])*sqrt(([0]+.5)/([0]-.5))))))");

/* generating function for the radius of the leaves */
TF1 * Dipole::r = new TF1("r(U(0,1))","1./sqrt(pow((1.+1./[0]**2),1.-x)-1.)");

/* to-integrate lambda function */
TF2 * Dipole::l = new TF2("lambda","exp(-(x*x + 1. + x*x - 2.*x*cos(y))/(2.*[0]*[0]))*1./(x*(1.+x*x-2.*x*cos(y)))");

Double_t Dipole::l1(Double_t *x, Double_t *parm) {
    
    l->FixParameter(0,R);
    TF12 * l12 = new TF12("l12",l,x[0],"y");
    if (x[0] > .5) return l12->Integral(acos(1./(2.*x[0])),TMath::Pi(),precision);
    else return l12->Integral(0.,TMath::Pi(),precision);
    
}

TF1 * Dipole::lr = new TF1("lr",l1,0.,TMath::Infinity(),0);

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

Double_t Dipole::lambda() {
    
    if (!IR) {
        if (delta/this->size() <= .5)
            return log(1./3.) + log((this->size()*this->size())/(delta*delta) - 1.) + 2./TMath::Pi() * f2->Integral(.5,TMath::Infinity(),precision);
        else
            return 2./TMath::Pi() * f2->Integral(delta/this->size(),TMath::Infinity(),precision);
    }
    else return 2./TMath::Pi() * lr->Integral(delta/this->size(),TMath::Infinity(),precision);
    
}

Dipole * Dipole::split() {
    
    r->FixParameter(0,delta/this->size());
    
    Double_t radius, ratio, theta;
    
    int quadrant;
    
    Bool_t gen = false;
    
    while(!gen) {
        
        radius = r->Eval(1-sim->Rndm());
        
        if (radius <= .5) ratio = (f1->Eval(radius))/(g->Eval(radius));
        else ratio = (f2->Eval(radius))/(g->Eval(radius));
        
        if (sim->Rndm() <= ratio) {
            
            if (radius <= .5) {
                t1->FixParameter(0,radius);
                theta = t1->Eval(1-sim->Rndm());
            } else {
                t2->FixParameter(0,radius);
                theta = t2->Eval(1-sim->Rndm());
            }
            
            if (!IR) gen = true;
            else if (sim->Rndm() <= exp(-(radius*radius + 1. + radius*radius - 2.*radius*cos(theta))/(2.*R*R))) gen = true;
            
            if (gen) {
                
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
                
                Double_t rapidity = this->y - log(1.-sim->Rndm())/this->lambda();
                
                return new Dipole(this->px,this->py,vxd,vyd,rapidity);
                
            }
            
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
TF1 * P_n = new TF1("P_n", "1./x * [3]*[2]/([0]*[0]) * exp(-log(x)*log(x)/(4*[1]))");
TF1 * n_bar = new TF1("n_bar", "1./x * exp([0]*4.*log(2.)) * exp(-log(x*x)*log(x*x)/([0]*56.*[1])) * 1./sqrt([0]*14.*TMath::Pi()*[1])");

void buildTree() {
    
    Dipole * base = new Dipole(0.,0.,1.,0.,0.);
    
    Double_t px,py,vx,vy,y;
    UInt_t ancestor;
    Double_t size;
    
    UInt_t from, to;
    
    Bool_t complete;
    Dipole * p, * d;
    
    UInt_t random_leaves[k_leaves];
    UInt_t r;
    
    UInt_t num;
    
    n_bar->FixParameter(0, y_max);
    n_bar->FixParameter(1, zeta_3);
    
    fit_Bessel->SetParameter(0,N);
    fit_Bessel->FixParameter(1,y_max);
    
    P_n->FixParameter(0, r_s);
    P_n->FixParameter(1, y_max);
    P_n->FixParameter(2, N);
    
    //TH1 * h_size = new TH1D("h_size", "Size of leaves", 10000, 0, 2);
    //TH1 * h_rap = new TH1D("h_rap", "Rapidity of leaves", 10000, 0, y_max);
    TH1 * h_anc = new TH1D("h_anc", TString::Format("Size of ancestors (N = %d, y_max = %.12g, delta = %.12g, k = %d)", N, y_max, delta, k_leaves), 50, -2, 2);
    TH1 * h_fluct = new TH1I("h_fluct", TString::Format("Fluctuations of multiplicity (N = %d, y_max = %.12g, delta = %.12g)", N, y_max, delta), 20, 0, 10*n_bar->Eval(delta));
    
    //BinLogX(h_anc);
    
    TCanvas *c1 = new TCanvas("Modèle des dipoles");
    
    //TFile f("tree.root","recreate");
    
    for (UInt_t n = 0; n < N; n++) {
        
        printf("%.12g %%\n",100.*n/N);
        
        TTree * t = new TTree("t","Binary tree of evolution");
        t->Branch("px",&px,"px/D");
        t->Branch("py",&py,"py/D");
        t->Branch("vx",&vx,"vx/D");
        t->Branch("vy",&vy,"vy/D");
        t->Branch("size",&size,"size/D");
        t->Branch("rapidity",&y,"rapidity/D");
        t->Branch("ancestor",&ancestor,"ancestor/i");
        
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
        l->Branch("ancestor",&ancestor,"ancestor/i");
        
        complete = false;
        from = 0;
        to = 0;
        
        while (!complete) {
            
            complete = true;
            
            to = t->GetEntries();
            
            for (UInt_t i = from; i < to; i++) {
                
                complete = false;
                
                t->GetEvent(i);
                
                p = new Dipole(px,py,vx,vy,y);
                d = p->split();
                
                if (d->y > y_max) {
                    
                    l->Fill();
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
        
        if (compute_ancestors) {
            
            if (l->GetEntries() >= k_leaves && l->GetEntries() >= m_factor*n_bar->Eval(delta)) {
                
                for(UInt_t k = 0; k < k_leaves; k++) {
                    
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
                
            }
            
        }
        
        if (compute_fluctuations) {
            
            num = l->GetEntries(TString::Format("size >= %.12g", r_s));
            h_fluct->Fill(num);
            
        }
        
        //t->Write();
        //l->Write();
        
    }
    
    
    //c1->Divide(2,2);
    
    /*
    //c1->cd(1);
    //gPad->SetLogx();
    gPad->SetLogy();
    //h_size = h_size->GetCumulative(kFALSE);
    h_size->Draw("E1");
    //gStyle->SetOptFit(1011);
    //h_size->Fit("fitBessel", "", "", delta, 0.1);
    
    //c1->cd(3);
    h_anc->GetXaxis()->SetTitle("Radius r");
    h_anc->GetYaxis()->SetTitle("Number of ancestors of radius r");
    gPad->SetLogx();
    //gPad->SetLogy();
    h_anc->Draw();
    //ancestors->Draw("size",NULL);
    */
    
    //c1->cd(4);
    gPad->SetLogy();
    h_fluct->Draw("E1");
    h_fluct->Fit("P_n", "", "", n_bar->Eval(delta), 10*n_bar->Eval(delta));
    
}
