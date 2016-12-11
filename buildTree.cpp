const Double_t rho = 0.000001; // cutoff UV
const Double_t R = 2.; // cutoff IR
const Double_t ymax = 0.35;
const Double_t rs = 0.05;
const Int_t nreal = 10000;

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

Dipole * Dipole::invertWithRespectTo(Dipole * p) {
    return new Dipole(this->px+this->vx,this->py+this->vy,p->vx-this->vx,p->vy-this->vy,this->y);
}

// générateur de hasard
TRandom * sim = new TRandom;
// f pour r < 1/2
TF1 * Dipole::f1 = new TF1("f(r) pour r<=1/2","TMath::Pi()/(x*(1.-x**2))");
// f pour r >= 1/2
TF1 * Dipole::f2 = new TF1("f(r) pour r>1/2","2./(x*abs(1.-x**2))*atan(abs(1.-x)/(1.+x)*sqrt((x+0.5)/(x-0.5)))");
// g (rejet)
TF1 * Dipole::g = new TF1("g(r)","5.*TMath::Pi()/(3.*x*(1.+x**2))");
// theta pour r < 1/2
TF1 * Dipole::t1 = new TF1("theta(r) pour r<=1/2","2.*atan((1.-[0])/(1.+[0])*1./(tan(x*TMath::Pi()/2.)))");
// theta pour r >= 1/2
TF1 * Dipole::t2 = new TF1("theta(r) pour r>1/2","2.*atan(abs(1.-[0])/(1.+[0])*1./(tan(x*atan(abs(1.-[0])/(1.+[0])*sqrt(([0]+0.5)/([0]-0.5))))))");
// génération du rayon r
TF1 * Dipole::r = new TF1("r(U(0,1))","1./sqrt(pow((1.+1./[0]**2),1.-x)-1.)");


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
    
    r->FixParameter(0,rho/this->size());
    
    if (rho/this->size() <= 0.5)
        lambda = log(1./3.) + log((this->size()*this->size())/(rho*rho) - 1.) + (2./TMath::Pi())*Dipole::f2->Integral(0.5,TMath::Infinity(),0.001);
    else
        lambda = (2./TMath::Pi())*Dipole::f2->Integral(rho/this->size(),TMath::Infinity(),0.001);
    
    
    Double_t radius, ratio, theta;
    
    int quadrant;
    
    Bool_t gen = false;
    
    while(!gen) {
        
        radius = r->Eval(1-sim->Rndm());
        
        if (radius <= 0.5) ratio = (f1->Eval(radius))/(g->Eval(radius));
        else ratio = (f2->Eval(radius))/(g->Eval(radius));
        
        //if (sim->Rndm() <= ratio && sim->Rndm() <= exp(-(radius*radius)/(2*R*R))) {
        if (sim->Rndm() <= ratio) {
            // on le garde
            gen = true;
            
            if (radius <= .5) {
                t1->FixParameter(0,radius);
                theta = t1->Eval(1-sim->Rndm());
            } else {
                t2->FixParameter(0,radius);
                theta = t2->Eval(1-sim->Rndm());
            }
            
            Double_t phi = atan2(this->vy,this->vx);
            radius *= this->size();
            
            Double_t vxd, vyd;
            
            // quadrant : 0, 1, 2 ou 3
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


TF1 * fit1 = new TF1("fit1","[0]*exp(2.*sqrt(-2.*[1]*log(x)))");
TF1 * fit2 = new TF1("fit2","[0]/x*exp(4.*log(2.)*[1])");
TF1 * fit3 = new TF1("fit3","[0]/(x*x)*exp(2.*sqrt(2.*[1]*log(x)))");
TF1 * fitBessel = new TF1("fitBessel", "[0]*TMath::BesselI0(2.*sqrt(-2.*[1]*log(x)))");
TF1 * Pn = new TF1("Pn", "1./x * [0]/([1]*[1]) * exp(-log(x)*log(x)/(4*[2]))", 3, 50);


TH1 * h1 = new TH1D("h1", "Size of leaves", 10000, -7, 2);

TH1 * h2 = new TH1D("h2", "Rapidity of leaves", 10000, 0, ymax);

TH1 * h3 = new TH1D("h3", "Size of ancestors", 10000, -7, 2);

TH1 * hfluct = new TH1I("hfluct", "Fluctuations", 100, 0, 100);


void buildTree() {
    
    Dipole * base = new Dipole(0.,0.,1.,0.,0.);
    
    Double_t px,py,vx,vy,y;
    Int_t ancestor;
    Double_t size;
    
    TTree * leaves = new TTree("leaves", "Tree containing leaves");
    leaves->Branch("size",&size,"size/D");
    leaves->Branch("rapidity",&y,"rapidity/D");
    leaves->Branch("ancestor",&ancestor,"ancestor/I");
    
    TTree * ancestors = new TTree("ancestors", "Tree containing ancestors");
    ancestors->Branch("size",&size,"size/D");
    ancestors->Branch("rapidity",&y,"rapidity/D");
    
    Int_t from, to;
    
    Bool_t complete;
    Dipole * p, * d;
    Int_t leaf1, leaf2;
    Int_t num;
    
    fitBessel->SetParameter(0,nreal);
    fitBessel->FixParameter(1,ymax);
    
    BinLogX(h1);
    BinLogX(h3);
    
    TCanvas *c1 = new TCanvas("Test");
    
    //TFile f("tree.root","recreate");
    
    for (Int_t k = 0; k < nreal; k++) {
        
        printf("réalisation : %i\n",k);
        
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
                
                if (d->y > ymax) {
                    
                    l->Fill();
                    leaves->Fill();
                    h1->Fill(size);
                    h2->Fill(y);
                    
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
        
        leaf1 = (int)l->GetEntries()*(1-sim->Rndm());
        leaf2 = (int)l->GetEntries()*(1-sim->Rndm());
        
        if (l->GetEntries() >= 2) {
            while (leaf2 == leaf1) leaf2 = (int)l->GetEntries()*(1-sim->Rndm());
        }
        
        while (leaf1 != leaf2) {
            
            if (leaf1 < leaf2) {
                
                t->GetEvent(leaf2);
                leaf2 = ancestor;
                
            } else {
                
                t->GetEvent(leaf1);
                leaf1 = ancestor;
                
            }
            
        }
        
        ancestors->Fill();
        h3->Fill(size);
        num = l->GetEntries(TString::Format("size >= %.12g", rs));
        hfluct->Fill(num);
        
        delete t;
        delete l;
        
        
    }
    
    //t.Write();
    //f.Close();
    
    c1->Divide(2,2);
    
    c1->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    h1 = h1->GetCumulative(kFALSE);
    h1->Draw();
    gStyle->SetOptFit(1011);
    h1->Fit("fitBessel", "", "", rho, 0.1);
    
    c1->cd(2);
    gPad->SetLogy();
    h2->Draw();
    
    c1->cd(3);
    gPad->SetLogx();
    gPad->SetLogy();
    h3->Draw();
    //ancestors->Draw("size",NULL);
    
    c1->cd(4);
    gPad->SetLogx();
    gPad->SetLogy();
    hfluct->Draw(); 
    Pn->FixParameter(1, rs);
    Pn->FixParameter(2, ymax);
    hfluct->Fit("Pn", "", "", 6, 100);
    
}
