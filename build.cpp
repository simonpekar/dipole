#include "constants.h"

#include <Math/Interpolator.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TRandom.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <TSystem.h>

using namespace std;

class Dipole {
public:
    Double_t px;
    Double_t py;
    Double_t vx;
    Double_t vy;
    Double_t y;
    
    Dipole(Double_t px, Double_t py, Double_t vx, Double_t vy, Double_t y);
    ~Dipole();
    
    Double_t size();
    void invertWithRespectTo(Dipole * p);
    Dipole * split();
};

Double_t lambda(Double_t x01);
Double_t getLambda(Double_t x01);

ROOT::Math::Interpolator interpolator;
map<Double_t, Double_t> lookup_table;
Double_t x01_min, x01_max;


Dipole::~Dipole(void) { }

/* random generator in )0,1) */
TRandom * sim = new TRandom;

/* generating function for the creation of leaves (r < 1/2) */
TF1 * f1 = new TF1("f(r) pour r<=1/2","TMath::Pi()/(x*(1.-x**2))");

/* generating function for the creation of leaves (r >= 1/2) */
TF1 * f2 = new TF1("f(r) pour r>1/2","2./(x*abs(1.-x**2))*atan(abs(1.-x)/(1.+x)*sqrt((x+.5)/(x-.5)))");

/* rejecting function */
TF1 * g = new TF1("g(r)","5.*TMath::Pi()/(3.*x*(1.+x**2))");

/* generating function for the angle of leaves (r < 1/2) */
TF1 * t1 = new TF1("theta(r) pour r<=1/2","2.*atan((1.-[0])/(1.+[0])*1./(tan(x*TMath::Pi()/2.)))");

/* generating function for the angle of leaves (r >= 1/2) */
TF1 * t2 = new TF1("theta(r) pour r>1/2","2.*atan(abs(1.-[0])/(1.+[0])*1./(tan(x*atan(abs(1.-[0])/(1.+[0])*sqrt(([0]+.5)/([0]-.5))))))");

/* generating function for the radius of the leaves */
TF1 * r = new TF1("r(U(0,1))","1./sqrt(pow((1.+1./[0]**2),1.-x)-1.)");

/* to integrate lambda function */
TF2 * l = new TF2("lambda","1./(x*(1.+x*x-2.*x*cos(y)))");

/* cut-off types */
TF2 * rigid = new TF2("rigid", "(x*x + 1. + x*x - 2.*x*cos(y) < 2.*[0]*[0])");                      /* type 1 */
TF2 * gaussian = new TF2("gaussian","exp(-(x*x + 1. + x*x - 2.*x*cos(y))/(2.*[0]*[0]))");           /* type 2 */
TF2 * exponential = new TF2("exponential","exp(-(sqrt(x*x + 1. + x*x - 2.*x*cos(y)))/([0]))");      /* type 3 */


Double_t l1(Double_t *x, Double_t *parm) {
    
    Double_t result;
    TF2 * mult;
    
    switch (IR_type) {
        case 0:
            mult = new TF2("mult","lambda");
            break;
            
        case 1:
            rigid->FixParameter(0,R);
            mult = new TF2("mult","lambda*rigid");
            break;
            
        case 2:
            gaussian->FixParameter(0,R);
            mult = new TF2("mult","lambda*gaussian");
            break;
            
        case 3:
            exponential->FixParameter(0,R);
            mult = new TF2("mult","lambda*exponential");
            break;
            
        default:
            cerr << "unknown cut-off" << endl;
            break;
    }
    
    TF12 * l12 = new TF12("l12",mult,x[0],"y");
    
    if (x[0] > .5) result = l12->Integral(acos(1./(2.*x[0])),TMath::Pi());
    else result = l12->Integral(0.,TMath::Pi(),1.);
    
    delete mult;
    delete l12;
    
    return result;
    
}

TF1 * lr = new TF1("lr",l1,0.,TMath::Infinity(),0);

void Dipole::invertWithRespectTo(Dipole * p) {
    
    Double_t px2 = this->px+this->vx;
    Double_t py2 = this->py+this->vy;
    Double_t vx2 = p->vx-this->vx;
    Double_t vy2 = p->vy-this->vy;
    this->px = px2;
    this->py = py2;
    this->vx = vx2;
    this->vy = vy2;
    
}

Dipole::Dipole(Double_t px, Double_t py, Double_t vx, Double_t vy, Double_t y) {
    
    this->px = px;
    this->py = py;
    this->vx = vx;
    this->vy = vy;
    this->y = y;
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
            
            
            switch (IR_type) {
                case 0:
                    gen = true;
                    break;
                    
                case 1:
                    rigid->FixParameter(0,R);
                    gen = (rigid->Eval(radius,theta) != 0.);
                    
                case 2:
                    gaussian->FixParameter(0,R);
                    gen = (sim->Rndm() <= gaussian->Eval(radius,theta));
                    
                case 3:
                    exponential->FixParameter(0,R);
                    gen = (sim->Rndm() <= exponential->Eval(radius,theta));
                    
                default:
                    break;
            }
            
            if (gen) {
                
                radius *= this->size();
                
                Double_t vxd, vyd;
                
                quadrant = (int)4*(1-sim->Rndm());
                
                Double_t phi = atan2(this->vy,this->vx);
                    
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
                
                Double_t rapidity = this->y - log(1.-sim->Rndm())/getLambda(delta/this->size());
                
                return new Dipole(this->px,this->py,vxd,vyd,rapidity);
                
                
            }
            
        }
        
    }
    
    return NULL;
    
}

Double_t lambda(Double_t x01) {
    
    Double_t result;
    
    if (IR_type == 0) {
        
        if (x01 <= .5) result = log(1./3.) + log(1/(x01*x01) - 1.) + 2./TMath::Pi() * f2->Integral(.5,TMath::Infinity());
        else result = 2./TMath::Pi() * f2->Integral(x01,TMath::Infinity());
        
    } else {
        
        result = 2./TMath::Pi() * lr->Integral(x01,TMath::Infinity(),0.001);
        
    }
    
    return result;
    
}

void writeLookupTable()
{
    int n = 200;
    
    Double_t x01 = 0.0;
    Double_t step = 0.000000000001;
    Double_t l;
    ofstream lut("table");
    if (lut.is_open())
    {
        lut << delta << " " << n << "\n";
        for (int i = 1; i <= n; ++i)
        {
            if (i%10 == 0) step *= 10. ;
            x01 += step;
            l = lambda(x01);
            lut << x01 << " " << l << "\n";
        }
        lut.close();
    }
}


void loadLookupTable()
{
    cout << "Reading lookup table..." << endl;
    ifstream lut("table");
    Double_t rho;
    int n;
    vector<Double_t> x01, l;
    if (lut.is_open())
    {
        lut >> rho >> n;
        Double_t a, b;
        while(lut >> a >> b)
        {
            lookup_table.insert(pair<Double_t, Double_t>(a, b));
            x01.push_back(a);
            l.push_back(b);
        }
        lut.close();
    }
    interpolator.SetData(x01, l);
    x01_min = x01.front();
    x01_max = x01.back();
    cout << "\033[1;32m Done. \033[0m" << endl;
}

void printLookupTable()
{
    cout << "Printing lookup table..." << endl;
    for (auto const& x: lookup_table)
    {
        cout << x.first << " " << x.second << endl;
    }
    cout << "done." << endl;
}

void setInterpolatorData()
{
    vector<Double_t> x01, l;
    for (auto const& x: lookup_table)
    {
        x01.push_back(x.first);
        l.push_back(x.second);
    }
    interpolator.SetData(x01, l);
}

Double_t getLambda(Double_t x01)
{
    if (lookup_table.size() <= 0)
    {
        loadLookupTable();
    }
    // Interpolate
    if (x01_min > x01 || x01 > x01_max)
    {
        cout << "\033[1;31mWarning : x01 out of lookup table range.\033[0m \n Adding entry for " << x01 << endl;
        // Add new points
        lookup_table.insert(pair<Double_t, Double_t>(x01, lambda(x01)));
        x01_min = TMath::Min(x01_min, x01);
        x01_max = TMath::Max(x01_max, x01);
        setInterpolatorData();
    }
    //cout << delta << " " << x01 << " " << interpolator.Eval(x01) << endl;
    return interpolator.Eval(x01);
}


void build() {
    
    writeLookupTable();
    loadLookupTable();
    printLookupTable();
    
    Dipole * base = new Dipole(0.,0.,1.,0.,0.);
    
    Double_t px,py,vx,vy;
    UInt_t ancestor;
    Double_t size, y;
    
    UInt_t from, to;
    Bool_t complete;
    
    Dipole * p, * d;
    
    TFile * f;
    
    switch (IR_type) {
            
        case 0:
            f = new TFile(TString::Format("datasets/y=%.12g, N=%d, UV=%.12g, no IR.root",y_max,N,delta),"recreate");
            break;
            
        case 1:
            f = new TFile(TString::Format("datasets/y=%.12g, N=%d, UV=%.12g, IR=%.12g rigid.root",y_max,N,delta,R),"recreate");
            break;
            
        case 2:
            f = new TFile(TString::Format("datasets/y=%.12g, N=%d, UV=%.12g, IR=%.12g gaussian.root",y_max,N,delta,R),"recreate");
            break;
            
        case 3:
            f = new TFile(TString::Format("datasets/y=%.12g, N=%d, UV=%.12g, IR=%.12g exponential.root",y_max,N,delta,R),"recreate");
            break;
            
        default:
            cerr << "unknown cut-off" << endl;
            break;
            
    }
    
    for (UInt_t n = 0; n < N; n++) {
        
        
        cout << 100.*n/N << "%" << endl;
        
        TTree * t = new TTree(TString::Format("t%d",n),"evolution");
        t->Branch("px",&px,"px/D");
        t->Branch("py",&py,"py/D");
        t->Branch("vx",&vx,"vx/D");
        t->Branch("vy",&vy,"vy/D");
        t->Branch("size",&size,"size/D");
        t->Branch("rapidity",&y,"rapidity/D");
        t->Branch("ancestor",&ancestor,"ancestor/i");
        
        px = base->px;
        py = base->py;
        vx = base->vx;
        vy = base->vy;
        y = base->y;
        size = base->size();
        ancestor = 0;
        
        t->Fill();
        
        TTree * l = new TTree(TString::Format("l%d",n), "leaves");
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
                    
                } else {
                    
                    px = d->px;
                    py = d->py;
                    vx = d->vx;
                    vy = d->vy;
                    y = d->y;
                    size = d->size();
                    ancestor = i;
                    t->Fill();
                    
                    d->invertWithRespectTo(p);
                    
                    px = d->px;
                    py = d->py;
                    vx = d->vx;
                    vy = d->vy;
                    size = d->size();
                    t->Fill();
                    
                }
                
                delete p;
                delete d;
                
            }
            
            from = to;
            
        }
        
        //t->Write();
        l->Write();
        
    }
    
    f->Close();
    
}
