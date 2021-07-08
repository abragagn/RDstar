#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TF1.h"
#include "Math/GoFTest.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "Math/DistFunc.h"

using namespace std;
const TString filename = "decay_info.root";
const float pval = 0.05; // no cut: 0
const float tval = 0.5; // no cut: 1
const int   minEvt = 25; // gof test not working with less than 2 events
const int   midEvt = 50; // gof test not working with less than 2 events
const int   nbins = 1000; // fine binning for CDF distributions
const int   rebin = 100;  // rebin factor for plotting 
const bool  print = false;

void KSanalysis(int mode = 0) // 0 = double charm, 1 = single charm
{
    // gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    // Cleaning
    gSystem->Exec("rm -r plots");
    gSystem->Exec("mkdir plots");
    gSystem->Exec("mkdir plots/lowStat");
    gSystem->Exec("mkdir plots/bad");
    gSystem->Exec("mkdir plots/lowStat/bad");

    // Opening files
    auto f = new TFile(filename);
    auto t = (TTree*)f->Get("tree");
    t->AddFriend("decayIndexTree", filename);
    int ne = t->GetEntries();

    ifstream dictionary;
    dictionary.open("KSdictionary.txt");

    // Classes to compare (defined in decay_dictionary)
    vector<int> v_doubleHc{13, 14, 16, 17, 19, 20, 21, 25, 27, 29, 30, 31, 32, 33, 34, 35, 36, 38, 39, 40, 42, 43, 44, 45, 46, 48, 49, 50, 51 
                         , 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 68, 69, 70, 71, 72, 73, 74, 75, 76, 78, 79, 81, 82, 83
                         , 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 102, 103, 104, 105, 106, 107, 108, 109, 110
                         , 111, 112, 113, 114, 115, 116, 117, 118, 119, 120};

    vector<int> v_singleHc{1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 15, 18, 22, 23, 24, 26, 28, 37, 41, 47, 67, 77, 80, 101};

    vector<int> v_signal{0, 5};

    vector<int> v_classes;

    if(mode == 1) v_classes = v_singleHc;
    else          v_classes = v_doubleHc;

    int n_classes = v_classes.size();

    // Declare histograms (for plotting)
    TH1D *h_q2[n_classes];
    TH1D *h_m2_miss[n_classes];
    TH1D *h_e_star_mu3[n_classes];
    // TH2D *h2_q2_m2_miss[n_classes];
    // TH2D *h2_q2_e_star_mu3[n_classes];
    // TH2D *h2_m2_miss_e_star_mu3[n_classes];

    for(int i=0; i<n_classes; ++i){
        h_q2[i]         = new TH1D(TString::Format("h_q2_%i", i),           ";q^{2} [GeV^{2}];Events",          1400, 0., 14.);
        h_m2_miss[i]    = new TH1D(TString::Format("h_m2_miss_%i", i),      ";M^{2}_{miss} [GeV^{2}];Events",   1200, 0., 12.);
        h_e_star_mu3[i] = new TH1D(TString::Format("h_e_star_mu3_%i", i),   ";E_{#mu}* [GeV];Events",           1100, 0., 2.2);
        h_q2[i]->SetLineWidth(2);
        h_m2_miss[i]->SetLineWidth(2);
        h_e_star_mu3[i]->SetLineWidth(2);

        // h2_q2_m2_miss[i] = new TH2D(TString::Format("h2_q2_m2_miss%i", i), ";q^{2} [GeV^{2}];M^{2}_{miss} [GeV^{2}];Events", 100, 0., 12., 100, 0., 8.);
        // h2_q2_e_star_mu3[i] = new TH2D(TString::Format("h2_q2_e_star_mu3%i", i), ";q^{2} [GeV^{2}];E_{#mu}* [GeV];Events", 100, 0., 12., 100, 0., 2.2);
        // h2_m2_miss_e_star_mu3[i] = new TH2D(TString::Format("h2_m2_miss_e_star_mu3%i", i), ";M^{2}_{miss} [GeV^{2}];E_{#mu}* [GeV]Events", 100, 0., 8., 100, 0., 2.2);
    }

    // Declare vectors (fot gof tests)
    vector<double> v_q2[n_classes];
    vector<double> v_m2_miss[n_classes];
    vector<double> v_e_star_mu3[n_classes];

    // Produce dictionary vector
    TString park;
    TString dict[n_classes];
    int i_dict = 0;
    while(park.ReadToDelim(dictionary)){
        if(i_dict>=n_classes) break;
        dict[i_dict] = park;
        // cout<<dict[i_dict]<<endl;
        i_dict++;
    }


    // Filling
    float q2, m2_miss, e_star_mu3, has_double_charm, decayIndex;

    t->SetBranchAddress("q2", &q2);
    t->SetBranchAddress("m2_miss", &m2_miss);
    t->SetBranchAddress("e_star_mu3", &e_star_mu3);
    t->SetBranchAddress("has_double_charm", &has_double_charm); // currently not working as intended
    t->SetBranchAddress("decayIndex", &decayIndex);

    cout<<"Loop starting"<<endl;

    for(int i=0; i<ne; ++i){
        t->GetEntry(i);

        auto it = find(v_classes.begin(), v_classes.end(), int(decayIndex)); // find class position in v_classes, v.end() if not found
        if(it == v_classes.end()) continue; // skip events not beloging to one of the classes to study

        int index = distance(v_classes.begin(), it); // iterator -> index

        h_q2[index]->Fill(q2);
        h_m2_miss[index]->Fill(m2_miss);
        h_e_star_mu3[index]->Fill(e_star_mu3);
        // h2_q2_m2_miss[index]->Fill(q2,m2_miss);
        // h2_q2_e_star_mu3[index]->Fill(q2,e_star_mu3);
        // h2_m2_miss_e_star_mu3[index]->Fill(m2_miss,e_star_mu3);

        v_q2[index].push_back(q2);
        v_m2_miss[index].push_back(m2_miss);
        v_e_star_mu3[index].push_back(e_star_mu3);
    }

    // Get number of events for each class
    int v_ne[n_classes];
    for(int i=0; i<n_classes; ++i){
        v_ne[i] = h_q2[i]->GetEntries();
    }

    // Produce cumulative probability histograms
    TH1 *h_CDF_q2[n_classes];
    TH1 *h_CDF_m2_miss[n_classes];
    TH1 *h_CDF_e_star_mu3[n_classes];

    for(int i=0; i<n_classes; ++i){
        h_CDF_q2[i] = h_q2[i]->GetCumulative();
        h_CDF_m2_miss[i] = h_m2_miss[i]->GetCumulative();
        h_CDF_e_star_mu3[i] = h_e_star_mu3[i]->GetCumulative();

        h_CDF_q2[i]->Scale(1./v_ne[i]);
        h_CDF_m2_miss[i]->Scale(1./v_ne[i]);
        h_CDF_e_star_mu3[i]->Scale(1./v_ne[i]);

        // Cosmetics
        h_CDF_q2[i]->GetYaxis()->SetTitle("Cumulative Probability");
        h_CDF_m2_miss[i]->GetYaxis()->SetTitle("Cumulative Probability");
        h_CDF_e_star_mu3[i]->GetYaxis()->SetTitle("Cumulative Probability");

        h_CDF_q2[i]->SetMaximum(1.);
        h_CDF_m2_miss[i]->SetMaximum(1.);
        h_CDF_e_star_mu3[i]->SetMaximum(1.);

        h_CDF_q2[i]->SetMinimum(0.);
        h_CDF_m2_miss[i]->SetMinimum(0.);
        h_CDF_e_star_mu3[i]->SetMinimum(0.);
    }

    // Rebin distribution histograms (for plotting)
    for(int i=0; i<n_classes; ++i){
        h_q2[i]->Rebin(rebin);
        h_m2_miss[i]->Rebin(rebin);
        h_e_star_mu3[i]->Rebin(rebin);
    }

    // Compare distributions
    for(int i=0; i<n_classes; ++i){
        for(int j=i+1; j<n_classes; ++j){
            TString path = "plots/";

            if(v_ne[i] < minEvt || v_ne[j] < minEvt) continue;
            if(v_ne[i] < midEvt && v_ne[j] < midEvt) path += "lowStat/";

            // Create got tests (syntax: n evt1, array1, nevt2, array2)
            auto gof_q2         = new ROOT::Math::GoFTest(v_ne[i], &(v_q2[i])[0],           v_ne[j], &(v_q2[j])[0]);
            auto gof_m2_miss    = new ROOT::Math::GoFTest(v_ne[i], &(v_m2_miss[i])[0],      v_ne[j], &(v_m2_miss[j])[0]);
            auto gof_e_star_mu3 = new ROOT::Math::GoFTest(v_ne[i], &(v_e_star_mu3[i])[0],   v_ne[j], &(v_e_star_mu3[j])[0]);

            // KS pvalue (vs null hypothesis same distribution)
            double pKS_q2           = gof_q2->KolmogorovSmirnov2SamplesTest();
            double pKS_m2_miss      = gof_m2_miss->KolmogorovSmirnov2SamplesTest();
            double pKS_e_star_mu3   = gof_e_star_mu3->KolmogorovSmirnov2SamplesTest();

            // KS test statistics (aka max distrance between CDFs)
            double tKS_q2           = gof_q2->KolmogorovSmirnov2SamplesTest("t");
            double tKS_m2_miss      = gof_m2_miss->KolmogorovSmirnov2SamplesTest("t");
            double tKS_e_star_mu3   = gof_e_star_mu3->KolmogorovSmirnov2SamplesTest("t");

            double avg_pKS = (pKS_q2 + pKS_m2_miss + pKS_e_star_mu3)/3;
            double avg_tKS = (tKS_q2 + tKS_m2_miss + tKS_e_star_mu3)/3;

            bool good = avg_pKS >= pval && avg_tKS <= tval;

            if(!good) path += "bad/";

            if(print){
                cout<<endl<<"classes "<<v_classes[i]<<", "<<v_classes[j]<<" | nevt: "<<v_ne[i]<<", "<<v_ne[j];
                cout<<std::fixed<<std::setprecision(3)<<" | p-val: "<<pKS_q2<<" | KS-dist: "<<tKS_q2<<endl;
                cout<<"q2: "<<pKS_q2<<"  ("<<tKS_q2<<")"<<endl;
                cout<<"m2_miss: "<<pKS_m2_miss<<"  ("<<tKS_m2_miss<<")"<<endl;
                cout<<"e_star_mu3: "<<pKS_e_star_mu3<<"  ("<<tKS_e_star_mu3<<")"<<endl;
            }

            // PLOTTING
            auto c1 = new TCanvas("c1", "c1", 1600, 600);
            c1->Divide(2,2);

            c1->cd(1);
            gPad->Divide(2,1);
            gPad->cd(1);
            h_q2[i]->SetLineColor(kBlue);
            h_q2[j]->SetLineColor(kRed);
            h_q2[i]->Scale(1./h_q2[i]->Integral());
            h_q2[j]->Scale(1./h_q2[j]->Integral());
            h_q2[i]->GetYaxis()->SetTitleOffset(1.5);
            h_q2[i]->SetMaximum(1.2*max(h_q2[i]->GetBinContent(h_q2[i]->GetMaximumBin()), h_q2[j]->GetBinContent(h_q2[j]->GetMaximumBin())));
            h_q2[i]->Draw("HIST E");
            h_q2[j]->Draw("BE SAME");

            c1->cd(1);
            gPad->cd(2);
            h_CDF_q2[i]->SetLineColor(kBlue);
            h_CDF_q2[j]->SetLineColor(kRed);
            h_CDF_q2[i]->GetYaxis()->SetTitleOffset(1.2);
            h_CDF_q2[i]->Draw("HIST");
            h_CDF_q2[j]->Draw("HIST SAME");

            c1->cd(2);
            gPad->Divide(2,1);
            gPad->cd(1);
            h_m2_miss[i]->SetLineColor(kBlue);
            h_m2_miss[j]->SetLineColor(kRed);
            h_m2_miss[i]->Scale(1./h_m2_miss[i]->Integral());
            h_m2_miss[j]->Scale(1./h_m2_miss[j]->Integral());
            h_m2_miss[i]->GetYaxis()->SetTitleOffset(1.5);
            h_m2_miss[i]->SetMaximum(1.2*max(h_m2_miss[i]->GetBinContent(h_m2_miss[i]->GetMaximumBin()), h_m2_miss[j]->GetBinContent(h_m2_miss[j]->GetMaximumBin())));
            h_m2_miss[i]->DrawNormalized("HIST E");
            h_m2_miss[j]->DrawNormalized("BE SAME");

            c1->cd(2);
            gPad->cd(2);
            h_CDF_m2_miss[i]->SetLineColor(kBlue);
            h_CDF_m2_miss[j]->SetLineColor(kRed);
            h_CDF_m2_miss[i]->GetYaxis()->SetTitleOffset(1.2);
            h_CDF_m2_miss[i]->Draw("HIST");
            h_CDF_m2_miss[j]->Draw("HIST SAME");

            c1->cd(3);
            gPad->Divide(2,1);
            gPad->cd(1);
            h_e_star_mu3[i]->SetLineColor(kBlue);
            h_e_star_mu3[j]->SetLineColor(kRed);
            h_e_star_mu3[i]->Scale(1./h_e_star_mu3[i]->Integral());
            h_e_star_mu3[j]->Scale(1./h_e_star_mu3[j]->Integral());
            h_e_star_mu3[i]->GetYaxis()->SetTitleOffset(1.5);
            h_e_star_mu3[i]->SetMaximum(1.2*max(h_e_star_mu3[i]->GetBinContent(h_e_star_mu3[i]->GetMaximumBin()), h_e_star_mu3[j]->GetBinContent(h_e_star_mu3[j]->GetMaximumBin())));
            h_e_star_mu3[i]->DrawNormalized("HIST E");
            h_e_star_mu3[j]->DrawNormalized("BE SAME");

            c1->cd(3);
            gPad->cd(2);
            h_CDF_e_star_mu3[i]->SetLineColor(kBlue);
            h_CDF_e_star_mu3[j]->SetLineColor(kRed);
            h_CDF_e_star_mu3[i]->GetYaxis()->SetTitleOffset(1.2);
            h_CDF_e_star_mu3[i]->Draw("HIST");
            h_CDF_e_star_mu3[j]->Draw("HIST SAME");

            c1->cd(4);
            auto info = new TPaveText(.05,.1,.95,.8);
            info->AddText("Blue: " + dict[v_classes[i]] + TString::Format("   [%i]", v_classes[i]));
            info->AddText("Red: " + dict[v_classes[j]] + TString::Format("   [%i]", v_classes[j]));
            info->AddText(TString::Format("Number of entries(blue, red): (%i, %i)", v_ne[i], v_ne[j]));
            info->AddText(TString::Format("KS p-value (distance): %.3f (%.3f)", avg_pKS, avg_tKS));
            info->Draw();

            c1->Print(TString::Format(path + "plot_%i_%i.png", v_classes[i], v_classes[j]));

            // auto c2 = new TCanvas("c2","c2",800,600);
            // c2->Divide(2,2);
            // c2->cd(1);
            // h2_q2_m2_miss[i]->Draw();
            // c2->cd(2);
            // h2_q2_e_star_mu3[i]->Draw();
            // c2->cd(3);
            // h2_m2_miss_e_star_mu3[i]->Draw();
            // c2->Print(TString::Format("plots/2Dplot_%i.png", v_classes[i]));

            // double maxDist1 = 0.;
            // double maxDist2 = 0.;
            // double maxDist3 = 0.;

            // int maxDistbin1 = -1;
            // int maxDistbin2 = -1;
            // int maxDistbin3 = -1;

            // for(int ibin=1; ibin<=nbins; ++ibin){
            //     double dist1 = abs(h_CDF_q2[i]->GetBinContent(ibin) - h_CDF_q2[j]->GetBinContent(ibin));
            //     double dist2 = abs(h_CDF_m2_miss[i]->GetBinContent(ibin) - h_CDF_m2_miss[j]->GetBinContent(ibin));
            //     double dist3 = abs(h_CDF_e_star_mu3[i]->GetBinContent(ibin) - h_CDF_e_star_mu3[j]->GetBinContent(ibin));

            //     if(dist1>maxDist1){
            //         maxDist1 = dist1;
            //         maxDistbin1 = ibin;
            //     }

            //     if(dist2>maxDist2){
            //         maxDist2 = dist2;
            //         maxDistbin2 = ibin;
            //     }

            //     if(dist3>maxDist3){
            //         maxDist3 = dist3;
            //         maxDistbin3 = ibin;
            //     }
            // }
            // cout<<"Hand-made: "<<maxDist1<<" (bin: "<<maxDistbin1<<")"<<endl;
            // cout<<"Hand-made: "<<maxDist2<<" (bin: "<<maxDistbin2<<")"<<endl;
            // cout<<"Hand-made: "<<maxDist3<<" (bin: "<<maxDistbin3<<")"<<endl;

            delete info;
            delete c1;

            delete gof_q2;
            delete gof_m2_miss;
            delete gof_e_star_mu3;
        }
    }
}
