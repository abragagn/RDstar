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
const float pval = 0.10; // no cut: 0
const float tval = 0.20; // no cut: 1
const int   minEvt = 25; // gof test not working with less than 2 events
const int   midEvt = 50; // gof test not working with less than 2 events
const TString dir = "plots_grouping";

void KSgrouping()
{
    // gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    // Cleaning
    gSystem->Exec("rm -r " + dir);
    gSystem->Exec("mkdir " + dir);

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
    vector<int> v_classes = v_doubleHc;

    int n_classes = v_classes.size();

    // Declare histograms (for plotting)
    TH1D *h_q2[n_classes];
    TH1D *h_m2_miss[n_classes];
    TH1D *h_e_star_mu3[n_classes];

    for(int i=0; i<n_classes; ++i){
        h_q2[i]         = new TH1D(TString::Format("h_q2_%i", i),           ";q^{2} [GeV^{2}];Events",          14, 0., 14.);
        h_m2_miss[i]    = new TH1D(TString::Format("h_m2_miss_%i", i),      ";M^{2}_{miss} [GeV^{2}];Events",   12, 0., 12.);
        h_e_star_mu3[i] = new TH1D(TString::Format("h_e_star_mu3_%i", i),   ";E_{#mu}* [GeV];Events",           11, 0., 2.2);
        h_q2[i]->SetLineWidth(2);
        h_m2_miss[i]->SetLineWidth(2);
        h_e_star_mu3[i]->SetLineWidth(2);
        h_q2[i]->GetYaxis()->SetTitleOffset(1.5);
        h_m2_miss[i]->GetYaxis()->SetTitleOffset(1.5);
        h_e_star_mu3[i]->GetYaxis()->SetTitleOffset(1.5);
    }

    // Declare vectors (for gof tests)
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

    cout<<"Filling loop starting"<<endl;

    for(int i=0; i<ne; ++i){
        t->GetEntry(i);

        auto it = find(v_classes.begin(), v_classes.end(), int(decayIndex)); // find class position in v_classes, v.end() if not found
        if(it == v_classes.end()) continue; // skip events not beloging to one of the classes to study

        int index = distance(v_classes.begin(), it); // iterator -> index

        h_q2[index]->Fill(q2);
        h_m2_miss[index]->Fill(m2_miss);
        h_e_star_mu3[index]->Fill(e_star_mu3);

        v_q2[index].push_back(q2);
        v_m2_miss[index].push_back(m2_miss);
        v_e_star_mu3[index].push_back(e_star_mu3);
    }

    // Get number of events for each class
    int v_ne[n_classes];
    for(int i=0; i<n_classes; ++i){
        v_ne[i] = h_q2[i]->GetEntries();
    }

    vector<vector<int>> v_groups;
    vector<int> v_lowStat;

    cout<<"Comparison loop starting"<<endl;

    // Compare distributions
    for(int i=0; i<n_classes; ++i){

        auto res = find(v_lowStat.begin(), v_lowStat.end(), v_classes[i]);
        if(res != v_lowStat.end()) continue;

        if(v_ne[i] < minEvt){
            v_lowStat.push_back(v_classes[i]);
            continue;
        }

        bool already_used = false;

        for(auto inner_vec:v_groups){
            if(already_used) break;
            auto res = find(begin(inner_vec), end(inner_vec), v_classes[i]);
            if(res != end(inner_vec))
                already_used = true;
        }

        if(already_used) continue;
        
        vector<int> v;
        v.push_back(v_classes[i]);

        for(int j=i+1; j<n_classes; ++j){

            auto res = find(v_lowStat.begin(), v_lowStat.end(), v_classes[j]);
            if(res != v_lowStat.end()) continue;

            if(v_ne[j] < minEvt){
                v_lowStat.push_back(v_classes[j]);
                continue;
            }

            bool already_used = false;

            for(auto inner_vec:v_groups){
                if(already_used) break;
                auto res = find(begin(inner_vec), end(inner_vec), v_classes[j]);
                if(res != end(inner_vec))
                    already_used = true;
            }

            if(already_used) continue;
            

            // Create got tests (syntax: nevt1, array1, nevt2, array2)
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

            double avg_pKS = (pKS_q2 + pKS_m2_miss + pKS_e_star_mu3)/3; // atm useless as the three distances/p-values are identical
            double avg_tKS = (tKS_q2 + tKS_m2_miss + tKS_e_star_mu3)/3; // atm useless as the three distances/p-values are identical

            bool good = avg_pKS >= pval && avg_tKS <= tval;

            if(good) v.push_back(v_classes[j]);

            delete gof_q2;
            delete gof_m2_miss;
            delete gof_e_star_mu3;
        }
        v_groups.push_back(v);
    }


    // Printing
    cout<<"--- GROUPS ("<<v_groups.size()<<") ---"<<endl;

    for(int i=0; i<v_groups.size(); ++i) {
        for(int j=0; j<v_groups[i].size(); ++j)
            cout<<v_groups[i][j]<<" ";
        cout<<endl;
    }

    cout<<"--- LOW STAT SAMPLES ---"<<endl;
    for(auto it:v_lowStat)
        cout<<it<<" ";
    cout<<endl;

    // Plotting
    auto c1 = new TCanvas("c1", "c1", 1600, 1200);
    c1->Divide(4,3);
    for(int i=0; i<v_groups.size(); ++i) {
        c1->cd(i+1);

        for(int j=0; j<v_groups[i].size(); ++j){
            int index = distance(v_classes.begin(), find(v_classes.begin(), v_classes.end(), v_groups[i][j]));

            h_q2[index]->SetLineColor(j+1>=10 ? j+2 : j+1); //skip color codes 0 and 10 (white)
            if(j==0) h_q2[index]->Draw("");
            else     h_q2[index]->Draw("SAME");
        }
    }
    c1->Print(dir + "/q2.png");

    auto c2 = new TCanvas("c2", "c2", 1600, 1200);
    c2->Divide(4,3);
    for(int i=0; i<v_groups.size(); ++i) {
        c2->cd(i+1);

        for(int j=0; j<v_groups[i].size(); ++j){
            int index = distance(v_classes.begin(), find(v_classes.begin(), v_classes.end(), v_groups[i][j]));

            h_m2_miss[index]->SetLineColor(j+1>=10 ? j+2 : j+1); //skip color codes 0 and 10 (white)
            if(j==0) h_m2_miss[index]->Draw("");
            else     h_m2_miss[index]->Draw("SAME");
        }
    }
    c2->Print(dir + "/m2_miss.png");

    auto c3 = new TCanvas("c3", "c3", 1600, 1200);
    c3->Divide(4,3);
    for(int i=0; i<v_groups.size(); ++i) {
        c3->cd(i+1);

        for(int j=0; j<v_groups[i].size(); ++j){
            int index = distance(v_classes.begin(), find(v_classes.begin(), v_classes.end(), v_groups[i][j]));

            h_e_star_mu3[index]->SetLineColor(j+1>=10 ? j+2 : j+1); //skip color codes 0 and 10 (white)
            if(j==0) h_e_star_mu3[index]->Draw("");
            else     h_e_star_mu3[index]->Draw("SAME");
        }
    }
    c3->Print(dir + "/e_star_mu3.png");

}
