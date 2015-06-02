#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include "TFile.h"
#include "TList.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TBranch.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

using namespace std;

TRandom randomForSmearing;

TLorentzVector applyResolution(TLorentzVector l_gen);

void read_Ntuple(std::string filename, bool applyRes=false){
  
    std::string filenameT = filename + ".root";
    std::cout << "Processing " << filenameT << std::endl;
    TFile* inFile = new TFile(filenameT.c_str(),"OPEN");
    TTree* tree = inFile->Get("Events");
    
    
    TLorentzVector *H_4vec = new TLorentzVector(0.,0.,0.,0.);
    //float mH = 0.;
    TLorentzVector *V0_4vec = new TLorentzVector(0.,0.,0.,0.);
    int V0_id = 0;
    int V0_mid = 0;
    //float mV0 = 0.;
    TLorentzVector *V1_4vec = new TLorentzVector(0.,0.,0.,0.);
    int V1_id = 0;
    int V1_mid = 0;
    //float mV1 = 0.;
    TLorentzVector *V2_4vec = new TLorentzVector(0.,0.,0.,0.);
    int V2_id = 0;
    int V2_mid = 0;
    //float mV2 = 0.;
    std::vector <TLorentzVector> Daughter_4vec;
    std::vector <int> Daughter_id;
    std::vector <int> Daughter_mother;
    //std::vector <float> Daughter_mass;


    //Variables to be Stored
    tree->SetBranchAddress("H", &H_4vec);
    //tree->Branch("mH", &mH);
    tree->SetBranchAddress("V0", &V0_4vec);
    tree->SetBranchAddress("V0_id", V0_id);
    tree->SetBranchAddress("V0_mid", V0_mid);
    //tree->Branch("mV0", &mV0);
    tree->SetBranchAddress("V1", &V1_4vec);
    tree->SetBranchAddress("V1_id", V1_id);
    tree->SetBranchAddress("V1_mid", V1_mid);
    //tree->Branch("mV1", &mV1);
    tree->SetBranchAddress("V2", &V2_4vec);
    tree->SetBranchAddress("V2_id", &V2_id);
    tree->SetBranchAddress("V2_mid", V2_mid);
    //tree->Branch("mV2", &mV2);
    tree->SetBranchAddress("Daughters", &Daughter_4vec);
    tree->SetBranchAddress("Daughter_id", &Daughter_id);
    tree->SetBranchAddress("Daughter_mother", &Daughter_mother);
    //tree->Branch("Daughter_mass", &Daughter_mass);

    int ctr = 0;


    for(int event = 0; event < tree->GetEntries(); event++)
	{
		tree->GetEntry(event);
		std::cout << H_4vec.M() << std::endl;
	}

    std::cout << "Number of events processed: " << ctr << std::endl;
    char command[250];
    //sprintf(command,"rm %s.txt",filename.c_str());
    //gSystem->Exec(command);

    //fout.cd();
    //tree->Write();
    inFile.Close();
	
}

TLorentzVector applyResolution(TLorentzVector l_gen){

  float l_Perp, l_Theta, l_Phi;

  if(randomForSmearing.Uniform()<.9){
    l_Perp = l_gen.Perp()+(randomForSmearing.Gaus(0,0.012*1.15*l_gen.Perp()+0.00000*1.15* l_gen.Perp()* l_gen.Perp() ));
    l_Theta = l_gen.Theta()+(randomForSmearing.Gaus(0,0.001));
    l_Phi = l_gen.Phi()+(randomForSmearing.Gaus(0,0.001));
  }else{
    l_Perp = l_gen.Perp()+randomForSmearing.Gaus(-l_gen.Perp()*.04,l_gen.Perp()*.08);
    l_Theta = l_gen.Theta()+(randomForSmearing.Gaus(0,0.001));
    l_Phi = l_gen.Phi()+(randomForSmearing.Gaus(0,0.001));
  }
  float l_Px = l_Perp*cos(l_Phi);
  float l_Py = l_Perp*sin(l_Phi);
  float l_Pz = l_Perp/tan(l_Theta);
  float l_E = sqrt(l_Px*l_Px+l_Py*l_Py+l_Pz*l_Pz);

  TLorentzVector final_l;

  final_l.SetPxPyPzE(l_Px,l_Py,l_Pz,l_E);

  return (final_l);

}
