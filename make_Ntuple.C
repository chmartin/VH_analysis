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

void make_Ntuple(std::string filename){
  
    ifstream fin;
    std::string filenameT = filename + ".lhe";
    std::cout << "Processing " << filenameT << std::endl;
    fin.open(filenameT.c_str());
    if(!fin.is_open())
	{
		std::cout << "Error!" << std::endl;
		fin.close();
		return 0;
	} 
	
    char oname[250];
    sprintf(oname,"%s.root",filename.c_str());
    std::cout << "Creating " << oname << std::endl;
    TFile fout(oname, "RECREATE");
    TTree* tree = new TTree("Events", "Events");
    
    
    TLorentzVector H_4vec(0.,0.,0.,0.);
    //float mH = 0.;
    TLorentzVector V0_4vec(0.,0.,0.,0.);
    int V0_id = 0;
    int V0_mid = 0;
    //float mV0 = 0.;
    TLorentzVector V1_4vec(0.,0.,0.,0.);
    int V1_id = 0;
    int V1_mid = 0;
    //float mV1 = 0.;
    TLorentzVector V2_4vec(0.,0.,0.,0.);
    int V2_id = 0;
    int V2_mid = 0;
    //float mV2 = 0.;
    std::vector <TLorentzVector> Daughter_4vec;
    std::vector <int> Daughter_id;
    std::vector <int> Daughter_mother;
    //std::vector <float> Daughter_mass;


    //Variables to be Stored
    tree->Branch("H", &H_4vec);
    //tree->Branch("mH", &mH);
    tree->Branch("V0", &V0_4vec);
    tree->Branch("V0_id", &V0_id);
    tree->Branch("V0_mother", &V0_mid);
    //tree->Branch("mV0", &mV0);
    tree->Branch("V1", &V1_4vec);
    tree->Branch("V1_id", &V1_id);
    tree->Branch("V1_mother", &V1_mid);
    //tree->Branch("mV1", &mV1);
    tree->Branch("V2", &V2_4vec);
    tree->Branch("V2_id", &V2_id);
    tree->Branch("V2_mother", &V2_mid);
    //tree->Branch("mV2", &mV2);
    tree->Branch("Daughters", &Daughter_4vec);
    tree->Branch("Daughter_id", &Daughter_id);
    tree->Branch("Daughter_mother", &Daughter_mother);
    //tree->Branch("Daughter_mass", &Daughter_mass);

    char buffer[500];
    int ctr = 0;
    while(!fin.eof() && ctr < 100)
	{
		fin.getline(buffer, 500);
		//std::cout << buffer << std::endl;
		if (strcmp(buffer, "<event>") == 0)
		{
				//std::cout << "GOT HERE" << std::endl;
				fin.getline(buffer, 255);
                                int numParticles = 0;
                                numParticles = atoi(strtok(buffer," \t"));
				int numPtot = numParticles;
				int HiggsLine = 99;
				int V0Line = 99;
				int V1Line = 99;
				int V2Line = 99;
				int foundV1 = 0;
				Daughter_4vec.clear();
				Daughter_id.clear();
				Daughter_mother.clear();
				//std::cout << numParticles << std::endl;
				fin >> buffer;
				while(strcmp(buffer, "</event>") != 0)
                                        {
						//std::cout << buffer << std::endl;
						numParticles--;
						int ParticleLine = numPtot - numParticles;
						//std::cout << "ParticleLine: " << ParticleLine << std::endl;
						//code to get 4vec, id, mother
						TLorentzVector temp_4vec(0.,0.,0.,0.);
						int ID = atoi(buffer);
                                                int IN_OUT = 0;
                                                int Mid0 = 0, Mid1 = 0;
						int Cid0 = 0, Cid1 = 0;
						float px = 0., py = 0., pz = 0., E = 0., M = 0., trash0, trash1;
						fin >> IN_OUT >> Mid0 >> Mid1 >> Cid0 >> Cid1 >> px >> py >> pz >> E >> M >> trash0 >> trash1;
						//std::cout << ID << " " << IN_OUT << " " << Mid0 << " " << Mid1 << " " << Cid0 << " " << Cid1 << " " << px << " " << py << " " << pz << " " << E << " " << M << " " << trash0 << " " << trash1 << std::endl;
						temp_4vec.SetPxPyPzE(px,py,pz,E);
						//std::cout << "MassDiff: " << (temp_4vec.M()-M) << std::endl;
						if (IN_OUT == 2 && ID == 25)
							{
								H_4vec = temp_4vec;
								HiggsLine = ParticleLine;
								//std::cout << "Found a Higgs at line: " << HiggsLine << std::endl;
								//temp_4vec.Print();
							}
						else if (IN_OUT == 2 && (Mid0 == 1 || Mid1 == 1) && (ID == 23 || ID == 24 || ID == -24))
							{
				                		V0_4vec = temp_4vec;
								V0_id = ID;
								V0_mid = 99;
								V0Line = ParticleLine;
                                                                //std::cout << "Found an associated V at line: " << V0Line << std::endl;
								//temp_4vec.Print();
							}
						else if (IN_OUT == 2 && (Mid0 == HiggsLine || Mid1 == HiggsLine) && (ID == 23 || ID == 24 || ID == -24) && foundV1 == 0)
							{
				                		V1_4vec = temp_4vec;
								V1_id = ID;
								V1_mid = 25;
								V1Line = ParticleLine;
								foundV1 = 1;
                                                                //std::cout << "Found a decay V at line: " << V1Line << std::endl;
								//temp_4vec.Print();
							}                
						else if (IN_OUT == 2 && (Mid0 == HiggsLine || Mid1 == HiggsLine) && (ID == 23 || ID == 24 || ID == -24) && foundV1 == 1)
							{
								V2_4vec = temp_4vec;
								V2_id = ID;
								V2_mid = 25;
								V2Line = ParticleLine;                      
                                                                //std::cout << "Found a decay V at line: " << V2Line << std::endl;
								//temp_4vec.Print();
							}
						else if (IN_OUT == 1)
							{
								Daughter_4vec.push_back(temp_4vec);
								Daughter_id.push_back(ID);
								if( Mid0 == V2Line || Mid0 == V1Line || Mid1 == V1Line || Mid1 == V2Line)
									{
										Daughter_mother.push_back(25);
									}
								else if (Mid0 == V0Line || Mid1 == V0Line)
									{
										Daughter_mother.push_back(V0_id);
									}
								else std::cout << "ERROR!" << std::endl;
								//std::cout << Mid0 << " " << Mid1 << std::endl;
								//std::cout << V0Line << " " << V1Line << " " << V2Line << std::endl;
								//std::cout << "DaugherID: " << ID << std::endl;
								//temp_4vec.Print();
							}
						fin >> buffer;
						//std::cout << foundV1 << std::endl;
						//if (numParticles == 0)
						//	{
						//		std::cout << buffer << std::endl;
						//		break;
						//	}
						
					}
				if (numParticles == 0)
					{
						tree->Fill();
						ctr++;
					}
				else
					{
						std::cout << "Event: " << ctr << " did not match expected number of particles." << std::endl;
					}
				fin.getline(buffer, 255);
			}
		//std::cout << buffer << std::endl;
	}
    cout << "Number of events processed: " << ctr << endl;
    char command[250];
    //sprintf(command,"rm %s.txt",filename.c_str());
    //gSystem->Exec(command);

    fout.cd();
    tree->Write();
    fout.Close();
	
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
