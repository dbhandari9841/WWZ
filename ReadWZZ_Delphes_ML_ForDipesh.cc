#include <TH1F.h> //a class in the ROOT data analysis framework used to represent a one-dimensional histogram with floating-point bin contents.
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>//a file that defines the TROOT class in ROOT

#include <TFile.h>
#include <TTree.h> //a file in ROOT that defines a TTree as a columnar dataset
#include <TSystem.h>	//base class defining a generic interface to the underlying Operating System
#include <TChain.h> // a file in ROOT that defines a TChain, which is a collection of files that contain TTree objects
#include <TLorentzVector.h> //a four-vector class in ROOT that can be used to describe position and time, or momentum and energy
#include <TLegend.h> 
#include <iostream> 
#include <stdlib.h>
#include <stdio.h> 
#include <algorithm> //provides a collection of functions that perform algorithmic operations on data structures
#include <TGraphAsymmErrors.h>
#include <TVector2.h> //header file from the ROOT framework, a widely used for data analysis 
#include <TF1.h> //hep data analysis specifically
#include <TRandom.h>
#include <iostream>
#include <fstream> //specifying name
using namespace std;
//cdddd
double MUON_MASS = 105.6583745*10e-03;

bool sameVal(double a, double b)
{
   return fabs(a - b) < 1.000e-02;   //boolean function
}

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double M)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiM(pT, eta, phi, M);
  return object_p4;
}

/*typedef struct
{
  double pT;
  double eta;
  double phi;
  int charge;
  TLorentzVector muonLV;
} MuonInfo;
*/

typedef struct
{
  double pT;
  double eta;
  double phi;
  int charge;
  TLorentzVector leptonLV;
} LeptonInfo;

typedef struct
{
  double pT;
  double eta;
  double phi;
} PhotonInfo;

typedef struct
{
  double pT;
  double eta;
  double phi;
  double mass;
  int btag;
} JetInfo;

typedef struct
{ 
  double pT;
  double eta;
  double phi;
  double energy;
  double mass;
  int    pdgID;
  int    status;
  int    motherone;
  int    mothertwo;
} GenInfo; 

typedef struct
{
  TLorentzVector JetLV;
  float BTag_CSV;
} AnalysisJetInfo; 

typedef struct
{
 double at;
 double size;
}pdgid;

bool sortLeptonsInDescendingpT(LeptonInfo lep1, LeptonInfo lep2)
{
  return (lep1.pT > lep2.pT);
}

bool sortJetsInDescendingpT(JetInfo jet1, JetInfo jet2)
{
  return (jet1.pT > jet2.pT);
}

bool sortGenPartInDescendingpT(GenInfo genpart1, GenInfo genpart2)
{
  return (genpart1.pT > genpart2.pT);
}

bool sortJetVectorsInDescendingpT(AnalysisJetInfo jet1, AnalysisJetInfo jet2)
{
  return (jet1.JetLV.Pt() > jet2.JetLV.Pt());
}
/*
bool sortJetVectorsInDescendingpT(TLorentzVector jet1, TLorentzVector jet2)
{
  return (jet1.Pt() > jet2.Pt());
}*/

bool sortPhotonsInDescendingpT(PhotonInfo pho1, PhotonInfo pho2)
{
  return (pho1.pT > pho2.pT);
}

double mT(TLorentzVector p4, float met, float met_phi)
{
  float phi1 = p4.Phi();
  float phi2 = met_phi;
  float Et1  = p4.Et();
  float Et2  = met;
  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}


//int jet_Z1_idx1 = -999;


/*void electmuonfilt()
{
    int z_lep1;
    int z_lep2;
    int z_lep3;
    int z_lep4;
    int w_lep;
    z_lep1=z_lep2=z_lep3=z_lep4=w_lep=-999;
    double Mz = 91.1876;
    double pair1massDiff, pair2massDiff;
    pair1massDiff=pair2massDiff=0.0;
    double compare1 = 15;
    double compare2 = 15;

    vector<int>	pdgid=	ana.tx.getBranchLazy<vector<int>	>("Common_lep_pdgid");
    
    for(unsigned int i=0; i<pdgid.size(); i++)
    {
        //std::cout << "pdgid.at(i) = " << pdgid.size() << std::endl;
        for(unsigned int j=0; j<pdgid.size(); j++)
        {
            if(i!=j)//make sure not checking same lepton
            {
               if((pdgid.at(i)*pdgid.at(j))==-11 or (pdgid.at(i)*pdgid.at(j))==11) //check opposite sign pair
                {
                    pair1massDiff = ((ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[i]+ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[j]).M() -  Mz);
                    for(unsigned int k=0; k<pdgid.size(); k++)
                    {
                        for(unsigned int l=0; l<pdgid.size(); l++)
                        {
                            if(j!=l and j!=k and i!=k and i!=l)//make sure not checking same lepton
                            {
                                if((pdgid.at(k)*pdgid.at(l))==-121 or (pdgid.at(k)*pdgid.at(l))==-169)
                                {
                                    pair2massDiff = ((ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[k]+ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[l]).M() -  Mz);
                                    if(fabs(pair1massDiff) < compare1 and fabs(pair2massDiff) < compare2)
                                    {
                                        compare1 = fabs(pair1massDiff);
                                        compare2 = fabs(pair2massDiff);
                                        z_lep1=i;
                                        z_lep2=j;
                                        z_lep3=k;
                                        z_lep4=l;

                                        for(unsigned int m=0; m<pdgid.size(); m++)
                                        {
                                          if(i!=m and j!=m and k!=m and l!=m)//make sure not checking same lepton
                                          {
                                             //mT = sqrt(2*ana.tx.getBranchLazy<LorentzVector>("Common_met_p4").pt()*ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[m].Pt()*(1.0 - cos(ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[m].Phi() - ana.tx.getBranchLazy<LorentzVector>("Common_met_p4").phi())));
                                             w_lep=m;
                                          }//w-tag
                                        }//m-loop
                                     }//zz-tag
                                  }//second opposite sign pair
                               }//lepton index check
                            }//l-loop
                         }//k-loop
                      }//first opposite sign pair
                    }//lepton index check
                 }//j-loop
              }//i-loop

    lep_5Lep_Z1_idx1 = z_lep1;
    lep_5Lep_Z1_idx2 = z_lep2;
    lep_5Lep_Z2_idx1 = z_lep3;
    lep_5Lep_Z2_idx2 = z_lep4;
    lep_5Lep_W_idx = w_lep;
    //if(lep_5Lep_W_idx>0 and lep_5Lep_Z1_idx1>0 and lep_5Lep_Z1_idx2>0 and lep_5Lep_Z2_idx1>0 and lep_5Lep_Z2_idx2>0) std::cout << lep_5Lep_W_idx << std::endl;

    return;
}
*/

void select5leptons(vector<LeptonInfo> lepton_vector, vector<LeptonInfo> electron_vector, vector<LeptonInfo> muon_vector)
{
    int z1_lepton1;
    int z1_lepton2;
    int z2_lepton1;
    int z2_lepton2;
    int w_lepton;
    z1_lepton1=z1_lepton2=z2_lepton1=z2_lepton2=w_lepton=-999;
    double Mz = 91.1876;
    double pair1massDiff, pair2massDiff;
    pair1massDiff=pair2massDiff=0.0;
    double compare1 = 15;
    double compare2 = 15; 

    //vector<int>	pdgid=	ana.tx.getBranchLazy<vector<int>	>("Common_lep_pdgid");

    for(unsigned int a=0; a<lepton_vector.size(); a++)
    {
	  for(unsigned int b=0; b<lepton_vector.size(); b++)
          {
		 if (a!=b)//making sure not checking the same lepton
	         {
	              if((electron_vector.at(a).charge*electron_vector.at(b).charge)==-1  || (muon_vector.at(a).charge*muon_vector.at(b).charge)==-1)//check opposite sign pair
			{
				pair1massDiff=((lepton_vector.at(a).leptonLV+lepton_vector.at(b).leptonLV).M() -  Mz);
			    for(unsigned int c=0; c<lepton_vector.size(); c++)
			    { 
				   for(unsigned int d=0; d<lepton_vector.size(); d++)
				   {
						if(b!=d and b!=c and a!=c and a!=d)//making sure not checking same lepton
							{
								if ((electron_vector.at(c).charge*electron_vector.at(d).charge)==-1 || (muon_vector.at(c).charge*muon_vector.at(d).charge)==-1) //checking opposite sign pair 
										pair2massDiff=((lepton_vector.at(c).leptonLV+lepton_vector.at(d).leptonLV).M() -  Mz);
						if (fabs(pair1massDiff)< compare1 and fabs(pair2massDiff)< compare2)
						{
							compare1 = fabs(pair1massDiff);
							compare2 = fabs(pair2massDiff);
							z1_lepton1=a;
							z1_lepton2=b;
							z2_lepton1=c;
							z2_lepton2=d;

							for (unsigned int e=0; e<lepton_vector.size(); e++)
							{
								if(a!=e and b!=e and c!=e and d!=e)//make sure not checking same lepton 
							{//mT= sqrt(2*met*muon_vector.at(e).pT*(1.0 - cos(muon_vector.at(e).phi - muon_vector.phi)));
								w_lepton=e;
							 }//wtag
				}//m-loop
						       }//zz-tag
                                  }//second opposite sign pair
                               }//muon index check
                            }//l-loop
                         }//k-loop
                      }//first opposite sign pair
                    }//muon index check
                 }//j-loop
              //i-loop

       int lep_5lepton_Z1_idx1 = z1_lepton1;
       int lep_5lepton_Z1_idx2 = z1_lepton2;
       int lep_5lepton_Z2_idx1 = z2_lepton1;
       int lep_5lepton_Z2_idx2 = z2_lepton2;
       int lep_5lepton_W_idx = w_lepton;
	
	return;  
} 
/*
void select5muons(vector<MuonInfo> muon_vector)
{
    int z1_muon1; 
    int z1_muon2;
    int z2_muon1;
    int z2_muon2;
    int w_muon;
    z1_muon1=z1_muon2=z2_muon1=z2_muon2=w_muon=-999;
    double Mz=91.1876;
    double pair1massdiff=0.0;
    double pair2massdiff=0.0;
    double compare1 = 15;
    double compare2 = 15; 

    for(unsigned int a=0; a<muon_vector.size(); a++)
    {
	  for(unsigned int b=0; b<muon_vector.size(); b++)
          {
		 if (a!=b)//making sure not checking the same muon
	         {
	              if((muon_vector.at(a).charge*muon_vector.at(b).charge)==-1)//check opposite sign pair
			{
				pair1massdiff=((muon_vector.at(a).muonLV+muon_vector.at(b).muonLV).M() -  Mz);
			    for(unsigned int c=0; c<muon_vector.size(); c++)
			    { 
				   for(unsigned int d=0; d<muon_vector.size(); d++)
				   {
						if(b!=d and b!=c and a!=c and a!=d)//making sure not checking same meson
							{
								if (muon_vector.at(c).charge*muon_vector.at(d).charge==-1) //checking opposite sign pair 
										pair2massdiff=((muon_vector.at(c).muonLV+muon_vector.at(d).muonLV).M() -  Mz);
						if (fabs(pair1massdiff)< compare1 and fabs(pair2massdiff)< compare2)
						{
							compare1 = fabs(pair1massdiff);
							compare2 = fabs(pair2massdiff);
							z1_muon1=a;
							z1_muon2=b;
							z2_muon1=c;
							z2_muon2=d;

							for (unsigned int e=0; e<muon_vector.size(); e++)
							{
								if(a!=e and b!=e and c!=e and d!=e)//make sure not checking same muon 
							{//mT= sqrt(2*met*muon_vector.at(e).pT*(1.0 - cos(muon_vector.at(e).phi - muon_vector.phi)));
								w_muon=e;
							 }//wtag
				}//m-loop
						       }//zz-tag
                                  }//second opposite sign pair
                               }//muon index check
                            }//l-loop
                         }//k-loop
                      }//first opposite sign pair
                    }//muon index check
                 }//j-loop
              //i-loop

       int lep_5muon_Z1_idx1 = z1_muon1;
       int lep_5muon_Z1_idx2 = z1_muon2;
       int lep_5muon_Z2_idx1 = z2_muon1;
       int lep_5muon_Z2_idx2 = z2_muon2;
       int lep_5muon_W_idx = w_muon;
	
	return;
}
    */	
void select4Jets(vector<AnalysisJetInfo> Jet_vector)
{
    int z_jet1;
    int z_jet2;
    int z_jet3;
    int z_jet4;
    z_jet1=z_jet2=z_jet3=z_jet4=-999;
    double Mz = 91.1876;
    double pair1massDiff, pair2massDiff;
    pair1massDiff=pair2massDiff=0.0;
    double compare1 = 30;
    double compare2 = 30;

    for(unsigned int i=0; i<Jet_vector.size(); i++)
    {
      for(unsigned int j=0; j<Jet_vector.size(); j++)
      {
        if(i!=j)//make sure not checking same jet
        {
          pair1massDiff = ((Jet_vector.at(i).JetLV+Jet_vector.at(j).JetLV).M() -  Mz);
          for(unsigned int k=0; k<Jet_vector.size(); k++)
          {
            for(unsigned int l=0; l<Jet_vector.size(); l++)
            {
              if(j!=l and j!=k and i!=k and i!=l and k!=l)//make sure not checking same jet
              {
                pair2massDiff = ((Jet_vector.at(k).JetLV+Jet_vector.at(l).JetLV).M() -  Mz);
                if(fabs(pair1massDiff) < compare1 and fabs(pair2massDiff) < compare2)
		{
                  compare1 = fabs(pair1massDiff);
                  compare2 = fabs(pair2massDiff);
                  z_jet1=i;
                  z_jet2=j;
                  z_jet3=k;
                  z_jet4=l;
                  //z_muon=	
                }//zz-tag
              }//lepton index check
            }//l-loop
          }//k-loop
        }//lepton index check
      }//j-loop
    }//i-loop

   int jet_Z1_idx1 = z_jet1;
   int jet_Z1_idx2 = z_jet2;
   int jet_Z2_idx1 = z_jet3;
   int jet_Z2_idx2 = z_jet4;
    //if(jet_W_idx>0 and jet_Z1_idx1>0 and jet_Z1_idx2>0 and jet_Z2_idx1>0 and jet_Z2_idx2>0) std::cout << jet_W_idx << std::endl;

    return;
}


int ReadWZZ_Delphes_ML_ForDipesh(std::string infile, std::string outfile){

  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("LowPtSUSY_Tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  std::vector<double>   *ph_pt;
  std::vector<double>   *ph_phi;
  vector<double>   *ph_eta;
  Int_t           nPhotons;
  vector<double>   *el_pt;
  vector<double>   *el_phi;
  vector<double>   *el_eta;
  vector<int>     *el_charge;
  Int_t           nElectrons;
  vector<double>   *mu_pt;
  vector<double>   *mu_phi;
  vector<double>   *mu_eta;
  vector<int>     *mu_charge;
  Int_t           nMuons;
  vector<double>   *jet_pt;
  vector<double>   *jet_phi;
  vector<double>   *jet_eta;
  vector<double>   *jet_mass;
  vector<int>     *jet_btag;
  Int_t           nJets;
  Float_t         MET;
  Float_t         MET_Phi;
  vector<int>  *GenParticle_PDGId;
  vector<double>  *GenParticle_Pt;
  vector<double>  *GenParticle_Phi;
  vector<double>  *GenParticle_Eta;
  vector<double>  *GenParticle_Mass;
  vector<double>  *GenParticle_Energy;
  vector<int>  *GenParticle_Status;
  vector<int>  *GenParticle_MotherOne;
  vector<int>  *GenParticle_MotherTwo;

  ph_pt = 0;
  ph_phi = 0;
  ph_eta = 0;
  el_pt = 0;
  el_phi = 0;
  el_eta = 0;
  el_charge = 0;
  mu_pt = 0;
  mu_phi = 0;
  mu_eta = 0;
  mu_charge = 0;
  jet_pt = 0;
  jet_phi = 0;
  jet_eta = 0;
  jet_mass = 0;
  jet_btag = 0;
  GenParticle_PDGId = 0;
  GenParticle_Status = 0;
  GenParticle_Pt = 0;
  GenParticle_Phi = 0;
  GenParticle_Eta = 0;
  GenParticle_Mass = 0;
  GenParticle_Energy = 0;
  GenParticle_MotherOne = 0;
  GenParticle_MotherTwo = 0;

  tree->SetBranchAddress("ph_pt", &(ph_pt));
  tree->SetBranchAddress("ph_phi", &(ph_phi));
  tree->SetBranchAddress("ph_eta", &(ph_eta));
  tree->SetBranchAddress("nPhotons", &(nPhotons));
  tree->SetBranchAddress("el_pt", &(el_pt));
  tree->SetBranchAddress("el_eta", &(el_eta));
  tree->SetBranchAddress("el_phi", &(el_phi));
  tree->SetBranchAddress("el_charge", &(el_charge));
  tree->SetBranchAddress("nElectrons", &(nElectrons));
  tree->SetBranchAddress("mu_pt", &(mu_pt));
  tree->SetBranchAddress("mu_eta", &(mu_eta));
  tree->SetBranchAddress("mu_phi", &(mu_phi));
  tree->SetBranchAddress("mu_charge", &(mu_charge));
  tree->SetBranchAddress("nMuons", &(nMuons));
  tree->SetBranchAddress("jet_pt", &(jet_pt));
  tree->SetBranchAddress("jet_phi", &(jet_phi));
  tree->SetBranchAddress("jet_eta", &(jet_eta));
  tree->SetBranchAddress("jet_mass", &(jet_mass));
  tree->SetBranchAddress("jet_btag", &(jet_btag));
  tree->SetBranchAddress("nJets", &(nJets));
  tree->SetBranchAddress("MET", &(MET));
  tree->SetBranchAddress("MET_Phi", &(MET_Phi));
  tree->SetBranchAddress("GenParticle_PDGId", &(GenParticle_PDGId));
  tree->SetBranchAddress("GenParticle_Pt", &(GenParticle_Pt));
  tree->SetBranchAddress("GenParticle_Phi", &(GenParticle_Phi));
  tree->SetBranchAddress("GenParticle_Eta", &(GenParticle_Eta));
  tree->SetBranchAddress("GenParticle_Mass", &(GenParticle_Mass));
  tree->SetBranchAddress("GenParticle_Energy", &(GenParticle_Energy));
  tree->SetBranchAddress("GenParticle_Status", &(GenParticle_Status));
  tree->SetBranchAddress("GenParticle_MotherOne", &(GenParticle_MotherOne));
  tree->SetBranchAddress("GenParticle_MotherTwo", &(GenParticle_MotherTwo));

  TH1F *h_ST = new TH1F("h_ST", "ST (scalar sum of jet + lepton pT + MET); S_T [GeV]; Events/GeV", 100, 0, 5000.0);h_ST->Sumw2();
  TH1F *h_VVV_Mass = new TH1F("h_VVV_Mass", "h_VVV_Mass; Invariant Mass of VVV [GeV]; Events", 100, 0, 5000.0); h_VVV_Mass->Sumw2();
  TH1F *h_nleptons = new TH1F("h_nleptons", "h_leptons; Number of Leptons; Events", 10, -0.5, 9.5); h_nleptons->Sumw2();
  TH1F *h_nmuons = new TH1F("h_nmuons", "h_muons; Number of Muons; Events", 10, -0.5, 9.5); h_nmuons->Sumw2(); 
  //here

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  ofstream input_csv;
  input_csv.open ("input_csv.txt");
  ofstream input_csv_withJets;
  input_csv_withJets.open ("input_csv_withJets.txt");

  //lepton vector, sum up electron and muon vectors...
  for (int i=0; i<nEvents ; ++i)
  {
    tree->GetEvent(i);
    //vector that stores the electron properties
    std::vector<LeptonInfo> electrons;
    for (unsigned int j=0; j<el_pt->size(); ++j)
    {
      LeptonInfo electron;
      electron.pT=el_pt->at(j);
      electron.eta=el_eta->at(j);
      electron.phi=el_phi->at(j);
      electron.charge=el_charge->at(j);
      electron.leptonLV.SetPtEtaPhiM(electron.pT, electron.eta, electron.phi, 0.000511);
      electrons.push_back(electron);
    }
    // Now sorting this vector of structs
    std::sort (electrons.begin(), electrons.end(), sortLeptonsInDescendingpT);
   
    // filling the muon's properties into a vector of struct
    std::vector<LeptonInfo> muons;
    for (unsigned int j=0; j<mu_pt->size(); ++j)
    {
      LeptonInfo muon;
      muon.pT=mu_pt->at(j);
      muon.eta=mu_eta->at(j);
      muon.phi=mu_phi->at(j);
      muon.charge=mu_charge->at(j);
      muon.leptonLV.SetPtEtaPhiM(muon.pT, muon.eta, muon.phi, 0.1056);//replace with correct muon mass
      muons.push_back(muon);
    }
    // sorting this vector of structs
    std::sort (muons.begin(), muons.end(), sortLeptonsInDescendingpT); 
    
    //combining electrons and muons into a single vector
    std::vector<LeptonInfo> lepton_vector = electrons; //start with e-
    lepton_vector.insert(lepton_vector.end(), muons.begin(), muons.end()); //appending muons to the end of the electron vector

    //sorting the combined lepton vector by descending pT
    std::sort(lepton_vector.begin(), lepton_vector.end(), sortLeptonsInDescendingpT);
    

    int n_leptons = electrons.size() + muons.size(); h_nleptons-> Fill(n_leptons); //HERE
    //int n_electrons
    //int n_muons = muons.size(); h_nmuons-> Fill(n_muons);//HERE
    
    std::vector<JetInfo> jets;
    for (unsigned int j=0; j<jet_pt->size(); ++j)
    {
      JetInfo jet;
      jet.pT = jet_pt->at(j);
      jet.eta = jet_eta->at(j);
      jet.phi = jet_phi->at(j);
      jet.mass = jet_mass->at(j);
      jet.btag = jet_btag->at(j);
      jets.push_back(jet);
    }

    // Now sorting this vector of structs
    std::sort (jets.begin(), jets.end(), sortJetsInDescendingpT);   
  
    std::vector<GenInfo> genparts;
    std::vector<GenInfo> genmuons;
    std::vector<GenInfo> genpartsZ;
    std::vector<GenInfo> genneutrinos;
    for (unsigned int j=0; j<GenParticle_PDGId->size(); ++j)
    {
      GenInfo genpart;
      genpart.pT = GenParticle_Pt->at(j);
      genpart.eta = GenParticle_Eta->at(j);
      genpart.phi = GenParticle_Phi->at(j);
      genpart.energy = GenParticle_Energy->at(j);
      genpart.mass = GenParticle_Mass->at(j);
      genpart.pdgID = GenParticle_PDGId->at(j);
      genpart.status = GenParticle_Status->at(j);
      genpart.motherone = GenParticle_MotherOne->at(j);
      genpart.mothertwo = GenParticle_MotherTwo->at(j);
      if(abs(genpart.pdgID)==24 and abs(genpart.status)>20 and abs(genpart.status)<30) genparts.push_back(genpart);
      if(abs(genpart.pdgID)==23 and abs(genpart.status)>20 and abs(genpart.status)<30) genpartsZ.push_back(genpart);
      if(abs(genpart.pdgID)==13 and abs(genpart.status)>20 and abs(genpart.status)<30) genmuons.push_back(genpart);
      if(abs(genpart.pdgID)==14 and abs(genpart.status)>20 and abs(genpart.status)<30) genneutrinos.push_back(genpart);
      //for pythia-6
      //if(abs(genpart.pdgID)==24 and abs(genpart.status)==3) genparts.push_back(genpart);
    }

    std::sort (genparts.begin(), genparts.end(), sortGenPartInDescendingpT);
    std::sort (genmuons.begin(), genmuons.end(), sortGenPartInDescendingpT);
    std::sort (genpartsZ.begin(), genpartsZ.end(), sortGenPartInDescendingpT);
     
    TLorentzVector W1, Z1, Z2;
    if(genparts.size()  > 0.0) W1.SetPtEtaPhiE(genparts.at(0).pT, genparts.at(0).eta, genparts.at(0).phi, genparts.at(0).energy); 
    if(genpartsZ.size() > 0.0) Z1.SetPtEtaPhiE(genpartsZ.at(0).pT, genpartsZ.at(0).eta, genpartsZ.at(0).phi, genpartsZ.at(0).energy);
    if(genpartsZ.size() > 1.0) Z2.SetPtEtaPhiE(genpartsZ.at(1).pT, genpartsZ.at(1).eta, genpartsZ.at(1).phi, genpartsZ.at(1).energy);
    h_VVV_Mass->Fill((W1+Z1+Z2).M());

    double ST = 0.0; //jets are sorted. Don't care as far as ST is concerned.
    double HTb = 0.0; 
    vector<AnalysisJetInfo> Jet_vector;
    Jet_vector.clear();
    vector<AnalysisJetInfo> bJet_vector;
    bJet_vector.clear();
    for(unsigned int k=0; k<jets.size(); ++k)
    {
      AnalysisJetInfo Jet;
      AnalysisJetInfo bJet;
      if(fabs(jets.at(k).eta)<2.4 and jets.at(k).pT>30.0)
      {
        Jet.JetLV.SetPtEtaPhiM(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).mass);
        Jet.BTag_CSV = jets.at(k).btag;
        bool isGoodJet=true;
        for(unsigned int j=0; j<electrons.size(); ++j)
        {
          TLorentzVector Electron;
          Electron.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
          Electron.SetPtEtaPhiM(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, 0.0);
          double DRjet_el = Jet.JetLV.DeltaR(Electron);
          if(DRjet_el<0.5) isGoodJet=false;
        }
        for(unsigned int j=0; j<muons.size(); ++j)
        {
          TLorentzVector Muon;
          Muon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
          Muon.SetPtEtaPhiM(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, 0.0);
          double DRjet_mu = Jet.JetLV.DeltaR(Muon);
          if(DRjet_mu<0.5) isGoodJet=false;
        }  
        if(isGoodJet) Jet_vector.push_back(Jet);
      } 
      if(fabs(jets.at(k).eta)<2.4 and jets.at(k).pT>30.0 and jets.at(k).btag>0)
      {
        bJet.JetLV.SetPtEtaPhiM(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).mass);
        bJet.BTag_CSV = jets.at(k).btag;
        bool isGoodbJet=true;
        for(unsigned int j=0; j<electrons.size(); ++j)
        {
          TLorentzVector Electron;
          Electron.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
          Electron.SetPtEtaPhiM(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, 0.0);
          double DRjet_el = bJet.JetLV.DeltaR(Electron);
          if(DRjet_el<0.5) isGoodbJet=false;
        }
        for(unsigned int j=0; j<muons.size(); ++j)
        {
          TLorentzVector Muon;
          Muon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
          Muon.SetPtEtaPhiM(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, 0.0);
          double DRjet_mu = bJet.JetLV.DeltaR(Muon);
          if(DRjet_mu<0.5) isGoodbJet=false;
        } 
        if(isGoodbJet) bJet_vector.push_back(bJet);
      }//close four vector if
    }//close jet loop

    // Now sorting this vector of structs
    std::sort (Jet_vector.begin(), Jet_vector.end(), sortJetVectorsInDescendingpT);

    for(unsigned int m=0; m<Jet_vector.size(); m++)
    {
      ST += Jet_vector.at(m).JetLV.Pt();
    }

    for(unsigned int imu=0; imu<muons.size(); imu++)
    {
      ST += muons.at(imu).pT;
    }

    for(unsigned int iel=0; iel<electrons.size(); iel++)
    {
      ST += electrons.at(iel).pT;
    }
  
    ST += MET; 

    h_ST->Fill(ST);

    double inv_mass = (W1+Z1+Z2).M();
    //for (unsigned int j=0; j<jet_pt->size(); ++j)
    if (std::abs(pdgid.at(i)) == 11){
    input_csv << (W1+Z1+Z2).M() << "," << ST << electrons.at(i).pT << "," << electrons.at(i).eta << "," << electrons.at(i).phi << "," << 0 << "," << 0 << "," << 0 << std::endl;//for electron}
    
    else if(std::abs(pdgid.at(i)) == 13){ 
    input_csv << (W1+Z1+Z2).M() << "," << ST << 0 << "," << 0 << "," << 0 << "," << muons.at(i).pT << "," << muons.at(i).eta << "," << muons.at(i).phi << std::endl;//for muon}
    }
    //pT, azimuthal angle, eta
  }//event loop closed
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_VVV_Mass->Write();
  //h_ST->Write();  
  h_nleptons->Write();
  //h_nmuons->Write();//here
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0;
}
