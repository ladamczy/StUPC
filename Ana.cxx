// Run by: ./Ana file.list
// e.g. ./Ana /star/u/truhlar/star-upcDst/build/run17.list
// or you can open just root file
// ./Ana /star/data01/pwg_tasks/upc02/UPC_P20ic/18061097.root
// or you can open n-th root file in file.list
// ./Ana file.list index 


// Table of RP indecies and names
// RP_ID   0,    1,    2,   3,   4,   5,   6, 7
// RP_name E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D

// c++ headers
#include <iostream>
#include <string>    
#include <utility>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cstdlib>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include <TH2.h> 
#include <TF1.h> 
#include <TF2.h> 
#include <THStack.h> 
#include <TStyle.h> 
#include <TGraph.h> 
#include <TGraph2D.h> 
#include <TGraphErrors.h> 
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TGaxis.h> 
#include <TString.h> 
#include <TColor.h> 
#include <TLine.h> 
#include <TExec.h> 
#include <TFitResultPtr.h> 
#include <TFitResult.h> 
#include <TLatex.h> 
#include <TMath.h>
#include <TLorentzVector.h>

// picoDst headers
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"

using namespace std;

// enums are very usefull 
enum { kAll = 1, kCPT,  kRP, kOneVertex, kTPCTOF, 
    kTotQ, kMax};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

const double particleMass[nParticles] = { 0.13957, 0.49367, 0.93827}; // pion, kaon, proton in GeV /c^2 
const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const int CEPtriggers[] = { 570701, 570705, 570711, 590701, 590705, 590708};

TFile *infile, *outfile;
TChain *upcChain, *rpChain;
StUPCEvent *upcEvt;
StRPEvent *rpEvt;
TTree *upcTree, *recTree, *bcgTree;


TH1I *hTriggerBits, *hAnalysisFlow; // control plots
TH2F *hEtaPhi;
TH1D *hInvMass, *hPt;

TH1F* hNumberTrackPerBranch[nBranches]; // number of track per each brunch (west up, west down, east up, east down)
TH1F* hNumberTrack;

UInt_t runNumber;
Double_t invMass;
Double_t vertexeZ;
 


void Make();
void Init();
bool ConnectInput(int argc, char** argv);
TFile *CreateOutputTree(const string& out);
bool CheckTriggers();

//_____________________________________________________________________________
int main(int argc, char** argv) 
{
    //connect input file(s)
    if(!ConnectInput(argc, argv))
    {
        cout << "Wrong input parameters..." << endl; 
        return 1;
    }
    //open output file
    outfile = CreateOutputTree("AnaOutput.root"); 
    if(!outfile) 
    {
        cout << "Can not open output file." << endl; 
        return -1;
    }

    // Initiate histograms
    Init();

    //ask for number of events
    Long64_t nev = upcTree->GetEntries();
    cout<<"Proccesing "<<nev<<" events"<<endl;

    //event loop
    //nev = 10000; // use for debugging
    for(Long64_t iev=0; iev<nev; ++iev) 
    { //get the event
        upcTree->GetEntry(iev); 
        Make();
    } 

    //close the outputs
    outfile->Write();
    outfile->Close();
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main

void Make()
{
    hAnalysisFlow->Fill(kAll);
    if(!CheckTriggers())
        return;

    hAnalysisFlow->Fill(kCPT);

    // Vectors below will be filled with indices of good-quality tracks
    vector<int> rpTrackIdVec_perBranch[nBranches];
    vector<int> rpTrackIdVec_perSide[nSides];
    int numberOfTracks = 0;
    int numberOfTracksPerBranch[nBranches] = {0, 0, 0, 0};

    // Loop over all tracks reconstructed in Roman Pots  
    for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
    {
        // Get pointer to k-th track in Roman Pot data collection
        StUPCRpsTrack *trk = rpEvt->getTrack(k);
        trk->setEvent(rpEvt);
        // Get ID of a branch in which this k-th track was reconstructed
        int j = trk->branch();
        int side = j<2 ? E : W;

        ++numberOfTracks;
        ++numberOfTracksPerBranch[j];

        // Analyse proton tracks with at least 3 planes (per RP) used in the track reconstraction
        if( (trk->getTrackPoint(0) ? trk->getTrackPoint(0)->planesUsed()>=3 : false) &&
            (trk->getTrackPoint(1) ? trk->getTrackPoint(1)->planesUsed()>=3 : false))
        {
            rpTrackIdVec_perBranch[j].push_back( k );
            rpTrackIdVec_perSide[side].push_back( k );
        } 
    }

    // Fill histograms with number of reconstructed RP tracks 
    hNumberTrack->Fill(numberOfTracks);
    for(int i=0; i<nBranches; ++i) 
        hNumberTrackPerBranch[i]->Fill(numberOfTracksPerBranch[i]);

    // Alow only 1 proton on each side
    if( !(rpTrackIdVec_perSide[E].size()==1 && rpTrackIdVec_perSide[W].size()==1))
        return;

    hAnalysisFlow->Fill(kRP);

    // Skip all events with more vertecies than 1
    if( upcEvt->getNumberOfVertices() != 1) return;
    hAnalysisFlow->Fill(kOneVertex);

    // Get z-position of vertex
    vertexeZ = upcEvt->getVertex(0)->getPosZ();
    // get run number
    runNumber = upcEvt->getRunNumber();

    int totalCharge = 0;
    int nTOFtracks = 0;
    double eta, phi, pT;

    TLorentzVector state;
    for(int trackID = 0; trackID < upcEvt->getNumberOfTracks(); ++trackID)
    {
        const StUPCTrack* trk = upcEvt->getTrack(trackID);
        // Skip all tracks that are not primary or they are not matched with TOF
        if( !trk->getFlag(StUPCTrack::kPrimary) || !trk->getFlag(StUPCTrack::kTof)) 
            continue;
        nTOFtracks++;

        phi = trk->getPhi();
        eta = trk->getEta();  
        pT = trk->getPt();
        hEtaPhi->Fill(eta, phi);
        hPt->Fill(pT);

        totalCharge += static_cast<int>( trk->getCharge() );

        TLorentzVector track;
        trk->getLorentzVector(track, particleMass[Pion]);
        state += track;
    } 

    // Skip events with less than 2 TOF tracks
    if(nTOFtracks < 2)
        return;

    hAnalysisFlow->Fill(kTPCTOF);

    invMass = state.M();
    hInvMass->Fill(invMass);

    if(totalCharge) // Total charge is not zero => background event
    {
        bcgTree->Fill(); // Fill background Tree
    }else{
        hAnalysisFlow->Fill(kTotQ);
        recTree->Fill(); // Fill analysis (reco) Tree
    }

}

void Init()
{

    hAnalysisFlow = new TH1I("AnalysisFlow", "CutsFlow", kMax-1, 1, kMax);
    TString summaryLabels[] = { TString("All"), TString("CPT"), TString("RP"),
        TString("1 vertex"), TString("2+ TOF"), TString("TotQ=0")};
    for(int tb=1; tb<kMax; ++tb) 
        hAnalysisFlow->GetXaxis()->SetBinLabel(tb, summaryLabels[tb-1]);

    hTriggerBits = new TH1I("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
    for(int tb=0; tb<nTriggers; ++tb){
        TString label; label.Form("%d",triggerID[tb]);
        hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
    }

    TString branchName[nBranches] = { TString("EU"), TString("ED"), TString("WU"), TString("WD") };
    for(int i=0; i<nBranches; ++i)
        hNumberTrackPerBranch[i] = new TH1F("NumberTracksPerBranch_"+branchName[i],"Number of tracks in branch "+branchName[i], 8, -0.5, 7.5);
    hNumberTrack = new TH1F("NumberTracks", "Number of Tracks in RP", 40, -0.5, 39.5);
    
    
    hEtaPhi = new TH2F("hEtaPhi","Phi vs eta of TOF tracks",100,-2,2,100,-4,4);
    hPt = new TH1D("hPt", "Transverse momentum", 100, 0, 10);
    hInvMass = new TH1D("hInvMass", "Uncorrected inv. mass for pions", 64, 0.3, 3.5);
    
    outfile->cd();
}

bool ConnectInput(int argc, char** argv) 
{
    int fileId = -1;
    upcTree = 0x0;

    if( argc == 2 || argc == 3)
    {
        const string& input = argv[1];
        if(input.find(".root") != string::npos && argc == 2)
        {
            cout << "Input from root file: "<< input << endl;
            infile = TFile::Open(input.c_str(), "read");
            if(!infile)
            {
                cout<< "Couldn't open input root file..."<<endl;
                return false;
            } 
            upcTree = dynamic_cast<TTree*>( infile->Get("mUPCTree") );
        }else if(input.find(".list") != string::npos )
        {
            cout << "Input from chain" << endl;
            upcChain = new TChain("mUPCTree");
            ifstream instr(input.c_str());
            if (!instr.is_open())
            {
                cout<< "Couldn't open: "<<input.c_str()<<endl;
                return false;
            }

            string line;
            int lineId=0;
            if(argc == 3)
                fileId = atoi(argv[2]);

            while(getline(instr, line)) 
            {
                if(fileId==lineId || fileId== -1)
                {
                    upcChain->AddFile(line.c_str());
                    infile = TFile::Open(line.c_str(), "read");
                    if(!infile)
                    {
                        cout<< "Couldn't open: "<<line.c_str()<<endl;
                        return false;
                    } 
                }
                lineId++;
            }
            instr.close();
            upcTree = dynamic_cast<TTree*>( upcChain );
        }
    }
    
    if(!upcTree) 
        return false;

    upcTree->SetBranchAddress("mUPCEvent", &upcEvt);
    upcTree->SetBranchAddress("mRPEvent", &rpEvt); 

    return true;
}//ConnectInput


//_____________________________________________________________________________
TFile *CreateOutputTree(const string& out) {

    TFile *outputFile = TFile::Open(out.c_str(), "recreate");
    if(!outputFile) 
        return 0x0;

    //standard reconstructed tree
    recTree = new TTree("recTree", "recTree");

    // Invariant mass of state
    recTree->Branch("invMass", &invMass, "invMass/D");

    // Vertex info
    recTree->Branch("vertexZ", &vertexeZ, "vertexZ/D");

    // event info
    recTree->Branch("runNumber", &runNumber, "runNumber/i");

    // Setting background Tree
    bcgTree = recTree->CloneTree(0);
    bcgTree->SetName("Background");

    return outputFile;

}//CreateOutputTree

bool CheckTriggers()
{

    bool CPTtrigger = false;
    for(int var = 0; var < nTriggers; ++var)
    {
        if(upcEvt->isTrigger(triggerID[var]))
        {
            hTriggerBits->Fill(var);
            //Checked if it is CPT trigger
            for (int i = 0; i < *(&CEPtriggers + 1) - CEPtriggers; ++i)
                if(triggerID[var] == CEPtriggers[i])
                    CPTtrigger=true;
        }
    }

    return CPTtrigger;
}

