/*
Helper class to produce data cards and make plots of the data cards.
Originally Written by Merijn van de klundert, mvandekl@cern.ch

Editted for tt-Channel by Vinay. Changes mentioned by VK:
*/


#ifndef Configurable_tt_H
#define Configurable_tt_H


#include "TH2F.h"
#include "TSystem.h"
#include <iostream>
#include "TDirectory.h"
#include "TMath.h"
#include "TROOT.h"


using namespace std;
class Configurable_tt : public TObject{
  
 public:
  Configurable_tt(){};
  Configurable_tt(int category_,int sample_){
    cat=category_;
    sample=sample_;
  }
  ~Configurable_tt(){};

  //Static Members
  inline static const vector<TString> SampleNames{"data_obs","EmbedZTT", "ZL", "TTT","jetFakes", "VVT","EWKZ", 
      "ggH_sm_htt125", "ggH_ps_htt125", "ggH_mm_htt125","qqH_sm_htt125", "qqH_ps_htt125", "qqH_mm_htt125","QCD"};
  
  inline static const vector<TString> BaseCatNames{"ggH","qqH","ZTT","QCD","misc","NoCat"};
  inline static vector<TString> CatNames;
  inline static TString InputFileDir; 
  inline static TString OutputFileDir;
  inline static TString OutputDir;
  inline static TString OutputFileDataCard;
  inline static TString OutputDirPreFit;
  inline static int era;
  inline static int decaymode_1;
  inline static int decaymode_2;
  inline static int CPMethod;
  inline static int Observable;
  inline static bool do_Embedding;
  inline static TString ObservableDraw; 
  inline static TString ObservablePath;
  inline static TString ObservableXaxis;
  inline static TString ObservableYaxis; 
  inline static  double TTnorm = 1.0;
  inline static  double Wnorm  = 1.0;  
  inline static  TString topweight="1*";   
  inline static  TString qcdweight="1*";
  inline static  TString zptmassweight="1*";

  //for the histogram
  TH1F* OneDimHisto=NULL;
  TH1F* OneDimHistoSS=NULL;
  TH2F* TwoDimHisto=NULL;
  TH2F* TwoDimHistoSS=NULL;

  uint cat;
  int sample;
  TString CatName, SampleName;
  TString InputFileName;
  int TauChConstCut_1;
  int TauChConstCut_2;
  TString selcut;
  TString Isocut;
  TString AntiIsocut_1,AntiIsocut_2;
  TString FinalCut;

  TString EventsWeight;

  TString HistoName, HistoTitle;
  TString HistoName2D, HistoTitle2D;
  TString HistoNameSS, HistoTitleSS;
  TString HistoNameSS2D, HistoTitleSS2D;
  //fixed in InitialiseXandYAxis
  float xmin, xmax;
  int nbinsx;
  double * xbins=NULL; //initialise to NULL, avoids crash when delete and not assigned. 
  float ymin, ymax;
  int nbinsy;
  double * ybins=NULL;


  inline static void InitOutputBaseDir(){
    if(era>2018||era<2016){throw std::invalid_argument("Era given not valid, please choose 2016, 2017 or 2018");}
    if(decaymode_1>11||era<-1){throw std::invalid_argument("decaymode given not valid, please enter -1, 0,1, or 11");}
    if(CPMethod>1||era<0){throw std::invalid_argument("CPMethod given not valid, please choose 0 or 1");}

    OutputDir="";
    OutputDir+=era;
    OutputDir+="/";
    OutputDir+="DM_";
    OutputDir+=decaymode_1;
    OutputDir+=decaymode_2;
    OutputDir+="/";

  }
  inline static void InitOutputFileDataCard(){
    OutputFileDataCard=OutputDir;
    OutputFileDataCard+=ObservablePath;
    OutputFileDataCard+="/";
    OutputFileDataCard+="Datacard.root";
  }
  inline static void InitOutputPrefitPlotDir(){
    OutputDirPreFit="Plots_Prefit/";
    OutputDirPreFit+=OutputDir;
  }
  void Initialise(){
    SampleName=SampleNames[sample];
    InputFileName+=InputFileDir;
    InputFileName+=SampleNametoInputFile(SampleName);
    CatName=CatNames[cat];
    EventsWeight = "xsec_lumi_weight*idisoweight_1 * idisoweight_2 * trigweight*mcweight*zptweight*topptweight*";
    if(!do_Embedding)EventsWeight += "puweight*";
    if(CPMethod==0){
      TauChConstCut_1=40;
      TauChConstCut_2=40;
    }
    selcut="";
    if(decaymode_1!=-1 && decaymode_2!=-1){
      selcut+="&&tau_MVADM_1==";
      selcut+=decaymode_1;
      selcut+="&&tau_MVADM_2==";
      selcut+=decaymode_2;}
    selcut+="&&pt_1>";
    selcut+=TauChConstCut_1;
    selcut+="&&pt_2>";
    selcut+=TauChConstCut_2;
    Isocut="&&mva17_1 >0.5 &&mva17_2 >0.5 &&Prompt_pT > 50";
    AntiIsocut_1 = "&& mva17_1 < 0.5 && mva17_2 > 0.5 &&Prompt_pT > 50";
    if(decaymode_1!=-1 && decaymode_2!=-1){
      AntiIsocut_1+="&&tau_MVADM_1==";
      AntiIsocut_1+=decaymode_1;
      AntiIsocut_1+="&&tau_MVADM_2==";
      AntiIsocut_1+=decaymode_2;}
    AntiIsocut_1+="&&pt_1>";
    AntiIsocut_1+=TauChConstCut_1;
    AntiIsocut_1+="&&pt_2>";
    AntiIsocut_1+=TauChConstCut_2;
    AntiIsocut_2 = "&& mva17_2 < 0.5 && mva17_1 > 0.5 &&Prompt_pT > 50";
    if(sample==0){
      FinalCut="(os>0.5"+selcut+")";}
    else if(sample==1){//Ztt
      TString isZTT="((gen_match_1==5&&gen_match_2==5)||(gen_match_2==5&&gen_match_1==5))*";
      FinalCut=isZTT+EventsWeight+zptmassweight+"(os>0.5"+selcut+")";}
    else if(sample==2){//Zll
      TString isZLL="((gen_match_1==1&&gen_match_2==1)||((gen_match_2==2&&gen_match_1==2))*puweight*";
      FinalCut=isZLL+EventsWeight+zptmassweight+"(os>0.5"+selcut+")";}
    else if(sample==3){
      TString isZLL="((gen_match_1==1&&gen_match_2==1)||((gen_match_2==2&&gen_match_1==2))*";
      FinalCut=isZLL+EventsWeight+zptmassweight+"(os>0.5"+selcut+")";}
    else if(sample==4){
      FinalCut=EventsWeight+"(os>0.5"+AntiIsocut_1+")*ff_nom";}
    else{
      FinalCut=EventsWeight+"(os>0.5"+selcut+")";}
  

    InitialiseXandYAxis(cat);
    InitHistogram(); 
  }
  TString SampleNametoInputFile(TString sample){
    TString out="tt-NOMINAL_ntuple_";
    TString SampleAppendix="";
    if(sample.Contains("data_obs")){
      SampleAppendix="Data";}
    else if(sample.Contains("ggH")){
      SampleAppendix="ggH125";}
    else if(sample.Contains("qqH")){
      SampleAppendix="qqH125";}
    else if(sample.Contains("ZTT")){
      SampleAppendix="Embedded";}
    else if(sample.Contains("ZL")){
      SampleAppendix="DY";}
   
    else{
      SampleAppendix+=sample;}
    out+=SampleAppendix;
    out+=".root";
    return out;
  }
  void InitHistogram(){ 
    HistoName=CatName+SampleName;
    HistoName2D=CatName+SampleName+"_2D";
    OneDimHisto=new TH1F(HistoName,HistoName,nbinsx,xbins);
    OneDimHisto->Sumw2();
    OneDimHisto->SetXTitle(ObservableXaxis);
    OneDimHisto->SetYTitle(ObservableYaxis);
    TwoDimHisto=new TH2F(HistoName2D,HistoName2D,nbinsx,xbins,nbinsy,ybins);
    TwoDimHisto->Sumw2();
    TwoDimHisto->SetXTitle(ObservableXaxis);
    TwoDimHisto->SetYTitle(ObservableYaxis);
  }
  void Write1DHistoToDir(TDirectory *CatDir){
    CatDir->cd();    
    TH1F* tmphisto=new TH1F();
    OneDimHisto->Copy(*tmphisto);
    tmphisto->SetName(SampleName);
    tmphisto->Write();
    delete tmphisto;}
  void Write2DHistoToDir(TDirectory *CatDir){
    CatDir->cd();    
    TH2F* tmphisto=new TH2F();
    TwoDimHisto->Copy(*tmphisto);
    tmphisto->SetName(SampleName+"_2D");
    tmphisto->Write();
    delete tmphisto;}
  inline static TH1F * Unfold(TH2F * histInput) {
    int nBinsX = histInput->GetNbinsX();
    int nBinsY = histInput->GetNbinsY();
    int nBins = nBinsX * nBinsY;

    TString NameNew = TString(histInput->GetName())+TString("_unfolded");
    TH1F * histOutput = new TH1F(NameNew,"",nBins,0,float(nBins));
    int iBin = 1;
    for (int j=1; j<=nBinsY; ++j) {
      for (int i=1; i<=nBinsX; ++i) {
	histOutput->SetBinContent(iBin,histInput->GetBinContent(i,j));
	histOutput->SetBinError(iBin,histInput->GetBinError(i,j));
	iBin++;}}
    return histOutput;}
  inline static TH1D * ProjectonXaxis(TH2F * histInput){
    TH1D * histOutput=histInput->ProjectionX();
    return histOutput;}

  inline static void InitCatNames(){
    CatNames.clear();
    TString BaseCat="";  
    for(uint i=0; i<BaseCatNames.size();i++){
      BaseCat="tt_";
      BaseCat+=era;
      BaseCat+="_";
      BaseCat+=BaseCatNames[i];
      if(decaymode_1==-1 && decaymode_2==-1) BaseCat+="_tautau";
      if(decaymode_1==0 && decaymode_2 ==0) BaseCat+="_pipi";
      if((decaymode_1==1&&decaymode_2==0)||(decaymode_2==1&&decaymode_1==0)) BaseCat+="_pirho";
      if(decaymode_1==11 || decaymode_2 ==11) BaseCat+="_other";
      if(CPMethod==0) BaseCat+="_IP";
      if(CPMethod==1) BaseCat+="_DP";
      CatNames.push_back(BaseCat);
    }}
  inline static void InitObservableNames(){
    if(Observable>21||Observable<0){throw std::invalid_argument("For index given, no observable is in the lookup table, exiting");}
    ObservableYaxis="Occ.";//correct for nearly all
    ObservableDraw="predicted_prob:"; //we always draw two dimensionally w.r.t. the nn score, below we add the specific observable for x-axis. Adjust ObservableDraw below if you want to unroll in a different observable   
    if(Observable==0){ ObservablePath="pt_1";   ObservableXaxis="pt #mu [GeV]";}
    if(Observable==1){ ObservablePath="pt_2";   ObservableXaxis="pt #tau [GeV]";}
    if(Observable==2){ ObservablePath="jpt_1";   ObservableXaxis="Leading jet pT [GeV]"; }
    if(Observable==3){ ObservablePath="jpt_2";   ObservableXaxis="Trailing jet pT [GeV]"; }
    if(Observable==4){ ObservablePath="bpt_1";   ObservableXaxis="Leading b-jet pT [GeV]"; }
    if(Observable==5){ ObservablePath="bpt_2";   ObservableXaxis="Trailing b-jet pT [GeV]"; }   
    if(Observable==6){ ObservablePath="njets";   ObservableXaxis="Jet multiplicity"; }
    if(Observable==7){ ObservablePath="nbtag";   ObservableXaxis="B-jet multiplicity"; }
    if(Observable==8){ ObservablePath="mt_1";   ObservableXaxis="#mu transverse mass [GeV]"; }
    if(Observable==9){ ObservablePath="pt_tt";   ObservableXaxis="Di-#tau pT [GeV]"; }
    if(Observable==10){ ObservablePath="mjj";   ObservableXaxis="Jet pair invar. mass [GeV]"; }
    if(Observable==11){ ObservablePath="jdeta";   ObservableXaxis="Jet pair #Delta#eta"; }
    if(Observable==12){ ObservablePath="dijetpt";   ObservableXaxis="Di-jet pT [GeV]"; }
    if(Observable==13){ ObservablePath="m_vis";   ObservableXaxis="Visible mass [GeV]"; }
    if(Observable==14){ ObservablePath="eta_1";   ObservableXaxis="#mu #eta"; }
    if(Observable==15){ ObservablePath="eta_2";   ObservableXaxis="#tau #eta"; }
    if(Observable==16){ ObservablePath="met_rcmr";   ObservableXaxis="met [GeV]"; }
    if(Observable==17){ ObservablePath="pt_sv";   ObservableXaxis="pT SVfit [GeV]"; }
    if(Observable==18){ ObservablePath="";   ObservableXaxis=""; }//empty slots!
    if(Observable==19){ ObservablePath="";   ObservableXaxis=""; }
    //observables > 20 reserved for NN score and CP angles
    if(Observable==20){ ObservablePath="predicted_prob";   ObservableXaxis="NN output"; ObservableYaxis="dN/d(NN output)";}
    if(Observable==21){ ObservablePath="acotautau_0"; ObservablePath+=CPMethod;  ObservableXaxis="#phi CP"; }
   
    ObservableDraw+=ObservablePath;
  }
  void InitialiseXandYAxis(int cat){

    //we only declare here the number of bins and ranges for obs. 0-19. Below, axes are initialised quidistantly in case we don't do something special for obs. 20 and 21.
    if(Observable==0){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==1){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==2){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==3){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==4){ nbinsx  = 16; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==5){ nbinsx  = 16; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==6){ nbinsx  = 5; xmin=0;xmax=5; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==7){ nbinsx  = 5; xmin=0;xmax=5; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==8){ nbinsx  = 22; xmin=0;xmax=60; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==9){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==10){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==11){ nbinsx  = 8; xmin=0;xmax=5.6; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==12){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==13){ nbinsx  = 30; xmin=0;xmax=160;ymin=0; ymax=1; nbinsy=1;}
    if(Observable==14){ nbinsx  = 50; xmin=-2.5;xmax=2.5; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==15){ nbinsx  = 50; xmin=-2.5;xmax=2.5; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==16){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==17){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
    if(Observable==18){ ymin=0; ymax=1; nbinsy=1;}
    if(Observable==19){ ymin=0; ymax=1; nbinsy=1;}
    
    //below declare binnings for the NN axes for different categories.       
    if(Observable==20){
      ymin=0; ymax=1; nbinsy=1; //for plotting nn score, we just need one bin on y-axis
      AxesLUT(cat, nbinsx, xbins, decaymode_1,decaymode_2, CPMethod);
      
    }
    if(Observable==21){
      nbinsx  = 5; xmin=0;xmax=2*TMath::Pi();
      AxesLUT(cat,nbinsy, ybins, decaymode_1,decaymode_2, CPMethod);}
    
    if(xbins==NULL){
      xbins=new double[nbinsx+1];
      double step=(xmax-xmin)/nbinsx;
      for(int bc=0;bc<nbinsx+1;bc++){
	xbins[bc]=xmin+bc*step;}}
   
    if(ybins==NULL){
      ybins=new double[nbinsy+1];
      double step=(ymax-ymin)/nbinsy;     
      for(int bc=0;bc<nbinsy+1;bc++){
	ybins[bc]=ymin+bc*step;}}
  }
  void AxesLUT(int CAT, int & NBINS, double*& axis, int DM1, int DM2, int CM){
    //Here declare everything as function of CAT and whatever dependence you like to add; pt Higgs, ...!!
    vector<double> tmp;

    //pipi
    if(DM1==0 && DM2==0){
      if(CAT==0) tmp={0.2,0.3,0.4,0.5,0.6,0.7,1.0};//ggH
      if(CAT==1) tmp={0.2,0.4,0.6,0.7,0.8,0.9,1};//qqH
      if(CAT==2) tmp={0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0};//ZTT
      if(CAT==3) tmp={0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0};//ZL
      if(CAT==4) tmp={0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0};//QCD
      if(CAT==5) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1};
    }

    //pi+rho in dp method
    if(DM1==1&&CM==1){
      if(CAT==0) tmp={0, 0.25, 0.3,0.35, 1}; 
      if(CAT==1) tmp={0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1}; 
      if(CAT==2) tmp={0, 0.35, 0.45, 0.55, 0.65, 1}; 
      if(CAT==3) tmp={0, 0.25, 0.3, 0.4, 0.45, 0.55, 1}; 
      if(CAT==4) tmp={0, 0.25, 0.35,0.4, 0.45, 1}; 
      if(CAT==5) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1}; 
      if(CAT==6) tmp={0, 0.3, 0.4, 0.5, 0.6, 1}; 
      if(CAT==7) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.8, 1}; 
      if(CAT==8) tmp={0, 0.3, 0.6, 1}; //DNN score range in case we pick no category, i.e. it is void observable in for the DNN score or CP angles, but need to provide something.
    }

    if(DM1==11){//other category. We only need signal cats, rest is to avoid crash..
      if(CAT==0) tmp={0,0.3,1}; 
      if(CAT==1) tmp={0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1}; 
      if(CAT==2) tmp={0, 0.35, 0.45, 0.55, 0.65, 1}; 
      if(CAT==3) tmp={0, 0.25, 0.3, 0.4, 0.45, 0.55, 1}; 
      if(CAT==4) tmp={0, 0.25, 0.35,0.4, 0.45,0.5,0.6,1}; 
      if(CAT==5) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1}; 
      if(CAT==6) tmp={0, 0.3, 0.4, 0.5, 0.6, 1}; 
      if(CAT==7) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.8, 1}; 
      if(CAT==8) tmp={0, 0.3, 0.6, 1}; //DNN score range in case we pick no category, i.e. it is void observable in for the DNN score or CP angles, but need to provide something.
    }
   
    NBINS=tmp.size()-1;
    axis=new double[tmp.size()];
    for(uint index=0; index<tmp.size();index++){
      axis[index]=tmp[index];}
  }
  ClassDef (Configurable_tt,1)
    };
#endif
