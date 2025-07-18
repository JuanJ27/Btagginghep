#ifndef JETANALYZER_H
#define JETANALYZER_H

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <iostream>
#include "MyClass.C"

class JetAnalyzer {
public:
    // Constructor y Destructor
    JetAnalyzer(const std::vector<std::string>& inputFiles);
    ~JetAnalyzer();

    // Métodos principales
    void InitializeHistograms();
    void LoopEvents();
    void DrawHistograms();
    void SaveHistograms(const std::string& outputDir);

private:
    // Métodos auxiliares
    void ProcessEvent(Long64_t entry);

    // Miembros de datos
    TChain* fChain;
    MyClass* t;
    Long64_t nentries;

    // Histogramas
    // Histograma de jets por evento
    TH1F* hJetsPerEvent;

    // Histograma de pT de los jets
    TH1F* hJetPT[4];

    // Histograma de Eta de los jets
    TH1F* hJetEta[4];

    // Histograma de Phi de los jets
    TH1F* hJetPhi[4];

    // Delta R entre los primeros 4 jets, por parejas
    TH1F* hDeltaRPar[6];

    // Número de partículas cargadas por jet
    TH1F* hChargedParticles[4];

    // Número de partículas neutras por jet
    TH1F* hNeutralsParticles[4];

    // Fracción de pT cargado del jet
    TH1F* hChargedPTFraction[4];

    // Fracción de pT neutro del jet
    TH1F* hNeutralPTFraction[4];

    // pT promedio de las partículas de cada jet
    TH1F* hAveragePT[4];

    // Número de partículas con pT por debajo del pT promedio por cada jet
    TH1F* hParticlesBelowAvgPT[4];

    // Número de partículas con pT por encima del pT promedio por cada jet
    TH1F* hParticlesAboveAvgPT[4];

    // pT(par_max_pT) / pT(j_r)
    TH1F* hMaxPTRatio[4];

    // pT(par_min_pT) / pT(j_r)
    TH1F* hMinPTRatio[4];

    // pT(par_max_DR) / pT(j_r)
    TH1F* hMaxDRRatio[4];

    // pT(par_min_DR) / pT(j_r)
    TH1F* hMinDRRatio[4];

    // DeltaR(par_max_pT, j_r)
    TH1F* hDeltaRMaxPT[4];

    // DeltaR(par_min_pT, j_r)
    TH1F* hDeltaRMinPT[4];

    // DeltaR(par_max_DR, j_r)
    TH1F* hDeltaRMaxDR[4];

    // DeltaR(par_min_DR, j_r)
    TH1F* hDeltaRMinDR[4];

    // pT(par_max_pT) - pT(par_min_pT)
    TH1F* hPTDifference[4];

    // R al 50% del pT total del jet
    TH1F* hR50PercentPT[4];

    // R al 95% del pT total del jet
    TH1F* hR95PercentPT[4];

    // Histograma 2D de DeltaR vs Porcentaje acumulado de pT
    TH2F* hCumulativePT_vs_DeltaR[4];

    // Nuevos Histogramas 2D
    TH2F* hPT_vs_Eta[4]; // pT vs Eta para los primeros 4 jets
    TH2F* hCharged_vs_NeutralParticles[4]; // Número de Partículas Cargadas vs Neutras
    TH2F* hChargedPTFraction_vs_NeutralPTFraction[4]; // Fracción de pT Cargado vs Neutro
    TH2F* hAveragePT_vs_TotalParticles[4]; // pT Promedio vs Número Total de Partículas
    //TH2F* hDeltaR_vs_PTDifference[6]; // Delta R vs Diferencia en pT entre pares de jets
    TH2F* hMaxPTRatio_vs_DeltaRMaxPT[4]; // pT(par_max_pT) vs Delta R(par_max_pT, j_r)
    TH2F* hR50_vs_R95[4]; // R50% vs R95%
    
    // Nuevos histogramas relacionados con D0 y DZ
    TH2F* hCumulativePT_vs_DZTrack[4];   // Porcentaje acumulado de pT vs DZ del Track
    TH2F* hCumulativePT_vs_D0Track[4];    // Porcentaje acumulado de pT vs D0 del Track
    TH2F* hDeltaR_vs_DZTrack[4];          // DeltaR entre Track y Jet vs DZ del Track
    TH2F* hDeltaR_vs_D0Track[4];          // DeltaR entre Track y Jet vs D0 del Track

    // Vectores para almacenar jets y partículas
    std::vector<TLorentzVector> jets;
    std::vector<TLorentzVector> particles;
};

#endif // JETANALYZER_H
