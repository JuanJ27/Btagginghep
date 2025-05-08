#include "JetAnalyzer.h"
#include <algorithm>
#include <cmath>
#include <fstream>

JetAnalyzer::JetAnalyzer(const std::vector<std::string>& inputFiles) {
    // Inicializar TChain y agregar archivos
    fChain = new TChain("Delphes", "");
    for (const auto& file : inputFiles) {
        fChain->Add(file.c_str());
    }

    // Crear instancia de MyClass para acceder a las ramas del arbol
    t = new MyClass(fChain);
    nentries = t->fChain->GetEntries();

    // Inicializar histogramas
    InitializeHistograms();
}

JetAnalyzer::~JetAnalyzer() {
    // Liberar memoria de histogramas
    delete hJetsPerEvent;

    for (int i = 0; i < 4; i++) {
        delete hJetPT[i];
        delete hJetEta[i];
        delete hJetPhi[i];
        delete hChargedParticles[i];
        delete hNeutralsParticles[i];
        delete hChargedPTFraction[i];
        delete hNeutralPTFraction[i];
        delete hAveragePT[i];
        delete hParticlesBelowAvgPT[i];
        delete hParticlesAboveAvgPT[i];
        delete hMaxPTRatio[i];
        delete hMinPTRatio[i];
        delete hMaxDRRatio[i];
        delete hMinDRRatio[i];
        delete hDeltaRMaxPT[i];
        delete hDeltaRMinPT[i];
        delete hDeltaRMaxDR[i];
        delete hDeltaRMinDR[i];
        delete hPTDifference[i];
        delete hR50PercentPT[i];
        delete hR95PercentPT[i];
        delete hCumulativePT_vs_DeltaR[i];
    }

    for (int i = 0; i < 6; i++) {
        delete hDeltaRPar[i];
    }

    // Liberar memoria de TChain y MyClass
    delete t;
    delete fChain;
}

void JetAnalyzer::InitializeHistograms() {
    // Histograma de jets por evento
    hJetsPerEvent = new TH1F("hJetsPerEvent", "Numero de Jets por Evento", 10, 0, 10);

    // Inicializar histogramas para los primeros 4 jets
    for (int i = 0; i < 4; i++) {
        // pT de los jets
        hJetPT[i] = new TH1F(Form("hJetPT%d", i), "pT del Jet", 100, 5, 210);
        hJetPT[i]->SetLineColor(i+1);
        hJetPT[i]->SetLineWidth(2);

        // Eta de los jets
        hJetEta[i] = new TH1F(Form("hJetEta%d", i), "Eta del Jet", 100, -5.5, 5.5);
        hJetEta[i]->SetLineColor(i+1);
        hJetEta[i]->SetLineWidth(2);

        // Phi de los jets
        hJetPhi[i] = new TH1F(Form("hJetPhi%d", i), "Phi del Jet", 100, -3.5, 3.5);
        hJetPhi[i]->SetLineColor(i+1);
        hJetPhi[i]->SetLineWidth(2);

        // Numero de particulas cargadas por jet
        hChargedParticles[i] = new TH1F(Form("hChargedParticles%d", i), "Numero de particulas cargadas del Jet", 27, 0, 27);
        hChargedParticles[i]->SetLineColor(i+1);
        hChargedParticles[i]->SetLineWidth(2);

        // Numero de particulas neutras por jet
        hNeutralsParticles[i] = new TH1F(Form("hNeutralsParticles%d", i), "Numero de particulas neutras del Jet", 27, 0, 27);
        hNeutralsParticles[i]->SetLineColor(i+1);
        hNeutralsParticles[i]->SetLineWidth(2);

        // Fraccion de pT cargado del jet
        hChargedPTFraction[i] = new TH1F(Form("hChargedPTFraction%d", i), "Fraccion de pT cargado del Jet", 100, 0, 1);
        hChargedPTFraction[i]->SetLineColor(i+1);
        hChargedPTFraction[i]->SetLineWidth(2);

        // Fraccion de pT neutro del jet
        hNeutralPTFraction[i] = new TH1F(Form("hNeutralPTFraction%d", i), "Fraccion de pT neutro del Jet", 100, 0, 1);
        hNeutralPTFraction[i]->SetLineColor(i+1);
        hNeutralPTFraction[i]->SetLineWidth(2);

        // pT promedio de las particulas de cada jet
        hAveragePT[i] = new TH1F(Form("hAveragePT%d", i), "pT promedio de las particulas del Jet", 100, 0, 60);
        hAveragePT[i]->SetLineColor(i+1);
        hAveragePT[i]->SetLineWidth(2);

        // Numero de particulas con pT por debajo del pT promedio por cada jet
        hParticlesBelowAvgPT[i] = new TH1F(Form("hParticlesBelowAvgPT%d", i), "Numero de particulas con pT < pT promedio del Jet", 120, 0, 120);
        hParticlesBelowAvgPT[i]->SetLineColor(i+1);
        hParticlesBelowAvgPT[i]->SetLineWidth(2);

        // Numero de particulas con pT por encima del pT promedio por cada jet
        hParticlesAboveAvgPT[i] = new TH1F(Form("hParticlesAboveAvgPT%d", i), "Numero de particulas con pT > pT promedio del Jet", 60, 0, 60);
        hParticlesAboveAvgPT[i]->SetLineColor(i+1);
        hParticlesAboveAvgPT[i]->SetLineWidth(2);

        // pT(par_max_pT) / pT(j_r)
        hMaxPTRatio[i] = new TH1F(Form("hMaxPTRatio%d", i), "pT(par_max_pT) / pT(j_r) del Jet", 100, 0, 3.5);
        hMaxPTRatio[i]->SetLineColor(i+1);
        hMaxPTRatio[i]->SetLineWidth(2);

        // pT(par_min_pT) / pT(j_r)
        hMinPTRatio[i] = new TH1F(Form("hMinPTRatio%d", i), "pT(par_min_pT) / pT(j_r) del Jet", 100, 0, 0.15);
        hMinPTRatio[i]->SetLineColor(i+1);
        hMinPTRatio[i]->SetLineWidth(2);

        // pT(par_max_DR) / pT(j_r)
        hMaxDRRatio[i] = new TH1F(Form("hMaxDRRatio%d", i), "pT(par_max_DR) / pT(j_r) del Jet", 100, 0, 2.5);
        hMaxDRRatio[i]->SetLineColor(i+1);
        hMaxDRRatio[i]->SetLineWidth(2);

        // pT(par_min_DR) / pT(j_r)
        hMinDRRatio[i] = new TH1F(Form("hMinDRRatio%d", i), "pT(par_min_DR) / pT(j_r) del Jet", 100, 0, 2.5);
        hMinDRRatio[i]->SetLineColor(i+1);
        hMinDRRatio[i]->SetLineWidth(2);

        // DeltaR(par_max_pT, j_r)
        hDeltaRMaxPT[i] = new TH1F(Form("hDeltaRMaxPT%d", i), "DeltaR(par_max_pT, j_r) del Jet", 100, 0, 0.4);
        hDeltaRMaxPT[i]->SetLineColor(i+1);
        hDeltaRMaxPT[i]->SetLineWidth(2);

        // DeltaR(par_min_pT, j_r)
        hDeltaRMinPT[i] = new TH1F(Form("hDeltaRMinPT%d", i), "DeltaR(par_min_pT, j_r) del Jet", 100, 0, 0.4);
        hDeltaRMinPT[i]->SetLineColor(i+1);
        hDeltaRMinPT[i]->SetLineWidth(2);

        // DeltaR(par_max_DR, j_r)
        hDeltaRMaxDR[i] = new TH1F(Form("hDeltaRMaxDR%d", i), "DeltaR(par_max_DR, j_r) del Jet", 100, 0, 0.4);
        hDeltaRMaxDR[i]->SetLineColor(i+1);
        hDeltaRMaxDR[i]->SetLineWidth(2);

        // DeltaR(par_min_DR, j_r)
        hDeltaRMinDR[i] = new TH1F(Form("hDeltaRMinDR%d", i), "DeltaR(par_min_DR, j_r) del Jet", 100, 0, 0.4);
        hDeltaRMinDR[i]->SetLineColor(i+1);
        hDeltaRMinDR[i]->SetLineWidth(2);

        // pT(par_max_pT) - pT(par_min_pT)
        hPTDifference[i] = new TH1F(Form("hPTDifference%d", i), "pT(par_max_pT) - pT(par_min_pT) del Jet", 100, 0, 200);
        hPTDifference[i]->SetLineColor(i+1);
        hPTDifference[i]->SetLineWidth(2);

        // R al 50% del pT total del jet
        hR50PercentPT[i] = new TH1F(Form("hR50PercentPT%d", i), "DeltaR para 50% pT total del Jet", 100, 0, 0.4);
        hR50PercentPT[i]->SetLineColor(i+1);
        hR50PercentPT[i]->SetLineWidth(2);

        // R al 95% del pT total del jet
        hR95PercentPT[i] = new TH1F(Form("hR95PercentPT%d", i), "DeltaR para 95% pT total del Jet", 100, 0, 0.4);
        hR95PercentPT[i]->SetLineColor(i+1);
        hR95PercentPT[i]->SetLineWidth(2);

        // Histograma 2D de DeltaR vs Porcentaje acumulado de pT
        hCumulativePT_vs_DeltaR[i] = new TH2F(Form("hCumulativePT_vs_DeltaR%d", i+1),
                                              Form("DeltaR vs Porcentaje acumulado de pT del Jet %d", i+1),
                                              25, 0, 0.4,    // Eje X: DeltaR (0 - 0.4)
                                              25, 0, 1);   // Eje Y: Porcentaje acumulado de pT (0% - 100%)
        hCumulativePT_vs_DeltaR[i]->SetLineColor(i+1);
        hCumulativePT_vs_DeltaR[i]->SetLineWidth(2);

            // DeltaR vs D0
        hDeltaR_vs_D0[i] = new TH2F(Form("hDeltaR_vs_D0_%d", i+1),
                                    Form("DeltaR vs D0 del Track para Jet %d", i+1),
                                    50, -5, 5,     // Eje X: D0
                                    50, 0, 0.4);   // Eje Y: DeltaR
        hDeltaR_vs_D0[i]->GetXaxis()->SetTitle("D_{0} (cm)");
        hDeltaR_vs_D0[i]->GetYaxis()->SetTitle("#DeltaR(Track, Jet)");
        hDeltaR_vs_D0[i]->SetLineColor(i+1);
        hDeltaR_vs_D0[i]->SetLineWidth(2);

        // Cumulative PT vs D0
        hCumulativePT_vs_D0[i] = new TH2F(Form("hCumulativePT_vs_D0_%d", i+1),
                                        Form("Cumulative pT vs D0 del Track para Jet %d", i+1),
                                        50, -5, 5,     // Eje X: D0
                                        50, 0, 1);     // Eje Y: % Acumulado pT
        hCumulativePT_vs_D0[i]->GetXaxis()->SetTitle("D_{0} (cm)");
        hCumulativePT_vs_D0[i]->GetYaxis()->SetTitle("Cumulative pT fraction");
        hCumulativePT_vs_D0[i]->SetLineColor(i+1);
        hCumulativePT_vs_D0[i]->SetLineWidth(2);

    // pT vs Eta de los jets
        hPT_vs_Eta[i] = new TH2F(Form("hPT_vs_Eta%d", i), Form("pT vs Eta del Jet %d", i+1),
                                50, -5.5, 5.5, // Eje X: Eta
                                50, 0, 100);   // Eje Y: pT
        hPT_vs_Eta[i]->GetXaxis()->SetTitle("Eta");
        hPT_vs_Eta[i]->GetYaxis()->SetTitle("pT (GeV)");
        hPT_vs_Eta[i]->SetLineColor(i+1);
        hPT_vs_Eta[i]->SetLineWidth(2);

    // Numero de particulas cargadas vs neutras por jet
        hCharged_vs_NeutralParticles[i] = new TH2F(Form("hCharged_vs_NeutralParticles%d", i),
                                                Form("Cargadas vs Neutras del Jet %d", i+1),
                                                15, 0, 15, // Eje X: Cargadas
                                                15, 0, 15); // Eje Y: Neutras
        hCharged_vs_NeutralParticles[i]->GetXaxis()->SetTitle("Numero de Particulas Cargadas");
        hCharged_vs_NeutralParticles[i]->GetYaxis()->SetTitle("Numero de Particulas Neutras");
        hCharged_vs_NeutralParticles[i]->SetLineColor(i+1);
        hCharged_vs_NeutralParticles[i]->SetLineWidth(2);


    // Fraccion de pT Cargado vs Neutro
        hChargedPTFraction_vs_NeutralPTFraction[i] = new TH2F(Form("hChargedPTFraction_vs_NeutralPTFraction%d", i),
                                                                Form("Fraccion de pT Cargado vs Neutro del Jet %d", i+1),
                                                                10, 0, 1, // Eje X: Fraccion Cargado
                                                                10, 0, 1); // Eje Y: Fraccion Neutro
        hChargedPTFraction_vs_NeutralPTFraction[i]->GetXaxis()->SetTitle("Fraccion de pT Cargado");
        hChargedPTFraction_vs_NeutralPTFraction[i]->GetYaxis()->SetTitle("Fraccion de pT Neutro");
        hChargedPTFraction_vs_NeutralPTFraction[i]->SetLineColor(i+1);
        hChargedPTFraction_vs_NeutralPTFraction[i]->SetLineWidth(2);

    // pT Promedio vs Numero Total de Particulas
        hAveragePT_vs_TotalParticles[i] = new TH2F(Form("hAveragePT_vs_TotalParticles%d", i),
                                                Form("pT Promedio vs Total Particulas del Jet %d", i+1),
                                                30, 0, 30, // Eje X: Total Particulas
                                                40, 0, 20);  // Eje Y: pT Promedio
        hAveragePT_vs_TotalParticles[i]->GetXaxis()->SetTitle("Numero Total de Particulas");
        hAveragePT_vs_TotalParticles[i]->GetYaxis()->SetTitle("pT Promedio (GeV)");
        hAveragePT_vs_TotalParticles[i]->SetLineColor(i+1);
        hAveragePT_vs_TotalParticles[i]->SetLineWidth(2);

    /* Delta R vs Diferencia en pT entre pares de jets
        hDeltaR_vs_PTDifference[i] = new TH2F(Form("hDeltaR_vs_PTDifference%d", i),
                                            Form("DeltaR vs Diferencia en pT para Par %d", i+1),
                                            50, 0, 10, // Eje X: DeltaR
                                            100, 0, 200); // Eje Y: Diferencia en pT
        hDeltaR_vs_PTDifference[i]->GetXaxis()->SetTitle("DeltaR");
        hDeltaR_vs_PTDifference[i]->GetYaxis()->SetTitle("|pT1 - pT2| (GeV)");
        hDeltaR_vs_PTDifference[i]->SetLineColor(i+1);
        hDeltaR_vs_PTDifference[i]->SetLineWidth(2);*/

    // pT(par_max_pT) vs Delta R(par_max_pT, j_r)
        hMaxPTRatio_vs_DeltaRMaxPT[i] = new TH2F(Form("hMaxPTRatio_vs_DeltaRMaxPT%d", i),
                                                Form("pT(par_max_pT) / pT(j_r) vs DeltaR(par_max_pT, j_r) del Jet %d", i+1),
                                                10, 0, 2, // Eje X: pT Ratio
                                                10, 0, 0.4); // Eje Y: DeltaR
        hMaxPTRatio_vs_DeltaRMaxPT[i]->GetXaxis()->SetTitle("pT(par_{max}^{PT}) / pT(j_{r})");
        hMaxPTRatio_vs_DeltaRMaxPT[i]->GetYaxis()->SetTitle("DeltaR(par_{max}^{PT}, j_{r})");
        hMaxPTRatio_vs_DeltaRMaxPT[i]->SetLineColor(i+1);
        hMaxPTRatio_vs_DeltaRMaxPT[i]->SetLineWidth(2);

    // R50% vs R95%
        hR50_vs_R95[i] = new TH2F(Form("hR50_vs_R95%d", i),
                                Form("R50%% vs R95%% del Jet %d", i+1),
                                10, 0, 0.4, // Eje X: R50%
                                10, 0, 0.4); // Eje Y: R95%
        hR50_vs_R95[i]->GetXaxis()->SetTitle("R50% del pT total");
        hR50_vs_R95[i]->GetYaxis()->SetTitle("R95% del pT total");
        hR50_vs_R95[i]->SetLineColor(i+1);
        hR50_vs_R95[i]->SetLineWidth(2);

    }

    // Delta R entre los primeros 4 jets, por parejas
    for (int i = 0; i < 6; i++) {
        hDeltaRPar[i] = new TH1F(Form("hDeltaRPar%d", i), "Delta R entre los primeros 4 jets, por parejas", 50, 0, 10);
        hDeltaRPar[i]->SetLineColor(i+1);
        hDeltaRPar[i]->SetLineWidth(2);
    }
}

void JetAnalyzer::LoopEvents() {
    std::cout << "Total Entries: " << nentries << std::endl;
    Long64_t nTen = nentries / 10; // Para imprimir el porcentaje de avance

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        ProcessEvent(jentry);

        // Mostrar progreso
        if (jentry % nTen == 0)
            std::cout << 10 * (jentry / nTen) << "%-" << std::flush;
        if (jentry == nentries - 1)
            std::cout << "100%" << std::endl;
    }
}

void JetAnalyzer::ProcessEvent(Long64_t entry) {
    // Cargar el evento
    Long64_t ientry = t->LoadTree(entry);
    if (ientry < 0) return;
    t->fChain->GetEntry(entry);

    if (t->Jet_size == 0) return;

    // Limpiar vectores
    jets.clear();
    particles.clear();

    // Llenar vectores de jets
    for (Int_t i = 0; i < t->Jet_size; i++) {
        TLorentzVector jetVector;
        jetVector.SetPtEtaPhiM(t->Jet_PT[i], t->Jet_Eta[i], t->Jet_Phi[i], t->Jet_Mass[i]);
        jets.push_back(jetVector);
    }

    // Llenar vectores de particulas
    for (Int_t i = 0; i < t->Particle_size; i++) {
        TLorentzVector particleVector;
        particleVector.SetPtEtaPhiM(t->Particle_PT[i], t->Particle_Eta[i], t->Particle_Phi[i], t->Particle_Mass[i]);
        particles.push_back(particleVector);
    }

    // Histograma de jets por evento
    hJetsPerEvent->Fill(t->Jet_size);



    // Procesamiento detallado para cada jet
    for (Int_t i = 0; i < std::min(4, t->Jet_size); i++) {

        // Histograma de pT, Eta y Phi de los jets
        hJetPT[i]->Fill(t->Jet_PT[i]);
        hJetEta[i]->Fill(t->Jet_Eta[i]);
        hJetPhi[i]->Fill(t->Jet_Phi[i]);

        // Delta R entre los primeros 4 jets, por parejas
        for (int j = i + 1; j < std::min(4, t->Jet_size); j++) {
            double deltaR_par = jets[i].DeltaR(jets[j]);
            int index = i * 3 - (i * (i + 1)) / 2 + j - 1;
            if (index >= 0 && index < 6) {
                hDeltaRPar[index]->Fill(deltaR_par);
            }
        }
        
        // Histograma de Numero de particulas cargadas y neutras por jet
        hChargedParticles[i]->Fill(t->Jet_NCharged[i]);
        hNeutralsParticles[i]->Fill(t->Jet_NNeutrals[i]);

        // Variables del jet actual
        TLorentzVector jetVector = jets[i];
        Double_t jetPT = jetVector.Pt();

        // Inicializacion de variables
        TLorentzVector chargedPTSum(0, 0, 0, 0);
        TLorentzVector neutralPTSum(0, 0, 0, 0);
        Double_t sumPT = 0.0;
        Int_t countParticles = 0;

        TLorentzVector maxPTParticle(0, 0, 0, 0);
        TLorentzVector minPTParticle(1e9, 0, 0, 0);
        TLorentzVector maxDRParticle(0, 0, 0, 0);
        TLorentzVector minDRParticle(0, 0, 0, 0);
        Double_t maxDR = -1;
        Double_t minDR = 1e9;

        // Estructura para almacenar informacion de particulas
        struct ParticleInfo {
            TLorentzVector vector;
            Double_t deltaR;
            Int_t charge;
            Double_t pt;
        };

        std::vector<ParticleInfo> particlesInJet;

        // Bucle sobre particulas para el jet actual
        for (size_t j = 0; j < particles.size(); j++) {
            const TLorentzVector& particleVector = particles[j];
            Double_t deltaR = jetVector.DeltaR(particleVector);

            if (deltaR < 0.4) {
                // Almacenar informacion de la particula
                particlesInJet.push_back({particleVector, deltaR, t->Particle_Charge[j], particleVector.Pt()});

                // Sumar pT y contar particulas
                sumPT += particleVector.Pt();
                countParticles++;

                // Actualizar particulas con maximo y minimo pT
                if (particleVector.Pt() > maxPTParticle.Pt()) {
                    maxPTParticle = particleVector;
                }
                if (particleVector.Pt() < minPTParticle.Pt()) {
                    minPTParticle = particleVector;
                }

                // Actualizar particulas con maximo y minimo DeltaR
                if (deltaR > maxDR) {
                    maxDR = deltaR;
                    maxDRParticle = particleVector;
                }
                if (deltaR < minDR) {
                    minDR = deltaR;
                    minDRParticle = particleVector;
                }

                // Sumar pT cargado y neutro
                if (t->Particle_Charge[j] != 0) {
                    chargedPTSum += particleVector;
                } else {
                    neutralPTSum += particleVector;
                }
            }
        }

        struct TrackInfo {
        Double_t pt;
        Double_t d0;
        Double_t deltaR;
        };

        std::vector<TrackInfo> tracksInJet;
        Double_t totalTrackPT = 0.0;

        for (Int_t j = 0; j < t->Track_size; j++) {
            // Crear vector del Track
            TLorentzVector trackVector;
            trackVector.SetPtEtaPhiM(t->Track_PT[j], t->Track_Eta[j], t->Track_Phi[j], 0);  // Track mass ≈ 0

            Double_t deltaR = jetVector.DeltaR(trackVector);

            if (deltaR < 0.4) {  // Asociamos el Track al Jet si está dentro del cono
                tracksInJet.push_back({t->Track_PT[j], t->Track_D0[j], deltaR});
                totalTrackPT += t->Track_PT[j];

                // Llenar directamente DeltaR vs D0
                hDeltaR_vs_D0[i]->Fill(t->Track_D0[j], deltaR);
            }
        }

        // Calcular fracciones de pT
        Double_t totalPT = chargedPTSum.Pt() + neutralPTSum.Pt();
        if (totalPT == 0) totalPT = 1e-9; // Evitar division por cero

        Double_t chargedPTFraction = chargedPTSum.Pt() / totalPT;
        Double_t neutralPTFraction = neutralPTSum.Pt() / totalPT;
        Double_t averagePT = (countParticles > 0) ? (sumPT / countParticles) : 0;

        Double_t maxPTRatio = (jetPT > 0) ? (maxPTParticle.Pt() / jetPT) : 0;
        Double_t minPTRatio = (jetPT > 0) ? (minPTParticle.Pt() / jetPT) : 0;
        Double_t maxDRRatio = (jetPT > 0) ? (maxDRParticle.Pt() / jetPT) : 0;
        Double_t minDRRatio = (jetPT > 0) ? (minDRParticle.Pt() / jetPT) : 0;

        Double_t deltaRMaxPT = jetVector.DeltaR(maxPTParticle);
        Double_t deltaRMinPT = jetVector.DeltaR(minPTParticle);
        Double_t deltaRMaxDR = jetVector.DeltaR(maxDRParticle);
        Double_t deltaRMinDR = jetVector.DeltaR(minDRParticle);

        Double_t ptDifference = maxPTParticle.Pt() - minPTParticle.Pt();

        // Contar particulas por encima y por debajo del pT promedio
        Int_t particlesBelowAvgPT = 0;
        Int_t particlesAboveAvgPT = 0;

        for (const auto& pInfo : particlesInJet) {
            if (pInfo.pt < averagePT) {
                particlesBelowAvgPT++;
            } else {
                particlesAboveAvgPT++;
            }
        }

        // Llenar histogramas
        hChargedPTFraction[i]->Fill(chargedPTFraction);
        hNeutralPTFraction[i]->Fill(neutralPTFraction);
        hAveragePT[i]->Fill(averagePT);
        hParticlesBelowAvgPT[i]->Fill(particlesBelowAvgPT);
        hParticlesAboveAvgPT[i]->Fill(particlesAboveAvgPT);
        hMaxPTRatio[i]->Fill(maxPTRatio);
        hMinPTRatio[i]->Fill(minPTRatio);
        hMaxDRRatio[i]->Fill(maxDRRatio);
        hMinDRRatio[i]->Fill(minDRRatio);
        hDeltaRMaxPT[i]->Fill(deltaRMaxPT);
        hDeltaRMinPT[i]->Fill(deltaRMinPT);
        hDeltaRMaxDR[i]->Fill(deltaRMaxDR);
        hDeltaRMinDR[i]->Fill(deltaRMinDR);
        hPTDifference[i]->Fill(ptDifference);

        // Verificar que sumPT no sea cero para evitar divisiones por cero
        if (sumPT == 0.0) continue;

        // Ordenar particulas por DeltaR
        std::sort(particlesInJet.begin(), particlesInJet.end(), [](const ParticleInfo& a, const ParticleInfo& b) {
            return a.deltaR < b.deltaR;
        });

        // Calcular el porcentaje acumulado de pT y llenar el histograma 2D
        Double_t cumulativePT = 0.0;
        for (const auto& pInfo : particlesInJet) {
            cumulativePT += pInfo.pt;
            Double_t cumulativePTFraction = cumulativePT / sumPT;
            hCumulativePT_vs_DeltaR[i]->Fill(pInfo.deltaR, cumulativePTFraction);
        }

        // Ordenar por |D0| para el acumulativo
        std::sort(tracksInJet.begin(), tracksInJet.end(), [](const TrackInfo& a, const TrackInfo& b) {
            return fabs(a.d0) < fabs(b.d0);
        });

        // Calcular el porcentaje acumulado de pT vs D0
        cumulativePT = 0.0;
        for (const auto& track : tracksInJet) {
            cumulativePT += track.pt;
            Double_t cumulativeFraction = cumulativePT / totalTrackPT;
            hCumulativePT_vs_D0[i]->Fill(track.d0, cumulativeFraction);
        }

        // Calcular R para el 50% y 95% del pT total del jet
        Double_t aimPT50 = sumPT * 0.5;
        Double_t aimPT95 = sumPT * 0.95;
        cumulativePT = 0.0;
        Double_t r50PercentPT = 0.0;
        Double_t r95PercentPT = 0.0;

        for (const auto& pInfo : particlesInJet) {
            cumulativePT += pInfo.pt;
            if (cumulativePT >= aimPT50 && r50PercentPT == 0.0) {
                r50PercentPT = pInfo.deltaR;
            }
            if (cumulativePT >= aimPT95 && r95PercentPT == 0.0) {
                r95PercentPT = pInfo.deltaR;
                break;
            }
        }

        hR50PercentPT[i]->Fill(r50PercentPT);
        hR95PercentPT[i]->Fill(r95PercentPT);

        // 1. pT vs Eta
        hPT_vs_Eta[i]->Fill(t->Jet_Eta[i], t->Jet_PT[i]);

        // 2. Numero de particulas cargadas vs neutras
        hCharged_vs_NeutralParticles[i]->Fill(t->Jet_NCharged[i], t->Jet_NNeutrals[i]);

        // 3. Fraccion de pT cargado vs neutro
        hChargedPTFraction_vs_NeutralPTFraction[i]->Fill(chargedPTFraction, neutralPTFraction);

        // 4. pT Promedio vs Numero Total de Particulas
        Int_t totalParticles = t->Jet_NCharged[i] + t->Jet_NNeutrals[i];
        hAveragePT_vs_TotalParticles[i]->Fill(totalParticles, averagePT);

        /* 5. Delta R vs Diferencia en pT entre pares de jets
        // Este histograma requiere informacion sobre los pares de jets. Puedes llenarlo dentro del bucle que calcula hDeltaRPar
        // Por ejemplo:
        for (int j = i + 1; j < std::min(4, t->Jet_size); j++) {
            double deltaR_par = jets[i].DeltaR(jets[j]);
            double ptDifference = std::abs(t->Jet_PT[i] - t->Jet_PT[j]);
            int index = i * 3 - (i * (i + 1)) / 2 + j - 1;
            if (index >= 0 && index < 6) {
                hDeltaR_vs_PTDifference[index]->Fill(deltaR_par, ptDifference);
            }
        }*/

        // 6. pT(par_max_pT) vs Delta R(par_max_pT, j_r)
        hMaxPTRatio_vs_DeltaRMaxPT[i]->Fill(maxPTRatio, deltaRMaxPT);

        // 7. R50% vs R95%
        hR50_vs_R95[i]->Fill(r50PercentPT, r95PercentPT);

    }
}

void JetAnalyzer::DrawHistograms() {
    // Configuracion general
    gStyle->SetOptStat(0); // Desactivar la caja de estadisticas

    // Histograma de jets por evento
    TCanvas* cJets = new TCanvas("cJets", "Jets por evento", 600, 400);
    gPad->SetLogy(); // Escala logaritmica en el eje y
    hJetsPerEvent->Draw();
    cJets->SaveAs("plots/JetsPerEvent.png");
    delete cJets;

    // Histograma de pT de los jets
    TCanvas* c1 = new TCanvas("c1", "pT de los Jets", 600, 400);
    gPad->SetLogy();
    hJetPT[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hJetPT[i]->Draw("SAME");
    }
    TLegend* legend = new TLegend(0.8, 0.7, 0.9, 0.9);
    legend->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legend->AddEntry(hJetPT[i], Form("Jet %d", i+1), "l");
    }
    legend->Draw();
    c1->SaveAs("plots/PTjets.png");
    delete legend;
    delete c1;

    // Histograma de Eta de los jets
    TCanvas* cEta = new TCanvas("cEta", "Eta de los Jets", 600, 400);
    gPad->SetLogy();
    hJetEta[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hJetEta[i]->Draw("SAME");
    }
    TLegend* legendEta = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendEta->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendEta->AddEntry(hJetEta[i], Form("Jet %d", i+1), "l");
    }
    legendEta->Draw();
    cEta->SaveAs("plots/Etajets.png");
    delete legendEta;
    delete cEta;

    // Histograma de Phi de los jets
    TCanvas* cPhi = new TCanvas("cPhi", "Phi de los Jets", 600, 400);
    hJetPhi[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hJetPhi[i]->Draw("SAME");
    }
    TLegend* legendPhi = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendPhi->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendPhi->AddEntry(hJetPhi[i], Form("Jet %d", i+1), "l");
    }
    legendPhi->Draw();
    cPhi->SaveAs("plots/Phijets.png");
    delete legendPhi;
    delete cPhi;

    // Delta R entre los primeros 4 jets, por parejas
    TCanvas* cRPar = new TCanvas("cRPar", "R entre las primeras 4 parejas de jets", 600, 400);
    gPad->SetLogy();
    hDeltaRPar[0]->Draw();
    for (int i = 1; i < 6; i++) {
        hDeltaRPar[i]->Draw("SAME");
    }
    TLegend* legendRPar = new TLegend(0.7, 0.55, 0.9, 0.9);
    legendRPar->SetHeader("DeltaR pares de Jets", "C");
    legendRPar->AddEntry(hDeltaRPar[0], "Entre 1 y 2", "l");
    legendRPar->AddEntry(hDeltaRPar[1], "Entre 1 y 3", "l");
    legendRPar->AddEntry(hDeltaRPar[2], "Entre 1 y 4", "l");
    legendRPar->AddEntry(hDeltaRPar[3], "Entre 2 y 3", "l");
    legendRPar->AddEntry(hDeltaRPar[4], "Entre 2 y 4", "l");
    legendRPar->AddEntry(hDeltaRPar[5], "Entre 3 y 4", "l");
    legendRPar->Draw();
    cRPar->SaveAs("plots/RPar.png");
    delete legendRPar;
    delete cRPar;

    // Histograma de numero de particulas cargadas por jet
    TCanvas* cNCharged = new TCanvas("cNCharged", "Numero de particulas cargadas por jet", 600, 400);
    gPad->SetLogy();
    hChargedParticles[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hChargedParticles[i]->Draw("SAME");
    }
    TLegend* legendNCharged = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendNCharged->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendNCharged->AddEntry(hChargedParticles[i], Form("Jet %d", i+1), "l");
    }
    legendNCharged->Draw();
    cNCharged->SaveAs("plots/NCharged.png");
    delete legendNCharged;
    delete cNCharged;

    // Histograma de numero de particulas neutras por jet
    TCanvas* cNNeutrals = new TCanvas("cNNeutrals", "Numero de particulas neutras por jet", 600, 400);
    gPad->SetLogy();
    hNeutralsParticles[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hNeutralsParticles[i]->Draw("SAME");
    }
    TLegend* legendNNeutrals = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendNNeutrals->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendNNeutrals->AddEntry(hNeutralsParticles[i], Form("Jet %d", i+1), "l");
    }
    legendNNeutrals->Draw();
    cNNeutrals->SaveAs("plots/NNeutrals.png");
    delete legendNNeutrals;
    delete cNNeutrals;

    // Histograma de fraccion de pT cargado del jet
    TCanvas* cChargedPTFraction = new TCanvas("cChargedPTFraction", "Fraccion de pT cargado del Jet", 600, 400);
    gPad->SetLogy();
    hChargedPTFraction[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hChargedPTFraction[i]->Draw("SAME");
    }
    TLegend* legendChargedPTFraction = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendChargedPTFraction->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendChargedPTFraction->AddEntry(hChargedPTFraction[i], Form("Jet %d", i+1), "l");
    }
    legendChargedPTFraction->Draw();
    cChargedPTFraction->SaveAs("plots/ChargedPTFraction.png");
    delete legendChargedPTFraction;
    delete cChargedPTFraction;

    // Histograma de fraccion de pT neutro del jet
    TCanvas* cNeutralPTFraction = new TCanvas("cNeutralPTFraction", "Fraccion de pT neutro del Jet", 600, 400);
    gPad->SetLogy();
    hNeutralPTFraction[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hNeutralPTFraction[i]->Draw("SAME");
    }
    TLegend* legendNeutralPTFraction = new TLegend(0.1, 0.7, 0.2, 0.9);
    legendNeutralPTFraction->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendNeutralPTFraction->AddEntry(hNeutralPTFraction[i], Form("Jet %d", i+1), "l");
    }
    legendNeutralPTFraction->Draw();
    cNeutralPTFraction->SaveAs("plots/NeutralPTFraction.png");
    delete legendNeutralPTFraction;
    delete cNeutralPTFraction;

    // Histograma de pT promedio de las particulas del jet
    TCanvas* cAveragePT = new TCanvas("cAveragePT", "pT promedio de las particulas del Jet", 600, 400);
    gPad->SetLogy();
    hAveragePT[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hAveragePT[i]->Draw("SAME");
    }
    TLegend* legendAveragePT = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendAveragePT->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendAveragePT->AddEntry(hAveragePT[i], Form("Jet %d", i+1), "l");
    }
    legendAveragePT->Draw();
    cAveragePT->SaveAs("plots/AveragePT.png");
    delete legendAveragePT;
    delete cAveragePT;

    // Histograma de numero de particulas con pT por debajo del promedio
    TCanvas* cParticlesBelowAvgPT = new TCanvas("cParticlesBelowAvgPT", "Numero de particulas con pT < pT promedio del Jet", 600, 400);
    gPad->SetLogy();
    hParticlesBelowAvgPT[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hParticlesBelowAvgPT[i]->Draw("SAME");
    }
    TLegend* legendParticlesBelowAvgPT = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendParticlesBelowAvgPT->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendParticlesBelowAvgPT->AddEntry(hParticlesBelowAvgPT[i], Form("Jet %d", i+1), "l");
    }
    legendParticlesBelowAvgPT->Draw();
    cParticlesBelowAvgPT->SaveAs("plots/ParticlesBelowAvgPT.png");
    delete legendParticlesBelowAvgPT;
    delete cParticlesBelowAvgPT;

    // Histograma de numero de particulas con pT por encima del promedio
    TCanvas* cParticlesAboveAvgPT = new TCanvas("cParticlesAboveAvgPT", "Numero de particulas con pT > pT promedio del Jet", 600, 400);
    gPad->SetLogy();
    hParticlesAboveAvgPT[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hParticlesAboveAvgPT[i]->Draw("SAME");
    }
    TLegend* legendParticlesAboveAvgPT = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendParticlesAboveAvgPT->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendParticlesAboveAvgPT->AddEntry(hParticlesAboveAvgPT[i], Form("Jet %d", i+1), "l");
    }
    legendParticlesAboveAvgPT->Draw();
    cParticlesAboveAvgPT->SaveAs("plots/ParticlesAboveAvgPT.png");
    delete legendParticlesAboveAvgPT;
    delete cParticlesAboveAvgPT;

    // Histograma de pT(par_max_pT) / pT(j_r)
    TCanvas* cMaxPTRatio = new TCanvas("cMaxPTRatio", "pT(par_max_pT) / pT(j_r) del Jet", 600, 400);
    gPad->SetLogy();
    hMaxPTRatio[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hMaxPTRatio[i]->Draw("SAME");
    }
    TLegend* legendMaxPTRatio = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendMaxPTRatio->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendMaxPTRatio->AddEntry(hMaxPTRatio[i], Form("Jet %d", i+1), "l");
    }
    legendMaxPTRatio->Draw();
    cMaxPTRatio->SaveAs("plots/MaxPTRatio.png");
    delete legendMaxPTRatio;
    delete cMaxPTRatio;

    // Histograma de pT(par_min_pT) / pT(j_r)
    TCanvas* cMinPTRatio = new TCanvas("cMinPTRatio", "pT(par_min_pT) / pT(j_r) del Jet", 600, 400);
    gPad->SetLogy();
    hMinPTRatio[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hMinPTRatio[i]->Draw("SAME");
    }
    TLegend* legendMinPTRatio = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendMinPTRatio->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendMinPTRatio->AddEntry(hMinPTRatio[i], Form("Jet %d", i+1), "l");
    }
    legendMinPTRatio->Draw();
    cMinPTRatio->SaveAs("plots/MinPTRatio.png");
    delete legendMinPTRatio;
    delete cMinPTRatio;

    // Histograma de pT(par_max_DR) / pT(j_r)
    TCanvas* cMaxDRRatio = new TCanvas("cMaxDRRatio", "pT(par_max_DR) / pT(j_r) del Jet", 600, 400);
    gPad->SetLogy();
    hMaxDRRatio[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hMaxDRRatio[i]->Draw("SAME");
    }
    TLegend* legendMaxDRRatio = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendMaxDRRatio->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendMaxDRRatio->AddEntry(hMaxDRRatio[i], Form("Jet %d", i+1), "l");
    }
    legendMaxDRRatio->Draw();
    cMaxDRRatio->SaveAs("plots/MaxDRRatio.png");
    delete legendMaxDRRatio;
    delete cMaxDRRatio;

    // Histograma de pT(par_min_DR) / pT(j_r)
    TCanvas* cMinDRRatio = new TCanvas("cMinDRRatio", "pT(par_min_DR) / pT(j_r) del Jet", 600, 400);
    gPad->SetLogy();
    hMinDRRatio[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hMinDRRatio[i]->Draw("SAME");
    }
    TLegend* legendMinDRRatio = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendMinDRRatio->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendMinDRRatio->AddEntry(hMinDRRatio[i], Form("Jet %d", i+1), "l");
    }
    legendMinDRRatio->Draw();
    cMinDRRatio->SaveAs("plots/MinDRRatio.png");
    delete legendMinDRRatio;
    delete cMinDRRatio;

    // Histograma de DeltaR(par_max_pT, j_r)
    TCanvas* cDeltaRMaxPT = new TCanvas("cDeltaRMaxPT", "DeltaR(par_max_pT, j_r) del Jet", 600, 400);
    gPad->SetLogy();
    hDeltaRMaxPT[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hDeltaRMaxPT[i]->Draw("SAME");
    }
    TLegend* legendDeltaRMaxPT = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendDeltaRMaxPT->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendDeltaRMaxPT->AddEntry(hDeltaRMaxPT[i], Form("Jet %d", i+1), "l");
    }
    legendDeltaRMaxPT->Draw();
    cDeltaRMaxPT->SaveAs("plots/DeltaRMaxPT.png");
    delete legendDeltaRMaxPT;
    delete cDeltaRMaxPT;

    // Histograma de DeltaR(par_min_pT, j_r)
    TCanvas* cDeltaRMinPT = new TCanvas("cDeltaRMinPT", "DeltaR(par_min_pT, j_r) del Jet", 600, 400);
    gPad->SetLogy();
    hDeltaRMinPT[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hDeltaRMinPT[i]->Draw("SAME");
    }
    TLegend* legendDeltaRMinPT = new TLegend(0.1, 0.7, 0.2, 0.9);
    legendDeltaRMinPT->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendDeltaRMinPT->AddEntry(hDeltaRMinPT[i], Form("Jet %d", i+1), "l");
    }
    legendDeltaRMinPT->Draw();
    cDeltaRMinPT->SaveAs("plots/DeltaRMinPT.png");
    delete legendDeltaRMinPT;
    delete cDeltaRMinPT;

    // Histograma de DeltaR(par_max_DR, j_r)
    TCanvas* cDeltaRMaxDR = new TCanvas("cDeltaRMaxDR", "DeltaR(par_max_DR, j_r) del Jet", 600, 400);
    gPad->SetLogy();
    hDeltaRMaxDR[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hDeltaRMaxDR[i]->Draw("SAME");
    }
    TLegend* legendDeltaRMaxDR = new TLegend(0.1, 0.7, 0.2, 0.9);
    legendDeltaRMaxDR->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendDeltaRMaxDR->AddEntry(hDeltaRMaxDR[i], Form("Jet %d", i+1), "l");
    }
    legendDeltaRMaxDR->Draw();
    cDeltaRMaxDR->SaveAs("plots/DeltaRMaxDR.png");
    delete legendDeltaRMaxDR;
    delete cDeltaRMaxDR;

    // Histograma de DeltaR(par_min_DR, j_r)
    TCanvas* cDeltaRMinDR = new TCanvas("cDeltaRMinDR", "DeltaR(par_min_DR, j_r) del Jet", 600, 400);
    gPad->SetLogy();
    hDeltaRMinDR[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hDeltaRMinDR[i]->Draw("SAME");
    }
    TLegend* legendDeltaRMinDR = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendDeltaRMinDR->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendDeltaRMinDR->AddEntry(hDeltaRMinDR[i], Form("Jet %d", i+1), "l");
    }
    legendDeltaRMinDR->Draw();
    cDeltaRMinDR->SaveAs("plots/DeltaRMinDR.png");
    delete legendDeltaRMinDR;
    delete cDeltaRMinDR;

    // Histograma de pT(par_max_pT) - pT(par_min_pT)
    TCanvas* cPTDifference = new TCanvas("cPTDifference", "pT(par_max_pT) - pT(par_min_pT) del Jet", 600, 400);
    gPad->SetLogy();
    hPTDifference[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hPTDifference[i]->Draw("SAME");
    }
    TLegend* legendPTDifference = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendPTDifference->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendPTDifference->AddEntry(hPTDifference[i], Form("Jet %d", i+1), "l");
    }
    legendPTDifference->Draw();
    cPTDifference->SaveAs("plots/PTDifference.png");
    delete legendPTDifference;
    delete cPTDifference;

    // Histograma de R al 50% del pT total del jet
    TCanvas* cR50 = new TCanvas("cR50", "R al 50% del pT total del jet", 600, 400);
    gPad->SetLogy();
    hR50PercentPT[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hR50PercentPT[i]->Draw("SAME");
    }
    TLegend* legendR50 = new TLegend(0.8, 0.7, 0.9, 0.9);
    legendR50->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendR50->AddEntry(hR50PercentPT[i], Form("Jet %d", i+1), "l");
    }
    legendR50->Draw();
    cR50->SaveAs("plots/R50.png");
    delete legendR50;
    delete cR50;

    // Histograma de R al 95% del pT total del jet
    TCanvas* cR95 = new TCanvas("cR95", "R al 95% del pT total del jet", 600, 400);
    gPad->SetLogy();
    hR95PercentPT[0]->Draw();
    for (int i = 1; i < 4; i++) {
        hR95PercentPT[i]->Draw("SAME");
    }
    TLegend* legendR95 = new TLegend(0.1, 0.7, 0.2, 0.9);
    legendR95->SetHeader("Jets", "C");
    for (int i = 0; i < 4; i++) {
        legendR95->AddEntry(hR95PercentPT[i], Form("Jet %d", i+1), "l");
    }
    legendR95->Draw();
    cR95->SaveAs("plots/R95.png");
    delete legendR95;
    delete cR95;

    // Dibujar y guardar histogramas 2D de DeltaR vs Porcentaje acumulado de pT
    for (int i = 0; i < 4; i++) {
        TCanvas* cCumulativePT_vs_DeltaR = new TCanvas(Form("cCumulativePT_vs_DeltaR%d", i+1),
                                                       Form("Porcentaje acumulado de pT vs DeltaR del Jet %d", i+1),
                                                       600, 400);
        gStyle->SetNumberContours(10); // Numero de colores a utilizar
        gPad->SetLogz(); // Escala logaritmica en el eje z
        hCumulativePT_vs_DeltaR[i]->GetXaxis()->SetTitle("#DeltaR");
        hCumulativePT_vs_DeltaR[i]->GetYaxis()->SetTitle("Porcentaje acumulado de pT");
        hCumulativePT_vs_DeltaR[i]->Draw("SURF2");
        cCumulativePT_vs_DeltaR->SaveAs(Form("plots/CumulativePT_vs_DeltaR_Jet%d.png", i+1));
        hCumulativePT_vs_DeltaR[i]->Draw("CONT1");
        cCumulativePT_vs_DeltaR->SaveAs(Form("plots/CONT_CumulativePT_vs_DeltaR_Jet%d.png", i+1));
        delete cCumulativePT_vs_DeltaR;
    }

    for (int i = 0; i < 4; i++) {
        // DeltaR vs D0
        TCanvas* cDeltaR_vs_D0 = new TCanvas(Form("cDeltaR_vs_D0_%d", i+1), "DeltaR vs D0", 800, 600);
        hDeltaR_vs_D0[i]->Draw("COLZ");
        cDeltaR_vs_D0->SaveAs(Form("plots/DeltaR_vs_D0_Jet%d.png", i+1));
        delete cDeltaR_vs_D0;

        // Cumulative PT vs D0
        TCanvas* cCumulativePT_vs_D0 = new TCanvas(Form("cCumulativePT_vs_D0_%d", i+1), "Cumulative PT vs D0", 800, 600);
        hCumulativePT_vs_D0[i]->Draw("COLZ");
        cCumulativePT_vs_D0->SaveAs(Form("plots/CumulativePT_vs_D0_Jet%d.png", i+1));
        delete cCumulativePT_vs_D0;
    }


    // pT vs Eta de los jets
    for (int i = 0; i < 4; i++) {
        TCanvas* cPT_vs_Eta = new TCanvas(Form("cPT_vs_Eta%d", i), Form("pT vs Eta del Jet %d", i+1), 800, 600);
        gStyle->SetNumberContours(10); // Numero de colores a utilizar
        gPad->SetLogz(); // Escala logaritmica en el eje z
        hPT_vs_Eta[i]->GetXaxis()->SetTitle("Eta");
        hPT_vs_Eta[i]->GetYaxis()->SetTitle("pT (GeV)");
        hPT_vs_Eta[i]->Draw("COLZ");
        cPT_vs_Eta->SaveAs(Form("plots/PT_vs_Eta_Jet%d.png", i+1));
        delete cPT_vs_Eta;
    }

    // Numero de particulas cargadas vs neutras por jet
    for (int i = 0; i < 4; i++) {
        TCanvas* cCharged_vs_Neutral = new TCanvas(Form("cCharged_vs_Neutral%d", i), Form("Cargadas vs Neutras del Jet %d", i+1), 800, 600);
        gStyle->SetNumberContours(10); // Numero de colores a utilizar
        gPad->SetLogz(); // Escala logaritmica en el eje z
        hCharged_vs_NeutralParticles[i]->GetXaxis()->SetTitle("Numero de Particulas Cargadas");
        hCharged_vs_NeutralParticles[i]->GetYaxis()->SetTitle("Numero de Particulas Neutras");
        hCharged_vs_NeutralParticles[i]->Draw("COLZ");
        cCharged_vs_Neutral->SaveAs(Form("plots/Charged_vs_Neutral_Jet%d.png", i+1));
        hCharged_vs_NeutralParticles[i]->Draw("CONT1");
        cCharged_vs_Neutral->SaveAs(Form("plots/CONT_Charged_vs_Neutral_Jet%d.png", i+1));
        delete cCharged_vs_Neutral;
    }

    // Fraccion de pT cargado vs neutro
    for (int i = 0; i < 4; i++) {
        TCanvas* cChargedPTFraction_vs_NeutralPTFraction = new TCanvas(Form("cChargedPTFraction_vs_NeutralPTFraction%d", i),
                                                                        Form("Fraccion de pT Cargado vs Neutro del Jet %d", i+1),
                                                                        800, 600);
        gStyle->SetNumberContours(10); // Numero de colores a utilizar
        gPad->SetLogz(); // Escala logaritmica en el eje z
        hChargedPTFraction_vs_NeutralPTFraction[i]->GetXaxis()->SetTitle("Fraccion de pT Cargado");
        hChargedPTFraction_vs_NeutralPTFraction[i]->GetYaxis()->SetTitle("Fraccion de pT Neutro");
        hChargedPTFraction_vs_NeutralPTFraction[i]->Draw("COLZ");
        cChargedPTFraction_vs_NeutralPTFraction->SaveAs(Form("plots/ChargedPTFraction_vs_NeutralPTFraction_Jet%d.png", i+1));
        delete cChargedPTFraction_vs_NeutralPTFraction;
    }

    // pT Promedio vs Numero Total de Particulas
    for (int i = 0; i < 4; i++) {
        TCanvas* cAveragePT_vs_TotalParticles = new TCanvas(Form("cAveragePT_vs_TotalParticles%d", i),
                                                            Form("pT Promedio vs Total Particulas del Jet %d", i+1),
                                                            800, 600);
        gStyle->SetNumberContours(10); // Numero de colores a utilizar
        gPad->SetLogz(); // Escala logaritmica en el eje z
        hAveragePT_vs_TotalParticles[i]->GetXaxis()->SetTitle("Numero Total de Particulas");
        hAveragePT_vs_TotalParticles[i]->GetYaxis()->SetTitle("pT Promedio (GeV)");
        hAveragePT_vs_TotalParticles[i]->Draw("COLZ");
        cAveragePT_vs_TotalParticles->SaveAs(Form("plots/AveragePT_vs_TotalParticles_Jet%d.png", i+1));
        hAveragePT_vs_TotalParticles[i]->Draw("CONT1");
        cAveragePT_vs_TotalParticles->SaveAs(Form("plots/CONT_AveragePT_vs_TotalParticles_Jet%d.png", i+1));
        delete cAveragePT_vs_TotalParticles;
    }

    /* Delta R vs Diferencia en pT entre pares de jets
    for (int i = 0; i < 6; i++) {
        TCanvas* cDeltaR_vs_PTDifference = new TCanvas(Form("cDeltaR_vs_PTDifference%d", i), Form("DeltaR vs Diferencia en pT para Par %d", i+1), 800, 600);
        hDeltaR_vs_PTDifference[i]->GetXaxis()->SetTitle("DeltaR");
        hDeltaR_vs_PTDifference[i]->GetYaxis()->SetTitle("|pT1 - pT2| (GeV)");
        hDeltaR_vs_PTDifference[i]->Draw("COLZ");
        cDeltaR_vs_PTDifference->SaveAs(Form("plots/DeltaR_vs_PTDifference_Par%d.png", i+1));
        delete cDeltaR_vs_PTDifference;
    }*/

    // pT(par_max_pT) vs Delta R(par_max_pT, j_r)
    for (int i = 0; i < 4; i++) {
        TCanvas* cMaxPTRatio_vs_DeltaRMaxPT = new TCanvas(Form("cMaxPTRatio_vs_DeltaRMaxPT%d", i),
                                                        Form("pT(par_max_pT) / pT(j_r) vs DeltaR(par_max_pT, j_r) del Jet %d", i+1),
                                                        800, 600);
        gStyle->SetNumberContours(10); // Numero de colores a utilizar
        gPad->SetLogz(); // Escala logaritmica en el eje z
        hMaxPTRatio_vs_DeltaRMaxPT[i]->GetXaxis()->SetTitle("pT(par_{max}^{PT}) / pT(j_{r})");
        hMaxPTRatio_vs_DeltaRMaxPT[i]->GetYaxis()->SetTitle("DeltaR(par_{max}^{PT}, j_{r})");
        hMaxPTRatio_vs_DeltaRMaxPT[i]->Draw("COLZ");
        cMaxPTRatio_vs_DeltaRMaxPT->SaveAs(Form("plots/MaxPTRatio_vs_DeltaRMaxPT_Jet%d.png", i+1));
        hMaxPTRatio_vs_DeltaRMaxPT[i]->Draw("CONT1");
        cMaxPTRatio_vs_DeltaRMaxPT->SaveAs(Form("plots/CONT_MaxPTRatio_vs_DeltaRMaxPT_Jet%d.png", i+1));
        delete cMaxPTRatio_vs_DeltaRMaxPT;
    }

    // R50% vs R95%
    for (int i = 0; i < 4; i++) {
        TCanvas* cR50_vs_R95 = new TCanvas(Form("cR50_vs_R95%d", i), Form("R50%% vs R95%% del Jet %d", i+1), 800, 600);
        gStyle->SetNumberContours(10); // Numero de colores a utilizar
        gPad->SetLogz(); // Escala logaritmica en el eje z
        hR50_vs_R95[i]->GetXaxis()->SetTitle("R50% del pT total");
        hR50_vs_R95[i]->GetYaxis()->SetTitle("R95% del pT total");
        hR50_vs_R95[i]->Draw("COLZ");
        cR50_vs_R95->SaveAs(Form("plots/R50_vs_R95_Jet%d.png", i+1));
        hR50_vs_R95[i]->Draw("CONT1");
        cR50_vs_R95->SaveAs(Form("plots/CONT_R50_vs_R95_Jet%d.png", i+1));
        delete cR50_vs_R95;
    }

    std::cout << "Los histogramas se han dibujado y guardado correctamente." << std::endl;
}


void JetAnalyzer::SaveHistograms(const std::string& outputDir) {
    // Crear directorio de salida si no existe
    std::string command = "mkdir -p " + outputDir;
    system(command.c_str());

    // Guarda los histogramas en un archivo ROOT
    TFile outFile((outputDir + "/histograms.root").c_str(), "RECREATE");
    hJetsPerEvent->Write();

    for (int i = 0; i < 4; i++) {
        hJetPT[i]->Write();
        hJetEta[i]->Write();
        hJetPhi[i]->Write();
        hChargedParticles[i]->Write();
        hNeutralsParticles[i]->Write();
        hChargedPTFraction[i]->Write();
        hNeutralPTFraction[i]->Write();
        hAveragePT[i]->Write();
        hParticlesBelowAvgPT[i]->Write();
        hParticlesAboveAvgPT[i]->Write();
        hMaxPTRatio[i]->Write();
        hMinPTRatio[i]->Write();
        hMaxDRRatio[i]->Write();
        hMinDRRatio[i]->Write();
        hDeltaRMaxPT[i]->Write();
        hDeltaRMinPT[i]->Write();
        hDeltaRMaxDR[i]->Write();
        hDeltaRMinDR[i]->Write();
        hPTDifference[i]->Write();
        hR50PercentPT[i]->Write();
        hR95PercentPT[i]->Write();
        hCumulativePT_vs_DeltaR[i]->Write();
        hDeltaR_vs_D0[i]->Write();
        hCumulativePT_vs_D0[i]->Write();
        hPT_vs_Eta[i]->Write();
        hCharged_vs_NeutralParticles[i]->Write();
        hChargedPTFraction_vs_NeutralPTFraction[i]->Write();
        hAveragePT_vs_TotalParticles[i]->Write();
        hMaxPTRatio_vs_DeltaRMaxPT[i]->Write();
        hR50_vs_R95[i]->Write();
    }

    for (int i = 0; i < 6; i++) {
        hDeltaRPar[i]->Write();
        //hDeltaR_vs_PTDifference[i]->Write();
    }

    outFile.Close();
}
