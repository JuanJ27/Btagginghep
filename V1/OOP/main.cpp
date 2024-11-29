#include "JetAnalyzer.cpp"

int main() {
    // Lista de archivos de entrada
    std::vector<std::string> inputFiles = {
        "/home/juan27/Btagginghep/rootfiles/run_01.root",
        "/home/juan27/Btagginghep/rootfiles/run_02.root",
        "/home/juan27/Btagginghep/rootfiles/run_03.root",
        "/home/juan27/Btagginghep/rootfiles/run_04.root"
    };

    // Crear instancia de JetAnalyzer
    JetAnalyzer analyzer(inputFiles);

    // Procesar eventos
    analyzer.LoopEvents();

    // Dibujar histogramas
    analyzer.DrawHistograms();

    // Guardar histogramas
    analyzer.SaveHistograms("plots");

    std::cout << "El análisis ha finalizado correctamente." << std::endl;

    return 0;
}