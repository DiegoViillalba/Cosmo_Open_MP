// =============================================================================
// test_main.cpp
// Punto de entrada de la suite de tests unitarios.
//
// Los casos de test se registran automáticamente mediante constructores
// estáticos en cada archivo test_*.cpp. Este main solo los ejecuta.
//
// Uso:
//   ./run_tests              → corre todos los tests
//   ./run_tests cic          → solo tests cuyo nombre contiene "cic"
//   ./run_tests integrador   → solo tests del integrador
// =============================================================================

#include "test_framework.hpp"
#include <iomanip>
#include <algorithm>
#include <stdexcept>

int main(int argc, char* argv[]) {

    std::string filter = (argc > 1) ? argv[1] : "";

    const auto& tests = test_registry();
    int passed = 0, failed = 0, skipped = 0;

    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════╗\n";
    std::cout << "║       Suite de tests — pm_cosmo  (Ng="
              << std::setw(3) << PM_GRID_SIZE << ")        ║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n\n";

    for (const auto& t : tests) {
        if (!filter.empty() && t.name.find(filter) == std::string::npos) {
            ++skipped;
            continue;
        }

        std::cout << "  [ RUN ] " << t.name << "\n";
        bool ok = false;
        try {
            ok = t.fn();
        } catch (const std::exception& e) {
            std::cerr << "    EXCEPCION: " << e.what() << "\n";
            ok = false;
        }

        if (ok) {
            std::cout << "  [ OK  ] " << t.name << "\n";
            ++passed;
        } else {
            std::cout << "  [FAIL] " << t.name << "\n";
            ++failed;
        }
    }

    std::cout << "\n"
              << "──────────────────────────────────────────\n"
              << "  Pasados : " << passed  << "\n"
              << "  Fallados: " << failed  << "\n"
              << "  Saltados: " << skipped << "\n"
              << "──────────────────────────────────────────\n";

    if (failed == 0)
        std::cout << "\n  TODOS LOS TESTS PASARON.\n\n";
    else
        std::cout << "\n  " << failed << " TEST(S) FALLARON.\n\n";

    return (failed == 0) ? 0 : 1;
}
