#pragma once

// =============================================================================
// test_framework.hpp
// Infraestructura compartida del framework de tests unitarios.
//
// DISEÑO DEL MACRO REGISTER_TEST:
// Usar lambdas resuelve el problema de las comas en el cuerpo del test.
// La lambda se crea en el punto de registro y se almacena en el vector.
// El cuerpo del test va entre corchetes normales después del nombre:
//
//   REGISTER_TEST(mi_test) {
//       ASSERT_NEAR(1.0, 1.0, 1e-10);
//       return true;
//   }
//
// Se usa un constructor estático para registrar antes de main().
// =============================================================================

#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

// ── Estructura de un caso de test ──────────────────────────────────────────
struct TestCase {
    std::string            name;
    std::function<bool()>  fn;
};

// ── Registro global (singleton) ────────────────────────────────────────────
inline std::vector<TestCase>& test_registry() {
    static std::vector<TestCase> reg;
    return reg;
}

// ── Helper de registro ─────────────────────────────────────────────────────
inline void register_test(const std::string& name, std::function<bool()> fn) {
    TestCase tc;
    tc.name = name;
    tc.fn   = fn;
    test_registry().push_back(tc);
}

// ── Macro de registro ──────────────────────────────────────────────────────
// Declara una función normal (no lambda) para evitar todos los problemas
// con comas y macros. El cuerpo sigue al macro en el fuente.
//
//   REGISTER_TEST(nombre_test) { ... código ... return true; }
//
#define REGISTER_TEST(test_name)                                          \
    static bool _testbody_##test_name();                                  \
    static struct _AutoReg_##test_name {                                  \
        _AutoReg_##test_name() {                                          \
            register_test(#test_name, _testbody_##test_name);             \
        }                                                                  \
    } _autoreg_instance_##test_name;                                      \
    static bool _testbody_##test_name()

// ── Macros de aserción ──────────────────────────────────────────────────────

#define ASSERT_TRUE(expr) \
    if (!(expr)) { \
        std::cerr << "    FALLO linea " << __LINE__ \
                  << ": (" << #expr << ") es falso\n"; \
        return false; \
    }

#define ASSERT_FALSE(expr) \
    if ((expr)) { \
        std::cerr << "    FALLO linea " << __LINE__ \
                  << ": (" << #expr << ") es verdadero\n"; \
        return false; \
    }

#define ASSERT_NEAR(a, b, tol) \
    do { \
        double _a = (double)(a), _b = (double)(b), _t = (double)(tol); \
        if (std::abs(_a - _b) > _t) { \
            std::cerr << "    FALLO linea " << __LINE__ << ": " \
                      << "|" << #a << " - " << #b << "| = " \
                      << std::abs(_a - _b) << " > " << _t << "\n"; \
            return false; \
        } \
    } while(0)

#define ASSERT_EQ(a, b) \
    if ((a) != (b)) { \
        std::cerr << "    FALLO linea " << __LINE__ << ": " \
                  << #a << " != " << #b << "\n"; \
        return false; \
    }

#define ASSERT_GT(a, b) \
    if (!((a) > (b))) { \
        std::cerr << "    FALLO linea " << __LINE__ << ": " \
                  << "(" << #a << "=" << (a) << ") no > " << (b) << "\n"; \
        return false; \
    }

#define ASSERT_LT(a, b) \
    if (!((a) < (b))) { \
        std::cerr << "    FALLO linea " << __LINE__ << ": " \
                  << "(" << #a << "=" << (a) << ") no < " << (b) << "\n"; \
        return false; \
    }
