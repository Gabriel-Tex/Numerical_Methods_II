cmake_minimum_required(VERSION 3.28)
project(NumericalMethods)

add_executable(methods
    main.cpp
    src/function.cpp
    src/numericalDiff/finiteDiff.cpp
    src/numericalDiff/finiteDiffTaylor.cpp
    src/numericalDiff/finiteDiffNewton.cpp
    src/numericalDiff/menuDerivacao.cpp
    src/numericalInt/newtonCotes.cpp
    src/numericalInt/gaussLegendre.cpp
    src/numericalInt/gaussSpecial.cpp
    src/numericalInt/singularIntegral.cpp
    src/numericalInt/menuIntegracao.cpp
    src/numericalInt/multiDimIntegral.cpp
    src/eigenValues/matrizUtils.cpp
    src/eigenValues/potencia.cpp
    src/eigenValues/householder.cpp
    src/eigenValues/jacobi.cpp
    src/eigenValues/metodoQR.cpp
    src/eigenValues/menuAutovalores.cpp
)

target_include_directories(methods PRIVATE include)