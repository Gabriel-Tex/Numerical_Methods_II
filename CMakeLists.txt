cmake_minimum_required(VERSION 3.28)
project(NumericalMethods)

add_executable(methods
    main.cpp
    src/function.cpp
    src/numericalDiff/finiteDiff.cpp
    src/numericalDiff/finiteDiffTaylor.cpp
    src/numericalDiff/finiteDiffNewton.cpp
    src/numericalDiff/menuNumericalDiff.cpp
    src/numericalInt/newtonCotes.cpp
    src/numericalInt/gaussLegendre.cpp
    src/numericalInt/gaussSpecial.cpp
    src/numericalInt/singularIntegral.cpp
    src/numericalInt/multiDimIntegral.cpp
    src/numericalInt/menuNumericalInt.cpp
    src/eigenvalues/matrixUtils.cpp
    src/eigenvalues/powerMethod.cpp
    src/eigenvalues/householder.cpp
    src/eigenvalues/jacobiMethod.cpp
    src/eigenvalues/qrMethod.cpp
    src/eigenvalues/menuEigenvalues.cpp
    src/pvi/pviUtils.cpp
    src/pvi/explicitEuler.cpp
    src/pvi/implicitEuler.cpp
    src/pvi/rungeKutta.cpp
    src/pvi/predictorCorrector.cpp
    src/pvi/menuPVI.cpp
    src/pvc/pvcUtils.cpp
    src/pvc/finiteDiff.cpp
    src/pvc/finiteElements.cpp
    src/pvc/menuPVC.cpp
)

target_include_directories(methods PRIVATE include)