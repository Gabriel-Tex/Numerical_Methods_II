#include <iostream>
#include "include/function.h"
#include "include/diffLimits.h"
#include "include/diffTaylor.h"

using namespace std;

/*
Anotações:

    para h = 1^-10, as derivadas primeiras convergem bem
    para h = 1^-6, as derivadas segundas convergem bem
    para h = 1^-1 e h=1^-2, as derivadas terceiras convergem bem
    para h = 1^-1 e h=1^-2, as derivadas quartas convergem bem
    quanto maior o h (para h < 1), mais a convergência melhora para as derivadas
terceiras e quartas (erro propagado e muitas subtrações)

*/


// main somente para testes
int main()
{
    // x, h, filosofia, Fx, ordem da derivada
    float x, h, result;
    int order, option;

    cout << "Entre o valor de x: ";
    cin >> x;

    cout << "Entre o valor de h: ";
    cin >> h;

    /*
    cout << "Entre a ordem da derivada: ";
    cin >> order;
    
    cout << "Escolha a filosofia: \n 0 - Fw \n 1 - Bw \n 2 - Ct \n";
    cin >> option;

    switch (option)
    {
    case 0:
        result = forwardRec(x, h, order);
        break;
    case 1:
        result = backwardRec(x, h, order);
        break;
    case 2:
        result = centralRec(x, h, order);
        break;
    }
    */
    
    cout << "=====================================" << endl;
    cout << "firstForward: " << firstForwardDerivation(x, h) << endl;
    cout << "firstCentral: " << firstCentralDerivation(x, h) << endl;
    cout << "firstBackward: " << firstBackwardDerivation(x, h) << endl;
    cout << "TaylorLinear: " << firstTaylorDerivation_LinearErrorAndForward(x, h).derivateValue << endl;
    cout << "TaylorQuadratic: " << firstTaylorDerivation_QuadraticErrorAndForward(x, h).derivateValue << endl;
    cout << "TaylorCubic: " << firstTaylorDerivation_CubicErrorAndForward(x, h).derivateValue << endl;
    
    cout << "=====================================" << endl;
    cout << "secondForward: " << secondForwardDerivation(x, h) << endl;
    cout << "secondCentral: " << secondCentralDerivation(x, h) << endl;
    cout << "secondBackward: " << secondBackwardDerivation(x, h) << endl;
    cout << "secondForwardRec: " << forwardRec(x, h, 2) << endl;
    cout << "secondCentralRec: " << backwardRec(x, h, 2) << endl;
    cout << "secondBackwardRec: " << centralRec(x, h, 2) << endl;

    cout << "=====================================" << endl;
    cout << "thirdForward: " << thirdForwardDerivation(x, h) << endl;
    cout << "thirdCentral: " << thirdCentralDerivation(x, h) << endl;
    cout << "thirdBackward: " << thirdBackwardDerivation(x, h) << endl;
    cout << "thirdForwardRec: " << forwardRec(x, h, 3) << endl;
    cout << "thirdCentralRec: " << backwardRec(x, h, 3) << endl;
    cout << "thirdBackwardRec: " << centralRec(x, h, 3) << endl;

    cout << "=====================================" << endl;
    cout << "fourthForward: " << fourthForwardDerivation(x, h) << endl;
    cout << "fourthCentral: " << fourthCentralDerivation(x, h) << endl;
    cout << "fourthBackward: " << fourthBackwardDerivation(x, h) << endl;
    cout << "fourthForwardRec: " << forwardRec(x, h, 4) << endl;
    cout << "fourthCentralRec: " << backwardRec(x, h, 4) << endl;
    cout << "fourthBackwardRec: " << centralRec(x, h, 4) << endl;
    
    return 0;
}