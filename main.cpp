#include <iostream>
#include "include/function.h"
#include "include/diffLimits.h"
#include "include/diffTaylor.h"

using namespace std;

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

    cout << "resultado: " << result << endl;

    return 0;
}