#include <iostream>
#include "include/function.h"
#include "include/diff.h"

using namespace std;

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
        result = forward(x, h, order);
        break;
    case 1:
        result = forward(x, h, order);
        break;
    case 2:
        result = forward(x, h, order);
        break;
    default:
        cout << "Inválido.";
        break;
    }

    cout << "resultado: " << result << endl;

    return 0;
}