#include "../include/function.h"

float forward(float x, float h, int order){
    if(order == 0) return F(x);

    float left = forward(x+h, h, order - 1);
    float right = forward(x, h, order - 1);

    return ( left - right ) / h;
}

float backward(float x, float h, int order){
    if(order == 0) return F(x);

    float left = backward(x, h, order - 1);
    float right = backward(x-h, h, order - 1);

    return ( left - right ) / h;

}

float central(float x, float h, int order){
    if (order == 0) return F(x);

    float left  = central(x+h, h, order - 1);
    float right = central(x-h, h, order - 1);

    return (left - right) / (2*h);
}