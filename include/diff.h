#ifndef DIFF_H
#define DIFF_H

float forward(float x, float h, int order);

float backward(float x, float h, int order);

float central(float x, float h, int order);


#endif