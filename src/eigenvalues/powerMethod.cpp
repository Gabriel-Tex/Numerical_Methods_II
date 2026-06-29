#include "../../include/eigenvalues/powerMethod.h"
#include <cmath>

void regularPower(const Matrix &A, const Vector &v0, double eps,
                  double &lambda, Vector &x)
{
    int n = (int)A.size();
    double lambdaNew = 0.0;
    Vector vNew = v0;
    double lambdaOld;
    Vector xOld(n);

    do
    {
        lambdaOld = lambdaNew;
        Vector vOld = vNew;
        xOld = normalizeVec(vOld);
        vNew = matVecProduct(A, xOld);
        lambdaNew = dotProduct(xOld, vNew);
    } while (std::fabs((lambdaNew - lambdaOld) / lambdaNew) > eps);

    lambda = lambdaNew;
    x = xOld;
}

void inversePower(const Matrix &A, const Vector &v0, double eps,
                  double &lambda, Vector &x)
{
    int n = (int)A.size();
    Matrix L, U;
    luDecomposition(A, L, U);

    double lambdaBarNew = 0.0;
    Vector vNew = v0;
    double lambdaBarOld;
    Vector xOld(n);

    do
    {
        lambdaBarOld = lambdaBarNew;
        Vector vOld = vNew;
        xOld = normalizeVec(vOld);
        vNew = solveLU(L, U, xOld);
        lambdaBarNew = dotProduct(xOld, vNew);
    } while (std::fabs((lambdaBarNew - lambdaBarOld) / lambdaBarNew) > eps);

    lambda = 1.0 / lambdaBarNew;
    x = xOld;
}

void shiftedPower(const Matrix &A, const Vector &v0, double eps, double mu,
                  double &lambda, Vector &x)
{
    int n = (int)A.size();
    Matrix Ahat = A;
    for (int i = 0; i < n; i++)
        Ahat[i][i] -= mu;

    double lambdaHat;
    Vector xHat;
    inversePower(Ahat, v0, eps, lambdaHat, xHat);

    lambda = lambdaHat + mu;
    x = xHat;
}