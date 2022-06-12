#include <vector>
#include <math.h>
#include <iostream>
using namespace std;
typedef vector<vector<float> > matrix;

int main()
{
    const int nx = 41;
    const int ny = 41;
    const int nt = 500;
    const int nit = 50;
    float dx = 2 / (nx - 1);
    float dy = 3 / (ny - 1);
    const float dt = 0.01;
    const int rho = 1;
    const float nu = 0.02;

    matrix u(ny, vector<float>(nx));
    matrix v(ny, vector<float>(nx));
    matrix p(ny, vector<float>(nx));
    matrix b(ny, vector<float>(nx));

    matrix un(ny, vector<float>(nx));
    matrix vn(ny, vector<float>(nx));
    matrix pn(ny, vector<float>(nx));

    int nti, nyi, nxi, niti;
    // np.zeros
    for (nyi = 0; nyi < ny; nyi++)
    {
        for (nxi = 0; nxi < nx; nxi++)
        {
            u[nyi][nxi] = v[nyi][nxi] = p[nyi][nxi] = b[nyi][nxi] = 0.0;
        }
    }

    for (nti = 0; nti < nt; nti++)
    {
        for (nyi = 1; nyi < ny - 1; nyi++)
        {
            for (nxi = 1; nxi < nx - 1; nxi++)
            {
                b[nyi][nxi] = rho * (1 / dt *
                                         ((u[nyi][nxi + 1] - u[nyi][nxi - 1]) / (2 * dx) + (v[nyi + 1][nxi] - v[nyi - 1][nxi]) / (2 * dy)) -
                                     pow((u[nyi][nxi + 1] - u[nyi][nxi - 1]) / (2 * dx), 2.0) - 2 * ((u[nyi + 1][nxi] - u[nyi - 1][nxi]) / (2 * dy) * (v[nyi][nxi + 1] - v[nyi][nxi - 1]) / (2 * dx)) - pow((v[nyi + 1][nxi] - v[nyi - 1][nxi]) / (2 * dy), 2.0));
            }
        }
        // std::cout << "b\n";

        for (niti = 0; niti < nit; niti++)
        {
            // copy p to pn
            for (nyi = 0; nyi < ny; nyi++)
            {
                for (nxi = 0; nxi < nx; nxi++)
                {
                    pn[nyi][nxi] = p[nyi][nxi];
                }
            }

            for (nyi = 1; nyi < ny-1; nyi++)
            {
                for (nxi = 1; nxi < nx-1; nxi++)
                {
                    p[nyi][nxi] = (pow(dy, 2.0) * (pn[nyi][nxi + 1] + pn[nyi][nxi - 1]) +
                                   pow(dx, 2.0) * (pn[nyi + 1][nxi] + pn[nyi - 1][nxi]) -
                                   b[nyi][nxi] * pow(dx, 2.0) * pow(dy, 2.0)) /
                                  (2 * (pow(dx, 2.0) + pow(dy, 2.0)));
                }
            }

            // std::cout << "p\n";

            // p matrix slicing of x axis
            for (nyi = 0; nyi < ny; nyi++)
            {
                p[nyi][nx - 1] = p[nyi][nx - 2];
                p[nyi][0] = p[nyi][1];
            }

            // p matrix slicing of y axis
            for (nxi = 0; nxi < nx; nxi++)
            {
                p[0][nxi] = p[1][nxi];
                p[ny - 1][nxi] = 0.0;
            }
        }

        // copy u, v to un, vn
        for (nyi = 0; nyi < ny; nyi++)
        {
            for (nxi = 0; nxi < nx; nxi++)
            {
                un[nyi][nxi] = u[nyi][nxi];
                vn[nyi][nxi] = v[nyi][nxi];
            }
        }

        for (nyi = 1; nyi < ny - 1; nyi++)
        {
            for (nxi = 1; nxi < nx - 1; nxi++)
            {
                u[nyi][nxi] = un[nyi][nxi] - un[nyi][nxi] * dt / dx * (un[nyi][nxi] - un[nyi][nxi - 1]) - un[nyi][nxi] * dt / dy * (un[nyi][nxi] - un[nyi - 1][nxi]) - dt / (2 * rho * dx) * (p[nyi][nxi + 1] - p[nyi][nxi - 1]) + nu * dt / pow(dx, 2.0) * (un[nyi][nxi + 1] - 2 * un[nyi][nxi] + un[nyi][nxi - 1]) + nu * dt / pow(dy, 2.0) * (un[nyi + 1][nxi] - 2 * un[nyi][nxi] + un[nyi - 1][nxi]);
                v[nyi][nxi] = vn[nyi][nxi] - vn[nyi][nxi] * dt / dx * (vn[nyi][nxi] - vn[nyi][nxi - 1]) - vn[nyi][nxi] * dt / dy * (vn[nyi][nxi] - vn[nyi - 1][nxi]) - dt / (2 * rho * dx) * (p[nyi + 1][nxi] - p[nyi - 1][nxi]) + nu * dt / pow(dx, 2.0) * (vn[nyi][nxi + 1] - 2 * vn[nyi][nxi] + vn[nyi][nxi - 1]) + nu * dt / pow(dy, 2.0) * (vn[nyi + 1][nxi] - 2 * vn[nyi][nxi] + vn[nyi - 1][nxi]);
            }
        }

        // set zeros and ones to u, v
        for (nyi = 0; nyi < ny; nyi++)
        {
            u[nyi][0] = u[nyi][nx - 1] = v[nyi][0] = v[nyi][nx - 1] = 0.0;
        }

        for (nxi = 0; nxi < nx; nxi++)
        {
            u[0][nxi] = v[0][nxi] = v[ny - 1][nxi] = 0.0;
            u[ny - 1][0] = 1.0;
        }
    }

    return 0;
}