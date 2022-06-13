
__global__ void cavity(int nx, int ny, int nt, int nit, float dx, float dy, float dt, int rho, float nu, int N, float *u, float *v, float *p, float *b, float *un, float *vn, float *pn) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    bool top = (i < nx);
    bool btm = (i > (ny-1) * nx - 1);
    bool lft = (i % nx == 0);
    bool rgt = ((i + 1) % nx == 0);

    int niti;

    // not boundary values
    if (^(top || btm || lft || rgt)) 
    {
        b[i] = rho * (1 / dt *
                            ((u[i+1] - u[i-1]) / (2 * dx) + (v[i+nx] - v[i-nx]) / (2 * dy)) -
                     powf(((u[i+1] - u[i-1]) / (2 * dx)), 2.0) - 2 * ((u[i+nx] - u[i-nx]) / (2 * dy) * (v[i+1] - v[i-1]) / (2 * dx)) - powf(((v[i+nx] - v[i-nx]) / (2 * dy)), 2.0));
    }

    for (niti = 0; niti < nit; niti++) 
    {
        pn[i] = p[i];
    
        if (^(top || btm || lft || rgt)) 
        {
            p[i] = (powf(dy, 2.0) * (pn[i+1] + pn[i-1]) + powf(dx, 2.0) * (pn[i+nx] + pn[i-nx]) -
                    b[i] * powf(dx, 2.0) * powf(dy, 2.0)) / (2 * (powf(dx, 2.0) + powf(dy, 2.0)));
        }

        // p matrix slicing
        if (rgt) p[i] = p[i-1];
        if (lft) p[i] = p[i+1];
        if (top) p[i] = p[i+nx];
        if (btm) p[i] = 0.0;
    }

    // copy u, v to un, vn
    un[i] = u[i];
    vn[i] = v[i];

    // not boundary values
    if (^(top || btm || lft || rgt)) 
    {
        u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) - un[i] * dt / dy * (un[i] - un[i-nx])
                    - dt / (2 * rho * dx) * (p[i+1] - p[i-1]) + nu * dt / powf(dx, 2.0) * (un[i+1] - 2 * un[i] + un[i-1])
                    + nu * dt / powf(dy, 2.0) * (un[i+nx] - 2 * un[i] + un[i-nx]);
        v[i] = vn[i] - vn[i] * dt / dx * (vn[i] - vn[i - 1]) - vn[i] * dt / dy * (vn[i] - vn[i-nx])
                    - dt / (2 * rho * dx) * (p[i+nx] - p[i-nx]) + nu * dt / powf(dx, 2.0) * (vn[i+1] - 2 * vn[i] + vn[i-1])
                    + nu * dt / powf(dy, 2.0) * (vn[i+nx] - 2 * vn[i] + vn[i-nx]);
    }
    
    // set zeros and ones to u, v
    if (lft || rgt || top) u[i] = v[i] = 0.0;
    if (btm) 
    {
        u[i] = 1.0;
        v[i] = 0.0;
    }
}


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
    const int M = 16;
    const int N = nx * ny;

    float *u, *v, *p, *b, *un, *vn, *pn;
    int ni, nti;

    cudaMallocManaged(&u, ny * nx * sizeof(float));
    cudaMallocManaged(&v, ny * nx * sizeof(float));
    cudaMallocManaged(&p, ny * nx * sizeof(float));
    cudaMallocManaged(&b, ny * nx * sizeof(float));
    cudaMallocManaged(&un, ny * nx * sizeof(float));
    cudaMallocManaged(&vn, ny * nx * sizeof(float));
    cudaMallocManaged(&pn, ny * nx * sizeof(float));

    for (ni = 0; ni < ny * nx; ni++) u[ni] = v[ni] = p[ni] = b[ni] = 0.0;
    
    for (nti = 0; nti < nt; nti++) 
    {
        cavity<<<(N+M-1)/M,M>>>(nx, ny, nt, nit, dx, dy, dt, rho, nu, N, u, v, p, b, un, vn, pn);
        cudaDeviceSynchronize();
    }

    cudaFree(u);
    cudaFree(v);
    cudaFree(p);
    cudaFree(b);
    cudaFree(un);
    cudaFree(vn);
    cudaFree(pn);

    return 0;
}