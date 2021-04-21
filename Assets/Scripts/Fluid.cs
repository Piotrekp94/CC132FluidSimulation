using UnityEngine;

public class Fluid : MonoBehaviour
{
    private static int size;
    [SerializeField] private float dt;
    [SerializeField] private float diff;
    [SerializeField] private float visc;
    [SerializeField] private GameObject gameObject;
    private float[] density;

    private GameObject[] objects;

    private float[] s;

    private float[] Vx;

    private float[] Vxo;
    private float[] Vy;
    private float[] Vyo;
    private float[] Vz;
    private float[] Vzo;

    private void Start()
    {
        size = 16;
        objects = new GameObject[size * size * size];

        s = new float[size * size * size];
        density = new float[size * size * size];
        Vx = new float[size * size * size];
        Vy = new float[size * size * size];
        Vz = new float[size * size * size];
        Vxo = new float[size * size * size];
        Vyo = new float[size * size * size];
        Vzo = new float[size * size * size];

        GenerateCubes(size);
    }

    private void Update()
    {
        var N = size;
        var visc = this.visc;
        var diff = this.diff;
        var dt = this.dt;
        var Vx = this.Vx;
        var Vy = this.Vy;
        var Vz = this.Vz;
        var Vx0 = Vxo;
        var Vy0 = Vyo;
        var Vz0 = Vzo;
        var s = this.s;
        var density = this.density;

        diffuse(1, Vx0, Vx, visc, dt, 2, N);
        diffuse(2, Vy0, Vy, visc, dt, 2, N);
        diffuse(3, Vz0, Vz, visc, dt, 2, N);

        project(Vx0, Vy0, Vz0, Vx, Vy, 2, N);

        advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
        advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
        advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);

        project(Vx, Vy, Vz, Vx0, Vy0, 2, N);

        diffuse(0, s, density, diff, dt, 2, N);
        advect(0, density, s, Vx, Vy, Vz, dt, N);


        AddDensity(8, 8, 8, 3);
        AddVelocity(8, 8, 8, 1, 1, 0);
        Render();
    }

    private void GenerateCubes(int i)
    {
        for (var a = 0; a < i; a++)
        for (var b = 0; b < i; b++)
        for (var c = 0; c < i; c++)
            objects[IX(a, b, c)] = Instantiate(gameObject, new Vector3(a, b, c), Quaternion.identity, transform);
    }

    private void Render()
    {
        for (var a = 0; a < size; a++)
        for (var b = 0; b < size; b++)
        for (var c = 0; c < size; c++)
        {
            var adensity = density[IX(a, b, c)];
            var x = objects[IX(a, b, c)];

            if (adensity <= 0.1)
                x.SetActive(false);
            else
                x.SetActive(true);


            x.GetComponent<MeshRenderer>().material.color = new Color(1.0f, 0.0f, 0.0f, density[IX(a, b, c)]);
        }
    }

    private void AddDensity(int x, int y, int z, float amount)
    {
        var index = IX(x, y, z);
        density[index] += amount;
    }

    private void AddVelocity(int x, int y, int z, float amountX, float amountY, float amountZ)
    {
        var index = IX(x, y, z);
        Vx[index] += amountX;
        Vy[index] += amountY;
        Vz[index] += amountZ;
    }

    private static int IX(int x, int y, int z)
    {
        return x + y * size + z * size * size;
    }

    private static void diffuse(int b, float[] x, float[] x0, float diff, float dt, int iter, int N)
    {
        var a = dt * diff * (N - 2) * (N - 2);
        lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
    }

    private static void lin_solve(int b, float[] x, float[] x0, float a, float c, int iter, int N)
    {
        var cRecip = (float) (1.0 / c);
        for (var k = 0; k < iter; k++)
        {
            for (var m = 1; m < N - 1; m++)
            for (var j = 1; j < N - 1; j++)
            for (var i = 1; i < N - 1; i++)
                x[IX(i, j, m)] =
                    (x0[IX(i, j, m)]
                     + a * (x[IX(i + 1, j, m)]
                            + x[IX(i - 1, j, m)]
                            + x[IX(i, j + 1, m)]
                            + x[IX(i, j - 1, m)]
                            + x[IX(i, j, m + 1)]
                            + x[IX(i, j, m - 1)]
                     )) * cRecip;
            set_bnd(b, x, N);
        }
    }

    private static void set_bnd(int b, float[] x, int N)
    {
        for (var j = 1; j < N - 1; j++)
        for (var i = 1; i < N - 1; i++)
        {
            x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
            x[IX(i, j, N - 1)] = b == 3 ? -x[IX(i, j, N - 2)] : x[IX(i, j, N - 2)];
        }

        for (var k = 1; k < N - 1; k++)
        for (var i = 1; i < N - 1; i++)
        {
            x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
            x[IX(i, N - 1, k)] = b == 2 ? -x[IX(i, N - 2, k)] : x[IX(i, N - 2, k)];
        }

        for (var k = 1; k < N - 1; k++)
        for (var j = 1; j < N - 1; j++)
        {
            x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
            x[IX(N - 1, j, k)] = b == 1 ? -x[IX(N - 2, j, k)] : x[IX(N - 2, j, k)];
        }

        x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)]
                                  + x[IX(0, 1, 0)]
                                  + x[IX(0, 0, 1)]);
        x[IX(0, N - 1, 0)] = 0.33f * (x[IX(1, N - 1, 0)]
                                      + x[IX(0, N - 2, 0)]
                                      + x[IX(0, N - 1, 1)]);
        x[IX(0, 0, N - 1)] = 0.33f * (x[IX(1, 0, N - 1)]
                                      + x[IX(0, 1, N - 1)]
                                      + x[IX(0, 0, N - 1)]);
        x[IX(0, N - 1, N - 1)] = 0.33f * (x[IX(1, N - 1, N - 1)]
                                          + x[IX(0, N - 2, N - 1)]
                                          + x[IX(0, N - 1, N - 2)]);
        x[IX(N - 1, 0, 0)] = 0.33f * (x[IX(N - 2, 0, 0)]
                                      + x[IX(N - 1, 1, 0)]
                                      + x[IX(N - 1, 0, 1)]);
        x[IX(N - 1, N - 1, 0)] = 0.33f * (x[IX(N - 2, N - 1, 0)]
                                          + x[IX(N - 1, N - 2, 0)]
                                          + x[IX(N - 1, N - 1, 1)]);
        x[IX(N - 1, 0, N - 1)] = 0.33f * (x[IX(N - 2, 0, N - 1)]
                                          + x[IX(N - 1, 1, N - 1)]
                                          + x[IX(N - 1, 0, N - 2)]);
        x[IX(N - 1, N - 1, N - 1)] = 0.33f * (x[IX(N - 2, N - 1, N - 1)]
                                              + x[IX(N - 1, N - 2, N - 1)]
                                              + x[IX(N - 1, N - 1, N - 2)]);
    }

    private static void project(float[] velocX, float[] velocY, float[] velocZ, float[] p, float[] div, int iter, int N)
    {
        for (var k = 1; k < N - 1; k++)
        for (var j = 1; j < N - 1; j++)
        for (var i = 1; i < N - 1; i++)
        {
            div[IX(i, j, k)] = -0.5f * (
                velocX[IX(i + 1, j, k)]
                - velocX[IX(i - 1, j, k)]
                + velocY[IX(i, j + 1, k)]
                - velocY[IX(i, j - 1, k)]
                + velocZ[IX(i, j, k + 1)]
                - velocZ[IX(i, j, k - 1)]
            ) / N;
            p[IX(i, j, k)] = 0;
        }

        set_bnd(0, div, N);
        set_bnd(0, p, N);
        lin_solve(0, p, div, 1, 6, iter, N);

        for (var k = 1; k < N - 1; k++)
        for (var j = 1; j < N - 1; j++)
        for (var i = 1; i < N - 1; i++)
        {
            velocX[IX(i, j, k)] -= 0.5f * (p[IX(i + 1, j, k)]
                                           - p[IX(i - 1, j, k)]) * N;
            velocY[IX(i, j, k)] -= 0.5f * (p[IX(i, j + 1, k)]
                                           - p[IX(i, j - 1, k)]) * N;
            velocZ[IX(i, j, k)] -= 0.5f * (p[IX(i, j, k + 1)]
                                           - p[IX(i, j, k - 1)]) * N;
        }

        set_bnd(1, velocX, N);
        set_bnd(2, velocY, N);
        set_bnd(3, velocZ, N);
    }

    private static void advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float[] velocZ, float dt,
        int N)
    {
        float i0, i1, j0, j1, k0, k1;

        var dtx = dt * (N - 2);
        var dty = dt * (N - 2);
        var dtz = dt * (N - 2);

        float s0, s1, t0, t1, u0, u1;
        float tmp1, tmp2, tmp3, x, y, z;

        float Nfloat = N;
        float ifloat, jfloat, kfloat;
        int i, j, k;

        for (k = 1, kfloat = 1; k < N - 1; k++, kfloat++)
        for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++)
        for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++)
        {
            tmp1 = dtx * velocX[IX(i, j, k)];
            tmp2 = dty * velocY[IX(i, j, k)];
            tmp3 = dtz * velocZ[IX(i, j, k)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;
            z = kfloat - tmp3;

            if (x < 0.5f) x = 0.5f;
            if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
            i0 = Mathf.Floor(x);
            i1 = i0 + 1.0f;
            if (y < 0.5f) y = 0.5f;
            if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
            j0 = Mathf.Floor(x);
            j1 = j0 + 1.0f;
            if (z < 0.5f) z = 0.5f;
            if (z > Nfloat + 0.5f) z = Nfloat + 0.5f;
            k0 = Mathf.Floor(x);
            k1 = k0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;
            u1 = z - k0;
            u0 = 1.0f - u1;

            var i0i = (int) i0;
            var i1i = (int) i1;
            var j0i = (int) j0;
            var j1i = (int) j1;
            var k0i = (int) k0;
            var k1i = (int) k1;

            d[IX(i, j, k)] =
                s0 * (t0 * (u0 * d0[IX(i0i, j0i, k0i)]
                            + u1 * d0[IX(i0i, j0i, k1i)])
                      + t1 * (u0 * d0[IX(i0i, j1i, k0i)]
                              + u1 * d0[IX(i0i, j1i, k1i)]))
                + s1 * (t0 * (u0 * d0[IX(i1i, j0i, k0i)]
                              + u1 * d0[IX(i1i, j0i, k1i)])
                        + t1 * (u0 * d0[IX(i1i, j1i, k0i)]
                                + u1 * d0[IX(i1i, j1i, k1i)]));
        }

        set_bnd(b, d, N);
    }
}