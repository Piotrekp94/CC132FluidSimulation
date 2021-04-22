using System;
using UnityEngine;
using UnityEngine.UIElements;
using Random = UnityEngine.Random;

public class Fluid : MonoBehaviour
{
    private static int _size;
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
        _size = 16;
        var arraySize = _size * _size * _size;
        objects = new GameObject[arraySize];

        s = new float[arraySize];
        density = new float[arraySize];
        Vx = new float[arraySize];
        Vy = new float[arraySize];
        Vz = new float[arraySize];
        Vxo = new float[arraySize];
        Vyo = new float[arraySize];
        Vzo = new float[arraySize];

        GenerateCubes(_size);
    }

    private void Update()
    {
        const int iter = 2;

        Diffuse(1, Vxo, Vx, visc, dt, iter);
        Diffuse(2, Vyo, Vy, visc, dt, iter);
        Diffuse(3, Vzo, Vz, visc, dt, iter);

        project(Vxo, Vyo, Vzo, Vx, Vy, iter);
        advect(1, Vx, Vxo, Vxo, Vyo, Vzo, dt);
        advect(2, Vy, Vyo, Vxo, Vyo, Vzo, dt);
        advect(3, Vz, Vzo, Vxo, Vyo, Vzo, dt);

        project(Vx, Vy, Vz, Vxo, Vyo, iter);

        Diffuse(0, s, density, diff, dt, iter);
        advect(0, density, s, Vx, Vy, Vz, dt);
        int x = Random.Range(0, 16);
        int y = Random.Range(0, 16);
        int z = Random.Range(0, 16);
        AddDensity(x, y, z, 1000);
        AddVelocity(x, y, z, 10, 10, 10);
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
        for (var a = 0; a < _size; a++)
        for (var b = 0; b < _size; b++)
        for (var c = 0; c < _size; c++)
        {
            var adensity = density[IX(a, b, c)];
            var x = objects[IX(a, b, c)];

            // x.SetActive(!(adensity <= 0.01));

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
        return x + y * _size + z * _size * _size;
    }

    private static void Diffuse(int b, float[] x, float[] x0, float diff, float dt, int iter)
    {
        var a = dt * diff * (_size - 2) * (_size - 2);
        lin_solve(b, x, x0, a, 1 + 6 * a, iter);
    }

    private static void lin_solve(int b, float[] x, float[] x0, float a, float c, int iter)
    {
        var cRecip = (float) (1.0 / c);
        for (var k = 0; k < iter; k++)
        {
            for (var m = 1; m < _size - 1; m++)
            for (var j = 1; j < _size - 1; j++)
            for (var i = 1; i < _size - 1; i++)
                x[IX(i, j, m)] =
                    (x0[IX(i, j, m)]
                     + a * (x[IX(i + 1, j, m)]
                            + x[IX(i - 1, j, m)]
                            + x[IX(i, j + 1, m)]
                            + x[IX(i, j - 1, m)]
                            + x[IX(i, j, m + 1)]
                            + x[IX(i, j, m - 1)]
                     )) * cRecip;
            set_bnd(b, x);
        }
    }

    private static void set_bnd(int b, float[] x)
    {
        for (var j = 1; j < _size - 1; j++)
        for (var i = 1; i < _size - 1; i++)
        {
            x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
            x[IX(i, j, _size - 1)] = b == 3 ? -x[IX(i, j, _size - 2)] : x[IX(i, j, _size - 2)];
        }

        for (var k = 1; k < _size - 1; k++)
        for (var i = 1; i < _size - 1; i++)
        {
            x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
            x[IX(i, _size - 1, k)] = b == 2 ? -x[IX(i, _size - 2, k)] : x[IX(i, _size - 2, k)];
        }

        for (var k = 1; k < _size - 1; k++)
        for (var j = 1; j < _size - 1; j++)
        {
            x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
            x[IX(_size - 1, j, k)] = b == 1 ? -x[IX(_size - 2, j, k)] : x[IX(_size - 2, j, k)];
        }

        x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)]
                                  + x[IX(0, 1, 0)]
                                  + x[IX(0, 0, 1)]);
        x[IX(0, _size - 1, 0)] = 0.33f * (x[IX(1, _size - 1, 0)]
                                          + x[IX(0, _size - 2, 0)]
                                          + x[IX(0, _size - 1, 1)]);
        x[IX(0, 0, _size - 1)] = 0.33f * (x[IX(1, 0, _size - 1)]
                                          + x[IX(0, 1, _size - 1)]
                                          + x[IX(0, 0, _size - 1)]);
        x[IX(0, _size - 1, _size - 1)] = 0.33f * (x[IX(1, _size - 1, _size - 1)]
                                                  + x[IX(0, _size - 2, _size - 1)]
                                                  + x[IX(0, _size - 1, _size - 2)]);
        x[IX(_size - 1, 0, 0)] = 0.33f * (x[IX(_size - 2, 0, 0)]
                                          + x[IX(_size - 1, 1, 0)]
                                          + x[IX(_size - 1, 0, 1)]);
        x[IX(_size - 1, _size - 1, 0)] = 0.33f * (x[IX(_size - 2, _size - 1, 0)]
                                                  + x[IX(_size - 1, _size - 2, 0)]
                                                  + x[IX(_size - 1, _size - 1, 1)]);
        x[IX(_size - 1, 0, _size - 1)] = 0.33f * (x[IX(_size - 2, 0, _size - 1)]
                                                  + x[IX(_size - 1, 1, _size - 1)]
                                                  + x[IX(_size - 1, 0, _size - 2)]);
        x[IX(_size - 1, _size - 1, _size - 1)] = 0.33f * (x[IX(_size - 2, _size - 1, _size - 1)]
                                                          + x[IX(_size - 1, _size - 2, _size - 1)]
                                                          + x[IX(_size - 1, _size - 1, _size - 2)]);
    }

    private static void project(float[] velocX, float[] velocY, float[] velocZ, float[] p, float[] div, int iter)
    {
        for (var k = 1; k < _size - 1; k++)
        for (var j = 1; j < _size - 1; j++)
        for (var i = 1; i < _size - 1; i++)
        {
            div[IX(i, j, k)] = -0.5f * (
                velocX[IX(i + 1, j, k)]
                - velocX[IX(i - 1, j, k)]
                + velocY[IX(i, j + 1, k)]
                - velocY[IX(i, j - 1, k)]
                + velocZ[IX(i, j, k + 1)]
                - velocZ[IX(i, j, k - 1)]
            ) / _size;
            p[IX(i, j, k)] = 0;
        }

        set_bnd(0, div);
        set_bnd(0, p);
        lin_solve(0, p, div, 1, 6, iter);

        for (var k = 1; k < _size - 1; k++)
        for (var j = 1; j < _size - 1; j++)
        for (var i = 1; i < _size - 1; i++)
        {
            velocX[IX(i, j, k)] -= 0.5f * (p[IX(i + 1, j, k)]
                                           - p[IX(i - 1, j, k)]) * _size;
            velocY[IX(i, j, k)] -= 0.5f * (p[IX(i, j + 1, k)]
                                           - p[IX(i, j - 1, k)]) * _size;
            velocZ[IX(i, j, k)] -= 0.5f * (p[IX(i, j, k + 1)]
                                           - p[IX(i, j, k - 1)]) * _size;
        }

        set_bnd(1, velocX);
        set_bnd(2, velocY);
        set_bnd(3, velocZ);
    }

    private static void advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float[] velocZ, float dt)
    {
        float i0, i1, j0, j1, k0, k1;

        var dtx = dt * (_size - 2);
        var dty = dt * (_size - 2);
        var dtz = dt * (_size - 2);

        float s0, s1, t0, t1, u0, u1;
        float tmp1, tmp2, tmp3, x, y, z;

        float Nfloat = _size;
        float ifloat, jfloat, kfloat;
        int i, j, k;

        for (k = 1, kfloat = 1; k < _size - 1; k++, kfloat++)
        for (j = 1, jfloat = 1; j < _size - 1; j++, jfloat++)
        for (i = 1, ifloat = 1; i < _size - 1; i++, ifloat++)
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

        set_bnd(b, d);
    }
}