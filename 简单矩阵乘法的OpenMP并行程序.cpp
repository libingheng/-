// 李秉恒 2019152001 并行计算1
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
using namespace std;
const int n = 1000;
double a[n][n], b[n][n], c[n][n], d[n][n];
int t_num = 4; // thread num
int main()
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i][j] = rand() * 1.0 / RAND_MAX;
            b[i][j] = rand() * 1.0 / RAND_MAX;
            c[i][j] = 0;
            d[i][j] = 0;
        }
    }
    clock_t start = clock();
    // 串行执行:
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                d[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    clock_t end = clock();
    double duration = (end - start);
    cout << "串行计算耗时: " << duration << "ms" << endl;
    start = clock();

    // 并行计算执行:
    omp_set_num_threads(t_num);
#pragma omp parallel shared(a, b, d)
    {
#pragma omp for
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
    }
    end = clock();
    duration = (end - start);
    cout << t_num << " 线程并行计算耗时: " << duration << "ms" << endl;

    bool error = false;
    for (int i = 0; i < n && !error; i++)
    {
        for (int j = 0; j < n && !error; j++)
        {
            if (c[i][j] != d[i][j])
            {
                error = true;
                break;
            }
        }
    }
    if (!error)
    {
        cout << "correct" << endl;
    }
}
