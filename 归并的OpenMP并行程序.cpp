#include <omp.h>
//#include<stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>
#include <fstream>
#include <Windows.h>
using namespace std;
#define DATASET_SIZE 100000000

// constexpr auto DATASET_SIZE = 10;
#define THREAD_NUM 1

void print_vector(vector<double> &x)
{
    for (int i = 0; i < x.size(); i++)

    {
        std::cout << x[i] << " ";
    }
    std::cout << endl
              << endl;
}

bool compare(vector<double> &a, vector<double> &b)
{
    for (int i = 0; i < a.size(); i++)
    {
        if (a[i] != b[i])
            return false;
    }

    return true;
}

int cal_rank(vector<double> &a, double x)
{
    bool flag;
    long result;
    flag = false;
    for (long i = 0; i < DATASET_SIZE; i++)
    {
        if (flag)
            continue;
        if (a[i] < x && x < a[i + 1])
        {
            result = i + 1;
            cout << "result = " << result << endl;
            flag = true;
        }
    }

    return result;
}

int cal_rank_se(double a[], double x)
{
    int u = 0;
    int v = DATASET_SIZE + 1;
    int m;
    while (u != v)
    {
        m = (u + v) / 2;
        if (a[m] < x)
            u = m + 1;
        else
        {
            v = m;
        }
    }
    return m + 1;
}

bool cheak_vector(double x[], double y[])
{

    for (int i = 0; i < DATASET_SIZE * 2; i++)
    {
        if (x[i] != y[i])
            return false;
    }
    return true;
}

int main()
{

    const char *logfilename = "log.txt";
    vector<double> second;
    // vector<double> a(DATASET_SIZE);
    // vector<double> b(DATASET_SIZE);

    ////vector<double> backupa, backupb;
    // vector<double> result_se(DATASET_SIZE * 2);
    // vector<double> result_pa(DATASET_SIZE * 2);
    double *a = new double[DATASET_SIZE];
    double *b = new double[DATASET_SIZE];

    double *result_se = new double[DATASET_SIZE * 2];
    double *result_pa = new double[DATASET_SIZE * 2];

    // acopy.resize(DATASET_SIZE);
    for (long i = 0; i < DATASET_SIZE; i++)
    {
        a[i] = rand() / double(RAND_MAX);
        b[i] = rand() / double(RAND_MAX);
        // acopy[i] = a[i];
        // backupa[i] = a[i];
        // backupb[i] = b[i];
    }

    // 串行计算
    clock_t start, end;
    // 开始计时

    // 开始排序
    sort(a, a + DATASET_SIZE);
    sort(b, b + DATASET_SIZE);

    start = clock();
    cout << "[SERIAL] start sorting " << endl
         << endl;

    merge(a, a + DATASET_SIZE, b, b + DATASET_SIZE, result_se);

    // print_vector(result_se);

    end = clock();
    double duration = double(end - start) * 1000 / CLOCKS_PER_SEC;
    double chuan = duration;

    cout << " 串行用时 = " << chuan << " ms" << endl
         << endl;
    ofstream log;
    log.open(logfilename, ios::out | ios::app);
    log << "串行运算用时： " << chuan << " ms" << endl
        << endl;

    //并行计算
    // 开始计时
    start = clock();
    cout << "[PARALLEL] start sorting " << endl
         << endl;

    // create u, v
    vector<int> u(THREAD_NUM + 1);
    vector<int> v(THREAD_NUM + 1);
    u[0] = v[0] = 0;
    u[THREAD_NUM] = v[THREAD_NUM] = DATASET_SIZE;

    omp_set_num_threads(THREAD_NUM);
    // 并行计算 for循环 求u，v
#pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < THREAD_NUM; i++)
    {
        // cout << omp_get_thread_num() << endl;

        u[i] = i * DATASET_SIZE / THREAD_NUM;
        v[i] = cal_rank_se(b, a[u[i] - 1]);
    }

    // 并行计算for循环 合并
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < THREAD_NUM; i++)
    {
        //
        merge(a + u[i], a + u[i + 1], b + v[i], b + v[i + 1], result_pa + u[i] + v[i]);
    }
    // print_vector(a);
    // print_vector(b);

    // print_vector(result_pa);

    end = clock();

    duration = double(end - start) * 1000 / CLOCKS_PER_SEC;

    cout << " 并行用时 = " << duration << " ms" << endl;

    // log.open(logfilename, ios::out | ios::app);
    log << "并行运算用时： " << duration << " ms" << endl
        << endl;

    if (cheak_vector(result_se, result_pa))
    {
        cout << "结果正确" << endl;
    }

    delete[] a, b, result_se, result_pa;
    // int total = 0;
    // for (int k = 0; k < 5; k++)
    //{
    //	total += second[k];
    // }
    // total = total / 5;
    // double rate = chuan / total;
    // cout << "平均值 " << total << " ms" << endl;
    // log << "平均值 " << total << " ms" << endl;
    // cout << "加速比 " << rate << endl;
    // log << "加速比 " << rate << endl;
}
