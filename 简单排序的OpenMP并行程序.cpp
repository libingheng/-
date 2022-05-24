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
#define THREAD_NUM 64

void Sort_parallel(vector<double> &data)
{
    const intptr_t size = data.size();
    intptr_t stride = DATASET_SIZE / THREAD_NUM;
    //对data进行分线程排序
    if (stride != 1)
    {
#pragma omp parallel for shared(data, stride) schedule(dynamic)

        for (intptr_t i = 0; i < size; i += stride)
        {
            auto left = data.begin() + i;
            auto right = i + stride < size ? data.begin() + i + stride : data.end();
            std::sort(left, right);
        }
    }
    //并行归并排序
#pragma omp parallel
    {
        intptr_t _stride = stride;
        do
        {
            _stride *= 2; // 两个数组步长之和
#pragma omp for schedule(dynamic)
            for (intptr_t i = 0; i < size; i += _stride)
            {
                // 合并[i, i+_stride/2], [i+_stride/2, i+_stride]
                auto left = data.begin() + i;
                auto mid = (i + i + _stride) / 2 < size ? data.begin() + (i + i + _stride) / 2 : data.end();
                auto right = i + _stride < size ? data.begin() + i + _stride : data.end();
                inplace_merge(left, mid, right);
            }
        } while (_stride < size);
    }
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

int main()
{

    const char *logfilename = "log.txt";
    vector<double> second;
    vector<double> a, backup;
    backup.resize(DATASET_SIZE);
    a.resize(DATASET_SIZE);
    // acopy.resize(DATASET_SIZE);
    for (int i = 0; i < DATASET_SIZE; i++)
    {
        a[i] = rand() / double(RAND_MAX);
        // acopy[i] = a[i];
        backup[i] = a[i];
    }

    // 创建二维数组
    long size_each_vec = DATASET_SIZE / THREAD_NUM;
    vector<vector<double>> matrix(THREAD_NUM, vector<double>(size_each_vec));

    // 串行计算
    clock_t start, end;
    // 开始计时
    start = clock();
    cout << "start sorting " << endl
         << endl;
    sort(a.begin(), a.end(), less<double>());
    // for (int i = 0; i < DATA_SIZE; i++)
    //	cout << a[i] << endl;
    end = clock();
    double duration = double(end - start) * 1000 / CLOCKS_PER_SEC;
    double chuan = duration;
    cout << " 串行用时 = " << duration << " ms" << endl;
    // cout << a[0] << "  " << a[DATASET_SIZE / 2] << "  " << a[DATASET_SIZE - 1] << endl;
    ofstream log;
    log.open(logfilename, ios::out | ios::app);
    log << "串行运算Sort()函数用时： " << duration << " ms" << endl
        << endl;

    //并行计算
    // vector<double> re;
    // re.resize(DATASET_SIZE);
    omp_set_num_threads(THREAD_NUM);
    for (int times = 1; times < 6; times++)
    {
        cout << times << endl;
        vector<double> acopy(backup);
        start = clock();
        Sort_parallel(acopy);
        end = clock();

        duration = double(end - start) * 1000 / CLOCKS_PER_SEC;
        second.push_back(duration);
        log << THREAD_NUM << " 线程第 " << times << " 次并行运算用时： " << duration << " ms" << endl;
        cout << THREAD_NUM << " 线程第 " << times << " 次并行运算用时： " << duration << " ms" << endl
             << endl;
        bool same = compare(a, acopy);
        if (same)
            cout << "并行计算结果正确" << endl;
    }

    int total = 0;
    for (int k = 0; k < 5; k++)
    {
        total += second[k];
    }
    total = total / 5;
    double rate = chuan / total;
    cout << "平均值 " << total << " ms" << endl;
    log << "平均值 " << total << " ms" << endl;
    cout << "加速比 " << rate << endl;
    log << "加速比 " << rate << endl;
}
