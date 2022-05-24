#include <omp.h>
//#include<stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;
int t_num = 2;

int main()
{

    int range;           //范围
    int total;           //因子和
    bool correct = true; // 判断并行计算结果是否相等
    // cin >> range;
    // int* result1 = new int[range + 1];//串行计算结果
    // int* result2 = new int[range + 1];//并行计算结果
    range = 5000000;
    vector<int> result1;
    vector<int> result2;
    clock_t start, end;
    /*串行运算*/
    clock_t start = clock();
    for (int i = 2; i <= range; i++)
    {
        total = 1;
        for (int j = 2; j <= sqrt(i); j++)
        {
            if (i % j == 0) //
            {
                total += j;
                if (j != sqrt(i))
                    total += i / j;
            }
        }
        if (total == i)
            result1.push_back(i);
    }
    clock_t end = clock();

    cout << "正确答案：";
    for (int i = 0; i < result1.size(); i++)
    {
        cout << result1[i] << " ";
    }
    double duration1 = end - start;
    cout << endl
         << "串行计算耗时：" << duration1 << "ms" << endl;
    double duration1 = 27483;
    double totaltime = 0;
    for (int times = 1; times <= 5; times++)
    {
        /*并行运算*/
        start = clock();
        correct = true;
        result2.clear();

        omp_set_num_threads(t_num);
#pragma omp parallel
        {
#pragma omp for schedule(dynamic)
            for (int i = 2; i <= range; i++)
            {
                total = 1;
                for (int j = 2; j <= sqrt(i); j++)
                {
                    if (i % j == 0) //
                    {
                        total += j;
                        if (j != sqrt(i))
                            total += i / j;
                    }
                }
                if (total == i)
                    result2.push_back(i);
            }
        }

        end = clock();

        cout << "并行答案：";
        for (int i = 0; i < result2.size(); i++)
        {
            cout << result2[i] << " ";
        }
        if (result1.size() != result2.size())
        {
            correct = false;
            cout << result2.size() << endl;
            for (int k = 0; k < result2.size(); k++)
            {
                cout << result2[k] << " ";
            }
            cout << "num not correct" << endl;
        }
        else
        {
            sort(result2.begin(), result2.end()); // 避免乱序

            for (int i = 0; i <= result1.size(); i++) // 检验并行计算结果正确性
            {
                if (result1[i] != result2[i])
                {
                    correct = false;
                    break;
                }
            }
        }

        if (correct)
        {
            double duration2 = end - start;
            totaltime += duration2;
            cout << endl
                 << "thread并行计算耗时：" << duration2 << "ms, 加速比:" << duration1 / duration2 << endl;
        }
    }
    totaltime = totaltime / 5;
    cout << "平均时间 " << totaltime << endl;
    double rate = duration1 / totaltime;
    cout << "加速比 = " << rate << endl;
}
