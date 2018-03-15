#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <list>
#include <queue>
#include <utility>
#include <unordered_map>
#include <iostream>
using namespace std;

class Template
{
public:
	vector<bool> Eratosthenes_Sieve(int n);										//Eratosthenes筛法
	int gcd(int a, int b);														//欧几里得算法
	void ex_gcd(int a, int b, int &x, int &y);									//扩展欧几里得算法
	vector<int> Prime_Factor(int n);											//分解质因数
	int Longest_substring(string s);											//最长无重复子串
	void merge_sort(vector<int> &target);										//归并排序
	void quick_sort(vector<int> &target);										//快速排序
	vector<int> topological_sort(vector<list<int>> adjacency_list);				//拓扑排序
private:
	void merge_sort_recursive(vector<int> &target, std::vector<int> &copy, size_t start, size_t end);
	void quick_sort_recursive(vector<int> &target, int start, int end);
};

/*Eratosthenes筛法
**筛选出n以内的所有质数
**返回参数res中如果res[i] == false则i为质数
***解释:对于p<=n && p>1的所有整数p,标记所有1p,2p,3p,4p......则未标记的数即为质数
***    !res[i]:只需判断p为素数的情况，若p非素数则p与p的倍数在之前的循环已经标记过
***    i * i:内层循环只需从i*i开始因为之前的循环已经标记过i * x(x<i)的情况
*/

vector<bool> Template::Eratosthenes_Sieve(int n)
{
	vector<bool> res(n+1, false);
	res[0] = true;
	res[1] = true;
	for (int i = 2; i <= n; i++)
	{
		if (!res[i])
		{
			for (int j = i * i; j <= n; j += i)
			{
				res[j] = true;
			}
		}
	}
	return res;
}

/*欧几里得算法（辗转相除法）
**找出a与b的最大公约数
**返回参数为a与b的最大公约数
*/

int Template::gcd(int a, int b)
{
	return b == 0 ? a : gcd(b, a%b);
}

/*扩展欧几里得算法
**找出ax + by = gcd(a,b)的一个x,y整数解
**参数x,y即为上述整数的一对整数解
***解释：
***		推理1：当b = 0时ax + by = gcd(a,b) = a,此时x = 1,取y = 0;
***		推理2：设ax1 + by1 = gcd(a,b), bx2 + a%by2 = gcd(b,a%b) 由欧几里得算法递归可知gcd(a,b) = gcd(b,a%b)
***			  则可得等式a(x1) + b(y1) = a(y2) + b(x2 - (a/b)*y2)视a,b为未知数由等式恒等定理可得
***		   	  递推关系 x1 = y2 , y1 = x2 - (a/b)*y2;
*/

void Template::ex_gcd(int a, int b, int &x, int &y)
{
	if (b == 0)
	{
		x = 1;
		y = 0;
		return;
	}
	int x1, y1;
	ex_gcd(b, a%b, x1, y1);
	x = y1;
	y = x1 - (a / b)*y1;
}

/*分解质因数（唯一分解定理）
**将整数n用多个质数相乘的形式表示
**返回参数res中的元素即为n的质数因子
***算术基本定理可表述为：任何一个大于1的自然数 N,如果N不为质数，那么N可以唯一分解成有限个质数的乘积;
*/

vector<int> Template::Prime_Factor(int n)
{
	vector<int> res;
	for (int i = 2; i <= n; i++)
	{
		while (n%i == 0)
		{
			res.push_back(i);
			n /= i;
		}
	}
	return res;
}

/*最长无重复子串
**找出串s的最长无重复子串 例如"abcabcbb"的最长无重复子串为"abc"长度为3
**返回参数即为最长无重复子串的长度
***
*/

int Template::Longest_substring(string s)
{
	int length = s.length(), res = 0;
	unordered_map<char, int> hash_map; 
	int low = 0;
	for (int high = 0 ; high < length; high++)
	{
		auto itr = hash_map.find(s[high]);
		if (itr != hash_map.end())
		{
			low = max(itr->second, low);
		}
		bool flag;
		flag = hash_map.insert({ s[high], high + 1 }).second;
		if (!flag)
		{
			hash_map[s[high]] = high + 1;
		}
		res = max(res, high + 1 - low);
	}
	return res;
}

/*归并排序
**以归并排序的方法排序容器target
***
*/

void Template::merge_sort(vector<int> &target)
{
	vector<int> copy = target;
	merge_sort_recursive(target, copy, 0, target.size() - 1);
}
void Template::merge_sort_recursive(vector<int> &target, std::vector<int> &copy, size_t start, size_t end)
{
	if (start >= end) return;
	int mid = (end - start + 1) / 2 + start;
	merge_sort_recursive(target, copy, start, mid - 1);
	merge_sort_recursive(target, copy, mid, end);
	int start1 = start, start2 = mid, counter = start;
	while (start1 <= mid - 1 && start2 <= end)
		target[counter++] = copy[start1] < copy[start2] ? copy[start1++] : copy[start2++];
	while (start2 <= end)
		target[counter++] = copy[start2++];
	while (start1 <= mid - 1)
		target[counter++] = copy[start1++];
	for (int i = start; i <= end; i++)
	{
		copy[i] = target[i];
	}
}

/*快速排序
**以快速排序的方法排序容器vector
*/

void Template::quick_sort(vector<int> &target)
{
	quick_sort_recursive(target, 0, target.size() - 1);
}
void Template::quick_sort_recursive(vector<int> &target, int start, int end)
{
	if (start >= end)
		return;
	int pivot_element = target[end];
	int flag = start;
	for (int j = start; j <= end - 1; j++)
	{
		if (target[j] < pivot_element)
			std::swap(target[flag++], target[j]);
	}
	std::swap(target[flag], target[end]);
	quick_sort_recursive(target, start, flag - 1);
	quick_sort_recursive(target, flag + 1, end);
}

/*拓扑排序
**以带队列的方式对图进行拓扑排序
**返回参数为排序后顶点的顺序
***解释：0：将所有入度为0的顶点存入队列
***		1：不断的弹出队列中的顶点元素，每弹出一个顶点元素，标记此顶点并将计数器加1，然后通过邻接列表访问此顶点指向的所有顶点
***		3：将每个顶点的入度减1，若减1后入度为0则将此顶点存入队列 返回第0步
*/

vector<int> Template::topological_sort(vector<list<int>> adjacency_list)
{
	map<int, int> vertices_indgree;										//全部顶点的入度表 first为顶点名称 second为此顶点的入度
	vector<int> res;													//拓扑排序后所有顶点的下标表 first为顶点名称 second为此顶点所处的位置
	bool cycle_found = false;											//若检测到图中有环则cycle_found = true；

	//构建入度表
	for (int i = 0; i < adjacency_list.size(); i++)
	{
		vertices_indgree.insert({ i,0 });
	}
	for (int i = 0; i < adjacency_list.size(); i++)
	{
		for (auto itr = adjacency_list[i].begin(); itr != adjacency_list[i].end(); itr++)
		{
			vertices_indgree[*itr]++;
		}
	}
	//构建完毕

	queue<pair<int, int>> que;
	int counter = 0;
	for (int i = 0; i < vertices_indgree.size(); i++)
	{
		if (vertices_indgree[i] == 0)
		{
			que.push({ i,vertices_indgree[i] });
		}
	}
	while (!que.empty())
	{
		pair<int, int> vertice = que.front();
		que.pop();
		res.push_back(vertice.first);
		counter++;
		vertices_indgree[vertice.first] = -1;						//标记此顶点入度为-1以确保不会在访问此顶点
		for (auto itr = adjacency_list[vertice.first].begin(); itr != adjacency_list[vertice.first].end(); itr++)
		{
			if (--vertices_indgree[*itr] == 0)
			{
				que.push({ *itr,vertices_indgree[*itr] });
			}
		}
	}
	if (counter != vertices_indgree.size())
	{
		cycle_found = true;
	}
	
	return res;
}