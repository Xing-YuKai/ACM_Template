#pragma once
#include <vector>
using namespace std;

class Template
{
public:
	vector<bool> Eratosthenes_Sieve(int n);					//Eratosthenes筛法
	int gcd(int a, int b);									//欧几里得算法
	void ex_gcd(int a, int b, int &x, int &y);				//扩展欧几里得算法
	vector<int> Prime_Factor(int n);						//分解质因数
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
		return; p;.gy
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