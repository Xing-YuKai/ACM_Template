// ACM_Template.cpp: 定义控制台应用程序的入口点。
// 开发环境 Visual Studio 2017
// 若非此开发环境删除line 5 #include "stdafx.h"

#include "stdafx.h"
#include "Template.h"

int main()
{
	Template obj;
	/*unit test for Eratosthenes Sieve
	int n;
	cin >> n;
	vector<bool> prime_number_list = obj.Eratosthenes_Sieve(n);
	for (int i = 0; i <= n; i++)
	{
		if (!prime_number_list[i])
			cout << i << " ";
	}
	*/

	/*unit test for gcd
	int a, b;
	cin >> a >> b;
	cout << obj.gcd(a, b);
	*/

	/*unit test for ex_gcd
	int a, b, x, y;
	cin >> a >> b;
	obj.ex_gcd(a, b, x, y);
	cout << x << " " << y;
	*/
	
	/*unit test for Prime_Factor
	int n;
	cin >> n;
	vector<int> res = obj.Prime_Factor(n);
	for (int a : res)
	{
		cout << a << " ";
	}
	*/

	/*unit test for longest_substring
	string s;
	cin>>s;
	cout << obj.Longest_substring(s);
	*/

	/*unit test for topological_sort
	vector<list<int>> adjacency_list;
	adjacency_list.resize(7);
	adjacency_list[0] = { 1,3,2 };
	adjacency_list[1] = { 3,4 };
	adjacency_list[2] = { 5 };
	adjacency_list[3] = { 5,6,2 };
	adjacency_list[4] = { 3,6 };
	adjacency_list[5] = {};
	adjacency_list[6] = { 5 };
	vector<int> res = obj.topological_sort(adjacency_list);
	for (int i = 0; i < res.size(); i++)
	{
		cout << res[i] << " ";
	}
	*/

	/*unit test for union_find
	Template::union_find obj_UF(8);
	obj_UF.Union(2, 3);
	obj_UF.Union(1, 0);
	obj_UF.Union(0, 4);
	obj_UF.Union(5, 7);
	obj_UF.Union(6, 2);
	obj_UF.Union(1, 4);
	cout << "连通分量数目为：" << obj_UF.get_counter() << endl;
	cout << (obj_UF.connected(6, 3) ? "点6、3相连" : "Exception");
	*/

	/*unit test for manacher
	string s, sym;
	getline(cin, s);
	pair<int, int> res = obj.manacher(s);
	cout << res.first << endl;
	for (int i = res.first; i > 0; i--)
	{
		if (s[res.second - i] != (char)1)
			cout << s[res.second - i];
	}
	if (s[res.second] != (char)1)
		cout << s[res.second];
	for (int i = 1; i <= res.first; i++)
	{
		if (s[res.second + i] != (char)1)
			cout << s[res.second + i];
	}
	*/

	/*unit test for KMP
	string a, b;
	getline(cin, a);
	getline(cin, b);
	cout << obj.KMP(a, b);
	*/

	/*unit text for Dijkstra
	**将此段写在main函数外
	vector<int> Distance(7, -1);
	vector<bool> known(7, false);
	vector<vector<int>> adjacency_matrix = { { 0,2,0,1,0,0,0 },{ 0,0,0,0,3,10,0 },{ 4,0,0,0,5,0 },{ 0,0,2,0,2,8,4 },{ 0,0,0,0,0,0,6 },{ 0,0,0,0,0,0,0 },{ 0,0,0,0,0,1,0 } };
	struct cmp
	{
		bool operator()(int a, int b)
		{
			return Distance[a] > Distance[b];
		}
	};
	**将此段写在main函数外
	int source = 0;
	Distance[source] = 0;
	priority_queue<int, vector<int>, cmp> que;
	que.push(source);
	while (!que.empty())
	{
		int tmp = que.top();
		que.pop();
		known[tmp] = true;
		for (int i = 0; i < adjacency_matrix[tmp].size(); i++)
		{
			if (adjacency_matrix[tmp][i] == 0 || known[i])
				continue;
			if (Distance[i] == -1 || Distance[tmp] + adjacency_matrix[tmp][i] < Distance[i])
			{
				Distance[i] = Distance[tmp] + adjacency_matrix[tmp][i];
				que.push(i);
			}
		}
	}
	for (auto i : Distance)
		cout << i << endl;
	*/

	/*unit test for Bash_game
	if (obj.Bash_game(20, 3))
	{
		cout << "先手胜" << endl;
	}
	else
	{
		cout << "后手胜" << endl;
	}
	if (obj.Bash_game(20, 2))
	{
		cout << "先手胜";
	}
	else
	{
		cout << "后手胜";
	}
	*/

    return 0;
}

