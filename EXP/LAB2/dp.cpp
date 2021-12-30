/*
 * 多重背包问题：
 * INPUT:
 * 第一行n,W分别标是物品件数和背包容量
 * 之后n行表示每种物品的weight value number
 * OUTPUT:
 * 能装入背包的最大价值
 */

#include <iostream>
using namespace std;
#define MAX 200
int weight[MAX + 1];
int value[MAX + 1];
int num[MAX + 1];
int dp[200+1][10000+1];
int max(int a, int b) {
    return a > b ? a : b;
}
int main() {
    int n, m;
    cin >> n;
    cin >> m;
    for (int i = 1; i <= n; i++) {
        cin >> weight[i] >> value[i] >> num[i];
    }
    for (int i = 0; i < n; ++i) {
        dp[i][0] = 0;
        dp[0][i] = 0;
    }
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int max = -1;
            for (int k = 0; k <= num[i]&&k*weight[i]<=j; ++k) {
                int tmp = dp[i-1][j-k*weight[i]]+k*value[i];
                if (tmp>max){
                    max = tmp;
                }
            }
            dp[i][j] = max;
        }
    }
    cout << dp[n][m];
}
