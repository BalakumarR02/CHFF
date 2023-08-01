#include <bits/stdc++.h>

using namespace std;
using ll = long long;

void solve()
{
    ll n, i, j;

    cin >> n;
    vector<double> phi(n + 1, 0), b(n + 1, 0);
    vector<vector<double>> a(n + 1, vector<double>(n + 1, 0));
    // To get the matrix
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            cin >> a[i][j];
        }
    }
    // Getting the 'b' values
    for (i = 1; i <= n; i++)
    {
        cin >> b[i];
    }
    // Getting inital values
    double ini;
    cin >> ini;
    // initalizing phi
    for (i = 1; i <= n; i++)
    {
        phi[i] = ini;
    }
    double err = 1e9;
    ll iter = 0;
    // iterating till err is minimum or iteration greater tha 1e5-> to indiacte it diverges
    while (err > 1e-5 && iter <= 1e5)
    {
        err = 0.0;
        // Doing the gauss sidel iteration for all phi values
        for (i = 1; i <= n; i++)
        {
            double x = 0.0, val = 0.0;
            for (j = max(i - 1, 1LL); j <= min(i + 1, n); j++)
            {
                if (i == j)
                {
                    x = a[i][j];
                }
                else
                {
                    val += (-1 * a[i][j] * phi[j]);
                }
            }
            val += b[i];
            val /= x;
            err += (val - phi[i]) * (val - phi[i]);
            phi[i] = val;
            cout << iter << " " << fixed << setprecision(4) << val << " ";
        }
        err = sqrt(err);
        iter++;
        cout << err << "\n";
    }
    if (iter > 1e5)
    {
        cout << "Diverges!!!\n";
        return;
    }
    // printing values
    for (i = 1; i <= n; i++)
    {
        cout << phi[i] << " ";
    }

    cout << "\n";
}
int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    solve();

    return 0;
}

// sort(begin(v), end(v), [] (int a, int b) { return a > b; });           (Custom sort using lambda function)

// 5
// 3 -1 0 0 0
// -1 2 -1 0 0
// 0 -1 2 -1 0
// 0 0 -1 2 -1
// 0 0 0 -1 3
// 204 12 20 28 1236
// 200
