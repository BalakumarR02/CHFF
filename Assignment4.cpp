#include <bits/stdc++.h>

using namespace std;
using ll = long long;
#define lb lower_bound
#define ub upper_bound
#define all(v) (v).begin(), (v).end()
#define all1(v) (v).begin() + 1, (v).end()
#define allr(v) (v).rbegin(), (v).rend()
#define allr1(v) (v).rbegin() + 1, (v).rend()
#define sort0(v) sort(all(v))
#define fo(i, a, b) for (i = a; i <= b; i++)
#define fi(i, a, b) for (i = a; i >= b; i--)

// freopen("input.txt", "r", stdin);
// freopen("output.txt", "w", stdout);
#define rep(i, begin, end) for (__typeof(end) i = (begin) - ((begin) > (end)); i != (end) - ((begin) > (end)); i += 1 - 2 * ((begin) > (end)))
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<ll> vll;
typedef pair<ll, ll> pll;

#define sz(x) (ll) x.size()
#define pb push_back
#define ppb pop_back
#define mkp make_pair
#define inf 1000000000000000005

const ll mod = 1e9 + 7;
const ll mod1 = 998244353;
const ll N = 1e5; // limit for array size

template <class T>
std::vector<std::vector<T>> Multiply(std::vector<std::vector<T>> &a, std::vector<std::vector<T>> &b)
{
    const int n = a.size();    // a rows
    const int m = a[0].size(); // a cols
    const int p = b[0].size(); // b cols

    std::vector<std::vector<T>> c(n, std::vector<T>(p, 0));
    for (auto j = 0; j < p; ++j)
    {
        for (auto k = 0; k < m; ++k)
        {
            for (auto i = 0; i < n; ++i)
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}

void solve()
{
    ll n, i, j;

    cin >> n;
    vector<vector<double>> dx(n + 1, vector<double>(n + 1)), dy(n + 1, vector<double>(n + 1));

    dx[1][1] = 0, dx[1][2] = 1, dx[1][3] = 0.5, dx[1][4] = 0.5;
    dx[2][1] = dx[1][2], dx[2][2] = 0, dx[2][3] = 0.5, dx[2][4] = 0.5;
    dx[3][1] = dx[1][3], dx[3][2] = dx[2][3], dx[3][3] = 0, dx[3][4] = 0;
    fo(i, 1, 4)
    {
        dx[4][i] = dx[i][4];
    }

    dy[1][1] = 0, dy[1][2] = 0, dy[1][3] = (double)sqrt(3) / 2, dy[1][4] = (double)1 / (2 * sqrt(3));
    dy[2][1] = dy[1][2], dy[2][2] = 0, dy[2][3] = (double)sqrt(3) / 2, dy[2][4] = (double)1 / (2 * sqrt(3));
    dy[3][1] = dy[1][3], dy[3][2] = dy[2][3], dy[3][3] = 0, dy[3][4] = (double)1 / (sqrt(3));
    fo(i, 1, 4)
    {
        dy[4][i] = dy[i][4];
    }
    vector<vector<double>> M4(n, vector<double>(2)), M4T(2, vector<double>(n));
    fo(i, 1, 3)
    {
        M4[i][1] = dx[4][i];
        M4[i][2] = dy[4][i];

        M4T[1][i] = dx[4][i];
        M4T[2][i] = dy[4][i];
    }

    cout << "\n";
}
int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    ll t;

    cin >> t;
    while (t--)
    {
        solve();
    }
    return 0;
}

// sort(begin(v), end(v), [] (int a, int b) { return a > b; });           (Custom sort using lambda function)
