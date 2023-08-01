#include <bits/stdc++.h>
#include <fstream>

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
vector<double> TDMA1D(ll n, vector<vector<double>> &a)
{
    ll i, j, m;
    m = 4;
    vector<double> phi(n + 2, 0);

    vector<double> P(n + 1, 0), Q(n + 1, 0);
    fo(i, 1, n)
    {
        P[i] = a[i][2] / (a[i][1] - a[i][3] * P[i - 1]);
        Q[i] = (a[i][3] * Q[i - 1] + a[i][4]) / (a[i][1] - a[i][3] * P[i - 1]);
    }
    fi(i, n, 1)
    {
        phi[i] = P[i] * phi[i + 1] + Q[i];
    }
    return phi;
}

int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    ll n, m, i, j;
    cin >> n;
    vector<vector<double>> a(n + 1, vector<double>(5, 0));

    // To get the matrix

    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= 4; j++)
        {

            cin >> a[i][j];
        }
    }
    vector<double> phi = TDMA1D(n, a);
    fo(i, 1, n)
    {
        cout << phi[i] << " ";
    }

    return 0;
}

// sort(begin(v), end(v), [] (int a, int b) { return a > b; });           (Custom sort using lambda function)

// 5
// 20  5 0  1100
// 15  5 5  100
// 15  5 5  100
// 15  5 5  100
// 10  0 5  100
// 64.2276 36.9106 26.5041 22.6016 21.3008w