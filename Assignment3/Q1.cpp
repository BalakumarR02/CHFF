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

double x_vel(double x, double rho, double del)
{
    return rho * (x * x + 1) * del;
}
double y_vel(double x, double rho, double del)
{
    return rho * (x * x + 1) * del;
}
double m_sour(double x, double y)
{
    return 4 * (x + y);
}

void solve()
{
    ll n, i;

    cin >> n;

    cout << "\n";
}
int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    ll i, j, n, m, k;
    double L = 3.0, rho = 2.0, phi_s = 0.5, phi_b = 1.0;
    cin >> n;
    m = n;
    vector<vector<vector<double>>> a(n + 1, vector<vector<double>>(n + 1, vector<double>(7, 0.0)));
    vector<vector<double>> phi(n + 2, vector<double>(n + 2, 0.0));
    double del_x = L / n, del_y = L / m;

    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            double F_e = x_vel(j * del_x, rho, del_y), F_w = x_vel((j - 1) * del_x, rho, del_y), F_n = y_vel(i * del_y, rho, del_x), F_s = y_vel((i - 1) * del_y, rho, del_x);
            double ap = F_e + F_n, ae = 0, an = 0, aw = F_w, as = F_s, b = phi_s * m_sour((i - 0.5) * del_x, (j - 0.5) * del_x) * del_x * del_y;
            if (i == 1)
            {
                as = 0;
                b += (F_s * phi_b);
            }
            if (j == 1)
            {
                aw = 0;
                b += (F_w * 0);
            }
            a[i][j][1] = ap, a[i][j][2] = ae, a[i][j][3] = aw, a[i][j][4] = an, a[i][j][5] = as, a[i][j][6] = b;
        }
    }
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            phi[i][j] = a[i][j][2] * phi[i][j + 1] + a[i][j][3] * phi[i][j - 1] + a[i][j][4] * phi[i + 1][j] + a[i][j][5] * phi[i - 1][j] + a[i][j][6];

            phi[i][j] = (phi[i][j] / a[i][j][1]);
        }
    }
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {

            cout << phi[n - i + 1][j] << " ";
        }
        cout << "\n";
    }

    return 0;
}

// sort(begin(v), end(v), [] (int a, int b) { return a > b; });           (Custom sort using lambda function)

