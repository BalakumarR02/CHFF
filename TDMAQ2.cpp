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
    ll i, j;
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

void TDMA2D(ll n, ll m, vector<vector<vector<double>>> &a, ll sweep, ofstream &MyFile)
{
    ll i, j, k;
    double err = 1e9;

    ll iter = 0;
    vector<vector<double>> phi(n + 2, vector<double>(m + 2, 200));
    while (err >= 1e-6 && iter <= 1e5)
    {
        err = 0.0;
        fo(i, 1, m)
        {
            vector<vector<double>> tda(n + 1, vector<double>(4, 0));
            ll ind = i;
            for (j = 1; j <= n; j++)
            {
                double val = 0;
                ll f1 = j, f2 = ind;
                for (k = 1; k < 6; k++)
                {
                    if (k == sweep)
                    {
                        val += a[f1][f2][k] * phi[j][ind + 1];
                        continue;
                    }
                    else if (k == sweep + 1)
                    {
                        val += a[f1][f2][k] * phi[j][ind - 1];
                        continue;
                    }
                    tda[j][max(1LL, k - 2)] = a[f1][f2][k];
                }
                tda[j][4] = val + a[f1][f2][6];
            }
            vector<double> ph1 = TDMA1D(n, tda);

            fo(j, 1, n)
            {
                err += (phi[j][ind] - ph1[j - 1]) * (phi[j][ind] - ph1[j - 1]);
                phi[j][ind] = ph1[j - 1];
            }
        }
        err = sqrt(err);
        iter++;
    }
    if (iter > 1e5)
    {
        cout << "Diverges!!!\n";
        return;
    }
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= n; j++)
        {
            MyFile << setprecision(4) << fixed << phi[j][i] << " ";
        }
        MyFile << "\n";
    }
}
vector<vector<double>> TDMAQ2_x(ll n, ll m, vector<vector<vector<double>>> &a, ofstream &MyFile)
{
    ll i, j, k;
    double err = 1e9;

    ll iter = 0;
    vector<vector<double>> phi(n + 2, vector<double>(m + 2, 0));
    while (err >= 1e-9 && iter <= 1e5)
    {
        err = 0.0;
        for (i = 1; i <= m; i++)
        {
            vector<vector<double>> matr(n + 1, vector<double>(5, 0));
            for (j = 1; j <= n; j++)
            {
                matr[j][1] = a[j][i][1];
                matr[j][2] = a[j][i][4];
                matr[j][3] = a[j][i][5];
                matr[j][4] = a[j][i][6] + a[j][i][2] * phi[j][i + 1] + a[j][i][3] * phi[j][i - 1];
            }
            vector<double> ph1 = TDMA1D(n, matr);
            fo(j, 1, n)
            {
                err += (phi[j][i] - ph1[j]) * (phi[j][i] - ph1[j]);
                phi[j][i] = ph1[j];
            }
        }
        err = sqrt(err);
        iter++;
    }
    return phi;
}
int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    ofstream MyFile("cv51x51.txt");

    ll n, m, i, j, p = 6, sweep = 2;
    cin >> n >> m;
    vector<vector<vector<double>>> a(n + 1, vector<vector<double>>(m + 1, vector<double>(7, 0)));
    // To get the matrix
    double delx = ((double)1 / m), dely = ((double)1 / n), k = 1.0, h = 10.0, ueq = 0;
    ueq = (h * ((double)2 * k / (dely))) / (h + ((double)2 * k / (dely)));

    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= n; j++)
        {

            double cons = 0.0, an = k * (delx / dely), ae = k * (dely / delx), aw = k * (dely / delx), as = k * (delx / dely), ab = 0, ap = 0;

            if (j == 1)
            {
                an = 0;
                ab += (2 * k * delx / (dely));
            }
            else if (j == n)
            {
                as = 0;
                ap += (ueq * delx);
                cons += (ueq * delx) * 300;
            }
            if (i == 1)
            {
                aw = 0;
                ab += (2 * k * dely / (delx));
            }
            else if (i == m)
            {
                ae = 0;
                ab += (2 * k * dely / (delx));
            }
            double src = 0.0;
            src = ((100.0 + 500.0 * ((i * delx) + ((n - j + 1) * (dely)) - ((delx)))) * delx * dely);
            cons += src + (ab * 500.0);
            ap += an + ae + aw + as + ab;
            a[j][i][1] = ap, a[j][i][2] = ae, a[j][i][3] = aw, a[j][i][4] = as, a[j][i][5] = an, a[j][i][6] = cons;
        }
    }

    vector<vector<double>> phi = TDMAQ2_x(n, m, a, MyFile);
    for (i = 1; i <= n; i++)
    {
        MyFile << '[';
        for (j = 1; j <= m; j++)
        {
            MyFile << setprecision(4) << fixed << phi[i][j] << ", ";
        }
        MyFile << "],\n";
    }
    MyFile.close();

    return 0;
}

// sort(begin(v), end(v), [] (int a, int b) { return a > b; });           (Custom sort using lambda function)
