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

class Fin_variables
{
public:
    double A, Tf, Ab, Af, Ap, As, Auf, Atot, P, L1, W1, B, H, t, L, W, Tb, Tinf, N, S, qb;
    Fin_variables()
    {
        Tinf = 293;
        Tb = 373;
        H = 1000 * 0.001;
        L = 200 * 0.001;
        L1 = 1;
        W1 = 1;
        W = 1;
        t = 10 * 0.001;
        N = 6;
        S = 4;
        qb = 10;

        solve();
    }

    void change_fin_dim(double H, double L, double W, double t, double N, double S)
    {
        this->L = L;
        this->H = H;
        this->W = W;
        this->t = t;
        this->N = N;
        this->S = S;
        solve();
    }
    void solve()
    {
        P = 2 * (t + L);
        A = W1 * L1;
        Af = L * t;
        As = P * H;
        Ap = t * H;
        Auf = L * W - (N * Af);
        Atot = Auf + As;
        Tf = (Tinf + Tb) / 2;
    }
};

class Flow_Parameters : public Fin_variables
{
public:
    double k, ka, NuL, Nux, Rex, Pe, Pr, q, RaL, ReL, Tb, Tinf, mu, nu, h, qs, uinf, rho, hx;
    Flow_Parameters()
    {
        Tinf = 293;
        Tb = 373;
        qs = 0;
        k = 247;
        h = 80;
        uinf = 1.6;
        hx = 80;
        solve();
    }

    void solve()
    {
        Pr = 0.7202;
        rho = 1.06;
        nu = 1.896 * 1e-5;
        mu = 2.008 * 1e-5;
        ka = 0.02808;

        // Nux = 0.453
    }
    void solve_h(double x)
    {
        Rex = uinf * x / nu;
        Nux = 0.332 * pow(Rex, (double)1 / 2) * pow(Pr, (double)1 / 3);
        hx = Nux * ka / x;
    }
};
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

vector<vector<double>> TDMA2D(ll n, ll m, vector<vector<vector<double>>> &a, ll sweep)
{
    ll i, j, k;
    double err = 1e9;
    ll sw = sweep;
    if (sweep == 3 || sweep == 5)
    {
        sweep--;
    }
    if (sweep > 3)
    {
        // sweep--;
        swap(n, m);
    }
    ll iter = 0;
    vector<vector<double>> phi(n + 2, vector<double>(m + 2, 0));
    while (err >= 1e-6 && iter <= 1e5)
    {
        err = 0.0;
        fo(i, 1, m)
        {
            ll ind = (sw == 2 || sw == 4) ? i : m - i + 1;
            vector<vector<double>> tda(n + 1, vector<double>(5, 0));

            for (j = 1; j <= n; j++)
            {
                double val = 0;
                ll f1 = j, f2 = ind;
                if (sw > 3)
                {
                    f1 = ind;
                    f2 = j;
                }
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
                    tda[j][max(1LL, k - ((sw < 4) ? 2 : 0))] = a[f1][f2][k];
                }
                tda[j][4] = val + a[f1][f2][6];
            }

            vector<double> ph1 = TDMA1D(n, tda);

            fo(j, 1, n)
            {
                err += (phi[j][ind] - ph1[j]) * (phi[j][ind] - ph1[j]);
                phi[j][ind] = ph1[j];
            }
        }
        err = sqrt(err);

        iter++;
    }
    if (iter > 1e5)
    {
        cout << "Diverges!!!\n";
    }
    // for (i = 1; i <= m; i++)
    // {
    //     for (j = 1; j <= n; j++)
    //     {
    //         cout << setprecision(4) << fixed << phi[j][i] << " ";
    //     }
    //     cout << "\n";
    // }
    return phi;
}

void fin_1D_sample()
{
    ofstream MyFile("plate_fin1.txt");

    ll n, i, m, j;
    Fin_variables f1;
    Flow_Parameters c1;
    // Insulated B.C
    cin >> n >> m;

    vector<vector<vector<double>>> v(n + 1, vector<vector<double>>(m + 1, vector<double>(7)));
    double del_x = f1.L / n, del_y = f1.H / m;

    fo(i, 1, n)
    {

        double x_conv = del_x * (i - 0.5);
        c1.solve_h(x_conv);

        cout << c1.hx << "\n";
        fo(j, 1, m)
        {
            double ap = 0, aw = del_y * c1.k * f1.t / del_x, ae = del_y * c1.k * f1.t / del_x, an = del_x * c1.k * f1.t / del_y, as = del_x * c1.k * f1.t / del_y, ab = 0, b = 0;

            double P = 2 * del_y, Ueq = c1.h * (2 * c1.k / del_y) / (c1.h + (2 * c1.k / del_y)), msq = c1.h * 2 * (del_x + del_y) / (del_x * del_y);
            if (i == 1)
            {
                ab += 2 * as;
                as = 0;
                b += ab * f1.Tb;
            }
            else if (i == n)
            {
                an = 0;
            }

            if (j == 1)
            {
                ab += Ueq * f1.t * del_y;
                aw = 0;
                b += Ueq * c1.Tinf * f1.t * del_y;
            }
            else if (j == m)
            {
                ab += Ueq * f1.t * del_y;
                ae = 0;
                b += Ueq * c1.Tinf * f1.t * del_y;
            }

            ap = aw + ae + an + as + ab + (msq * del_x * del_y * f1.t);
            b += msq * del_x * del_y * f1.t * f1.Tinf;
            v[i][j][1] = ap, v[i][j][2] = ae, v[i][j][3] = aw, v[i][j][4] = an, v[i][j][5] = as, v[i][j][6] = b;
        }
    }
    vector<vector<double>> phi = TDMA2D(n, m, v, 4);

    // Temperature and heat flux at different points
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= n; j++)
        {
            MyFile << setprecision(4) << fixed << phi[j][i] << " ";
        }
        MyFile << "\n";
    }
    MyFile << "\n";
}
int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    fin_1D_sample();

    return 0;
}

// sort(begin(v), end(v), [] (int a, int b) { return a > b; });           (Custom sort using lambda function)
