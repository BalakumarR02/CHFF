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

class Fin_variables
{
public:
    double A, Tf, Ab, Af, Ap, As, Auf, Atot, P, L1, W1, B, H, t, L, W, Tb, Tinf, N, S, qb;
    Fin_variables()
    {
        Tinf = 293;
        Tb = 373;
        H = 1000 * 0.001;
        L = 500 * 0.001;
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

void fin_1D_sample()
{
    ll n, i, m;
    Fin_variables f1;
    Flow_Parameters c1;
    // Insulated B.C
    cin >> n;

    vector<vector<double>> v(n + 1, vector<double>(5));
    double del_x = f1.L / n;
    double x_conv = f1.L1 / 1000;
    c1.solve_h(x_conv);
    cout << c1.hx << "\n";
    fo(i, 1, n)
    {

        double ap = 0, aw = 1 / del_x, ae = 1 / del_x, ab = 0, b = 0;
        double msq = (c1.h) * (f1.P) / ((c1.k) * (f1.Af));
        if (i == 1)
        {
            ab = 2 * aw;
            aw = 0;
        }
        else if (i == n)
        {
            ae = 0;
        }
        ap = aw + ae + (msq * del_x) + ab;
        b += msq * del_x * f1.Tinf + ab * f1.Tb;
        v[i][1] = ap, v[i][2] = ae, v[i][3] = aw, v[i][4] = b;
    }
    vector<double> phi = TDMA1D(n, v);

    // Temperature and heat flux at different points
    fo(i, 1, n)
    {
        cout << phi[i] << " " << c1.hx * f1.P * del_x * (phi[i] - f1.Tinf) << "  ";
        if (i % 10 == 0)
            cout << "\n";
    }
    cout << "\n";
}
int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    fin_1D_sample();

    return 0;
}

// sort(begin(v), end(v), [] (int a, int b) { return a > b; });           (Custom sort using lambda function)
